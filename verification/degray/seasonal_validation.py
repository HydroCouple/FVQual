#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Register and reproduce seasonal DeGray validation. Standard library only.

No fitting is performed. September–December casts are reserved from future
parameter selection; they are NOT blind because historical work used the data.
Depth interpolation is within cell centres only, with no time/depth extrapolation.
"""
import argparse
import bisect
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import time

START, END, HOLDOUT = 64.5, 358.7, 273.5
CASES = {
    "baseline": (30, 300, 1, ""),
    "dt150": (30, 150, 1, ""),
    "dt75": (30, 75, 1, ""),
    "layers15": (15, 300, 1, ""),
    "layers60": (60, 300, 1, ""),
    "lanes2": (30, 300, 2, ""),
    "no_wind_stress": (30, 300, 1, "nowind"),
    "reduced_shear_mixing": (30, 300, 1, "nomix"),
    "no_sources": (30, 300, 1, "nosources"),
    "scalar_mixing_off": (30, 300, 1, "scalar_nomix"),
    "viscosity_off": (30, 300, 1, "viscosity_nomix"),
    "outlet115": (30, 300, 1, "selective115"),
}
INPUTS = ("bth1.csv", "met.npt", "qin_br1.npt", "tin_br1.npt", "qot_br1.npt", "vpr.npt", "wsc.npt")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def read_casts(path):
    """Strict CNO/DAYO/NRO parser; retain the supplied station and clock."""
    lines = Path(path).read_text().splitlines()
    casts, i = [], 0
    while i < len(lines):
        if "CNO" not in lines[i] or "DAYO" not in lines[i]:
            i += 1
            continue
        i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1
        cno, day, count = map(float, lines[i].split())
        if count != int(count) or count < 1 or not math.isfinite(day):
            raise ValueError("invalid cast header")
        i += 1
        while i < len(lines) and "ZO" not in lines[i]:
            i += 1
        if i == len(lines):
            raise ValueError("missing cast depth header")
        i += 1
        values = []
        while i < len(lines) and len(values) < 2 * count:
            if "CNO" in lines[i]:
                raise ValueError("truncated cast")
            values.extend(map(float, lines[i].split()))
            i += 1
        if len(values) != 2 * count or not all(map(math.isfinite, values)):
            raise ValueError("invalid cast values")
        if cno != 22:
            continue
        depth, temperature = values[::2], values[1::2]
        if depth[0] < 0 or any(b <= a for a, b in zip(depth, depth[1:])):
            raise ValueError("cast depths must be nonnegative and strictly increasing")
        if any(c["jday"] == day for c in casts):
            raise ValueError("duplicate temperature cast")
        casts.append({"station": "DeGray dam (only supplied station)", "jday": day,
                      "depth_m": depth, "temperature_c": temperature})
    if not casts:
        raise ValueError("no CNO 22 temperature casts")
    return sorted(casts, key=lambda c: c["jday"])


def read_profiles(path):
    profiles = {}
    for row in csv.DictReader(Path(path).open()):
        day, depth, value = (float(row[k]) for k in ("jday", "depth_m", "temperature_C"))
        if not all(map(math.isfinite, (day, depth, value))):
            raise ValueError("nonfinite model profile")
        profiles.setdefault(day, []).append((depth, value))
    for rows in profiles.values():
        rows.sort()
        if any(b[0] <= a[0] for a, b in zip(rows, rows[1:])):
            raise ValueError("duplicate/nonmonotonic model depths")
    return profiles


def interpolate(rows, depth):
    if depth < rows[0][0] or depth > rows[-1][0]:
        return None
    i = bisect.bisect_left([r[0] for r in rows], depth)
    if i == 0:
        return rows[0][1]
    if i == len(rows):
        return rows[-1][1]
    lo, hi = rows[i-1], rows[i]
    return lo[1] + (hi[1]-lo[1]) * (depth-lo[0]) / (hi[0]-lo[0])


def sample(profiles, day, depth):
    """Exact time when present; otherwise bracket within <=1 daily interval."""
    days = sorted(profiles)
    index = bisect.bisect_left(days, day)
    if index < len(days) and abs(days[index] - day) < 1e-9:
        return interpolate(profiles[days[index]], depth)
    if index == 0 or index == len(days) or days[index] - days[index-1] > 1.000001:
        return None
    a, b = sample(profiles, days[index-1], depth), sample(profiles, days[index], depth)
    if a is None or b is None:
        return None
    weight = (day-days[index-1]) / (days[index]-days[index-1])
    return a * (1-weight) + b * weight


def metrics(errors):
    return {"n": len(errors), "bias_c": sum(errors)/len(errors),
            "rmse_c": math.sqrt(sum(e*e for e in errors)/len(errors)),
            "max_abs_c": max(map(abs, errors))} if errors else {"n": 0}


def compare(profiles, casts):
    points, excluded = [], []
    for cast in casts:
        day = cast["jday"]
        if day < START or day > END:
            excluded.append({"jday": day, "reason": "outside model interval", "n": len(cast["depth_m"])})
            continue
        split = "reserved_review" if day >= HOLDOUT else "development"
        for z, observed in zip(cast["depth_m"], cast["temperature_c"]):
            modeled = sample(profiles, day, z)
            if modeled is None:
                excluded.append({"jday": day, "depth_m": z, "reason": "outside model time/depth support"})
                continue
            points.append({"jday": day, "split": split, "depth_m": z,
                           "observed_c": observed, "modeled_c": modeled, "error_c": modeled-observed})
    summary = {split: metrics([p["error_c"] for p in points if split == "all" or p["split"] == split])
               for split in ("all", "development", "reserved_review")}
    summary["casts"] = [{"jday": day, "split": "reserved_review" if day >= HOLDOUT else "development",
                         **metrics([p["error_c"] for p in points if p["jday"] == day])}
                        for day in sorted(set(p["jday"] for p in points))]
    summary["bands"] = [{"depth_m": [lo, hi], **metrics([p["error_c"] for p in points if lo <= p["depth_m"] < hi])}
                        for lo, hi in zip((0, 5, 10, 20, 30), (5, 10, 20, 30, 60))]
    summary["excluded"] = excluded
    return summary, points


def budget_metrics(path):
    rows = list(csv.DictReader(Path(path).open()))
    daily = list(csv.DictReader((Path(path).parent/"degray_daily.csv").open()))
    result = {}
    for name in ("water", "heat"):
        residual = [abs(float(r[f"budget_{name}_residual"])) for r in rows]
        # Relative to actual inventory as an independently computed diagnostic.
        scales = [abs(float(r["volume_m3" if name == "water" else "heat_content_J"])) for r in daily]
        if len(rows) != len(daily) or any(float(a["jday"]) != float(b["jday"]) for a,b in zip(rows,daily)):
            raise ValueError("budget and inventory clocks differ")
        result[name] = {"max_abs_residual": max(residual), "max_relative_inventory_residual": max(r/max(scale,1) for r,scale in zip(residual,scales))}
    result["independent_water"] = {"max_abs_residual": max(abs(float(r["water_independent_residual_m3"])) for r in daily)}
    result["source_heat_in"] = {"max_abs_ledger_difference_J": max(abs(float(b["budget_heat_source_in"])-float(d["source_heat_in_independent_J"])) for b,d in zip(rows,daily))}
    audit_key = "source_heat_out_independent_J" if "source_heat_out_independent_J" in daily[0] else "source_heat_out_audit_J"
    result["source_heat_out"] = {"max_abs_ledger_difference_J": max(abs(float(b["budget_heat_source_out"])-float(d[audit_key])) for b,d in zip(rows,daily)),
                                  "independent_placement": json.loads((Path(path).parent/"driver.json").read_text()).get("outlet_placement", "uniform") == "uniform"}
    return result


def protocol(driver, inputs, observed, output, selected, end, heat_audit=False):
    casts = read_casts(observed)
    registration = {
        "schema": "fvqual.degray-validation/1.0", "year": 1980,
        "clock": "W2 Julian day, day 1 = January 1; noon .5; no timezone conversion",
        "station_count": 1, "casts": casts, "model_interval_jday": [START, end],
        "development_casts": [c["jday"] for c in casts if START <= c["jday"] < HOLDOUT and c["jday"] <= end],
        "reserved_review_casts": [c["jday"] for c in casts if HOLDOUT <= c["jday"] <= end],
        "excluded_preseason_casts": [c["jday"] for c in casts if c["jday"] < START],
        "split_disclosure": "Retrospective reservation for future tuning, not blind; all observations were used in earlier work.",
        "geometry": "Laterally averaged W2 storage curves; lanes2 infers lateral depth from widths, not surveyed 3D bathymetry. Treat as geometry sensitivity, not formal horizontal convergence.",
        "sampling": "Deepest dam lane; piecewise-linear depth between cell centres; no extrapolation. Exact daily times and bracketed <=1 day only.",
        "assumptions": {"bottom_elevation_m": 66.35, "latitude_degrees": 34.2,
                        "wind_height_m": 10, "evaporation_height_m": 2, "attenuation_per_m": .45,
                        "cloud_conversion": "tenths / 10", "solar": "computed clear-sky radiation with cloud correction",
                        "withdrawal": "uniform baseline; outlet115 uses workbook C158=115 m, C153=POINT, C168=3.14 radians with historical variable-flow files",
                        "no_wind_stress": "diagnostic only; weather heat exchange retained; not a fitted model"},
        "acceptance": {"numerical": "Complete requested interval, finite output, no barotropic floor violation, relative water/heat inventory residual <=1e-8",
                       "scientific": "Review seasonal/depth errors and refinement. No error threshold supplied; this workflow cannot declare field acceptance or choose tuned parameters."},
        "driver": str(driver), "driver_sha256": sha(driver), "script_sha256": sha(__file__), "core_sha256": sha(driver.parent.parent / "src/libfvqual-core.dylib"),
        "inputs": {name: {"path": str(inputs/name), "sha256": sha(inputs/name)} for name in INPUTS},
        "observations": {"path": str(observed), "sha256": sha(observed)},
        "heat_transport_audit": heat_audit,
        "heat_audit_scope": "Cumulative heat on moving sigma cells; last-step velocity/diffusivity/pressure snapshots, not time averages. Legacy reduced shear mixing changes both coefficients. scalar_mixing_off only zeros stable Richardson diffusivity; viscosity_off only zeros stable Richardson viscosity. Both retain unstable convective adjustment. No sources retains wind and weather heat exchange.",
        "cases": {key: {"layers": CASES[key][0], "max_dt_s": CASES[key][1], "max_lanes": CASES[key][2], "diagnostic": CASES[key][3]} for key in selected},
    }
    save(output / "registration.json", registration)
    return registration


def run_case(key, registration, driver, inputs, output, end):
    directory = output/key
    directory.mkdir(parents=True, exist_ok=True)
    layers, dt, lanes, diagnostic = CASES[key]
    env = dict(os.environ, FVQUAL_LAYERS=str(layers), FVQUAL_MAX_DT=str(dt), FVQUAL_LANES=str(lanes), FVQUAL_MIN_CELL="200" if lanes > 1 else "0")
    # Prevent unrelated inherited diagnostic settings from changing the experiment.
    env.pop("FVQUAL_TRACE_COLUMN", None)
    env.pop("FVQUAL_OUTLET_ELEVATION", None)
    env["FVQUAL_HEAT_AUDIT"] = "1" if registration.get("heat_transport_audit") else "0"
    if diagnostic == "selective115": env["FVQUAL_OUTLET_ELEVATION"] = "115"
    command = [str(driver), str(inputs), str(directory), str(end), "richardson", "" if diagnostic == "selective115" else diagnostic]
    manifest = {"case": key, "command": command, "environment": {k: env[k] for k in ("FVQUAL_LAYERS", "FVQUAL_MAX_DT", "FVQUAL_LANES", "FVQUAL_MIN_CELL", "FVQUAL_OUTLET_ELEVATION", "FVQUAL_HEAT_AUDIT") if k in env},
                "driver_sha256": registration["driver_sha256"], "core_sha256": registration["core_sha256"],
                "registration_sha256": sha(output/"registration.json"), "platform": platform.platform(), "status": "running"}
    save(directory/"provenance.json", manifest)
    print(f"starting {key}", flush=True)
    started = time.monotonic()
    with (directory/"driver.log").open("w") as log:
        completed = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
    manifest.update(wall_seconds=time.monotonic()-started, returncode=completed.returncode)
    if completed.returncode:
        manifest["status"] = "failed"
        save(directory/"provenance.json", manifest)
        raise RuntimeError(f"{key} failed; see {directory/'driver.log'}")
    if sha(driver) != registration["driver_sha256"] or sha(driver.parent.parent/"src/libfvqual-core.dylib") != registration["core_sha256"]:
        raise RuntimeError("driver/core changed during experiment")
    details = json.loads((directory/"driver.json").read_text())
    for artifact in directory.glob("degray_*.csv"):
        with artifact.open() as stream:
            for row in csv.DictReader(stream):
                if not all(math.isfinite(float(value)) for value in row.values()):
                    raise ValueError(f"nonfinite output in {artifact}")
    summary, points = compare(read_profiles(directory/"degray_profiles.csv"), registration["casts"])
    summary["budgets"] = budget_metrics(directory/"degray_budgets.csv")
    if (directory/"degray_momentum.csv").exists():
        with (directory/"degray_momentum.csv").open() as stream:
            momentum = list(csv.DictReader(stream))
        summary["momentum"] = {key: max(abs(float(r[key])) for r in momentum)
            for key in ("max_relative_impulse_residual", "max_relative_energy_residual")}
        summary["momentum"]["scope"] = "Implicit vertical operator only; mode correction measured separately"
    summary["driver"] = details
    summary["numerically_qualified"] = (abs(details["end_jday"]-end) < 1e-9 and not details["barotropic_floor_binding"] and
                                        all(summary["budgets"][name]["max_relative_inventory_residual"] <= 1e-8 for name in ("water", "heat")) and
                                        (not registration.get("heat_transport_audit") or details["maximum_relative_cell_heat_residual"] <= 1e-8) and
                                        all(summary.get("momentum",{}).get(key,0) <= 1e-8 for key in
                                            ("max_relative_impulse_residual", "max_relative_energy_residual")))
    save(directory/"metrics.json", summary)
    with (directory/"comparison.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(points[0]) if points else ["jday"])
        writer.writeheader(); writer.writerows(points)
    manifest.update(status="completed", outputs={p.name: sha(p) for p in directory.glob("*.csv")}, metrics=summary["all"])
    save(directory/"provenance.json", manifest)
    print(f"completed {key}: {summary['all']}", flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=[key for key in CASES if key not in ("outlet115", "reduced_shear_mixing", "no_sources", "scalar_mixing_off", "viscosity_off")])
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--end-day", type=float, default=END)
    parser.add_argument("--heat-audit", action="store_true")
    args = parser.parse_args()
    for name in ("driver", "inputs", "observed", "output"):
        setattr(args, name, getattr(args, name).resolve())
    registration = protocol(args.driver, args.inputs, args.observed, args.output, args.cases, args.end_day, args.heat_audit)
    summaries = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        tasks = {pool.submit(run_case, key, registration, args.driver, args.inputs, args.output, args.end_day): key for key in args.cases}
        for future in concurrent.futures.as_completed(tasks):
            summaries[tasks[future]] = future.result()
    save(args.output/"qualification.json", {"cases": summaries, "field_acceptance": "pending scientific review; no fitting performed"})


if __name__ == "__main__":
    main()
