#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build scientific figures and native Composer run artifacts from qualified CSV.

The GUI topology is W2's inferred chain, not surveyed reservoir bathymetry.
Plot dependencies: numpy and matplotlib; solver execution needs neither.
"""
import argparse
import csv
from contextlib import ExitStack
import datetime as dt
import json
from pathlib import Path
import uuid

from seasonal_validation import read_profiles, sample, metrics, save, sha


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def common_comparison(root, cases):
    points = {case: rows(root/case/"comparison.csv") for case in cases}
    lookup = {case: {(float(r["jday"]), float(r["depth_m"])): r for r in points[case]} for case in cases}
    common = set.intersection(*(set(table) for table in lookup.values()))
    summary = {}
    for case, table in lookup.items():
        summary[case] = {split: metrics([float(table[key]["error_c"]) for key in common
                                        if split == "all" or table[key]["split"] == split])
                         for split in ("all", "development", "reserved_review")}
    differences = {}
    baseline = lookup["baseline"]
    for case, table in lookup.items():
        differences[case] = metrics([float(table[key]["modeled_c"])-float(baseline[key]["modeled_c"]) for key in common])
    result = {"common_points": len(common), "common_depth_support": "Intersection of supported observed depths across every case; no extrapolation",
              "observation_errors": summary, "differences_from_baseline": differences,
              "qualification": "Changes on one finite grid/season do not establish asymptotic convergence or field acceptance."}
    save(root/"review/common_support.json", result)
    return result


def make_run(directory, caption, split_fields=False):
    """Lossless scalar/field conversion into SDK v1 recorded data contracts."""
    driver = json.loads((directory/"driver.json").read_text())
    geometry_rows = rows(directory/"degray_mesh.csv")
    columns, layers = driver["columns"], driver["layers"]
    polygons = {c: [] for c in range(columns)}
    nodes, coordinates, face_nodes, offsets = {}, [], [], [0]
    for row in geometry_rows:
        polygons[int(row["column"])].append((int(row["node_index"]), float(row["x"]), float(row["y"])))
    for c in range(columns):
        for _, x, y in sorted(polygons[c]):
            if (x, y) not in nodes:
                nodes[x, y] = len(coordinates); coordinates.append((x, y))
            face_nodes.append(nodes[x, y])
        offsets.append(len(face_nodes))
    fields = {}
    for row in rows(directory/"degray_field.csv"):
        day = float(row["jday"])
        fields.setdefault(day, {})[int(row["column"]), int(row["layer"])] = tuple(float(row[k]) for k in ("top_m", "bottom_m", "temperature_C"))
    days = sorted(fields)
    if not days or any(len(fields[day]) != columns*layers for day in days):
        raise ValueError("incomplete recorded DeGray field")
    diagnostic_fields = {
        "u_m_s": "m/s", "v_m_s": "m/s", "upward_relative_top_m_s": "m/s",
        "kz_top_m2_s": "m2/s", "kz_bottom_m2_s": "m2/s",
        "baroclinic_x_m_s2": "m/s2", "baroclinic_y_m_s2": "m/s2",
        "horizontal_advection_J": "J", "vertical_advection_J": "J",
        "vertical_diffusion_J": "J", "surface_J": "J", "sources_J": "J", "residual_J": "J"}
    diagnostics = {}
    if (directory/"degray_heat_transport.csv").exists():
        for r in rows(directory/"degray_heat_transport.csv"):
            diagnostics.setdefault(float(r["jday"]), {})[int(r["column"]),int(r["layer"])] = r
        if sorted(diagnostics) != days or any(len(diagnostics[d]) != columns*layers for d in days):
            raise ValueError("incomplete heat diagnostic field or mismatched clock")
        if all("km_top_m2_s" in row for field in diagnostics.values() for row in field.values()):
            diagnostic_fields.update(km_top_m2_s="m2/s",km_bottom_m2_s="m2/s")
    else:
        diagnostic_fields = {}
    bed = [fields[days[0]][c, layers-1][1] for c in range(columns)]
    node_bed = [[] for _ in coordinates]
    for c in range(columns):
        for n in face_nodes[offsets[c]:offsets[c+1]]:
            node_bed[n].append(bed[c])
    edges = set()
    for c in range(columns):
        polygon = face_nodes[offsets[c]:offsets[c+1]]
        edges.update(tuple(sorted((a,b))) for a,b in zip(polygon,polygon[1:]+polygon[:1]))
    topology = {"edge_nodes": [n for edge in sorted(edges) for n in edge], "schema": "hc.recorded-mesh/1.0", "node_x": [x for x, y in coordinates],
                "node_y": [y for x, y in coordinates], "node_z": [min(z) for z in node_bed],
                "face_nodes": face_nodes, "face_offsets": offsets,
                "reference": {"authority": "", "id": 0, "text": "", "units": 0,
                              "vertical_authority": "", "vertical_id": 0, "vertical_text": "", "vertical_units": 0}}
    save(directory/"gui_geometry.json", topology)
    # W2 day 1 is Jan 1, and .5 is noon. Convert once at this UI boundary.
    epoch = dt.date(1980, 1, 1).toordinal() + 1721424.5 - 1
    daily = {float(r["jday"]): r for r in rows(directory/"degray_daily.csv")}
    budgets = {float(r["jday"]): r for r in rows(directory/"degray_budgets.csv")}
    scalars = {"surface_temperature_C": "degC", "bottom_temperature_C": "degC", "thermocline_depth_m": "m",
               "surface_elevation_m": "m", "volume_m3": "m3", "heat_content_J": "J", "water_independent_residual_m3": "m3"}
    scalars.update({name: "m3" if name.startswith("budget_water") else "J" for name in next(iter(budgets.values())) if name != "jday"})
    momentum = {}
    if (directory/"degray_momentum.csv").exists():
        momentum = {float(row["jday"]): row for row in rows(directory/"degray_momentum.csv")}
        if sorted(momentum) != days: raise ValueError("momentum audit clock differs from field clock")
        scalars.update({name: "J" if name.endswith("_J") else "kg*m/s" if name.endswith("_kg_m_s") else "1"
            for name in next(iter(momentum.values())) if name != "jday"})
    interfaces_id = "degray_interfaces"
    def artifact(name):
        return "gui_"+name+".csv" if split_fields else "gui_fields.csv"
    with ExitStack() as stack:
        writers = {}
        names = ["temperature", interfaces_id, *diagnostic_fields, *scalars]
        shared = None
        for name in names:
            if split_fields or shared is None:
                shared = csv.writer(stack.enter_context((directory/artifact(name)).open("w")))
                shared.writerow(("time", "item", "index", "value"))
            writers[name] = shared
        for day in days:
            absolute = epoch + day
            for c in range(columns):
                for k in range(layers):
                    writers["temperature"].writerow((absolute, "temperature", c*layers+k, fields[day][c,k][2]))
                    for name in diagnostic_fields:
                        writers[name].writerow((absolute, name, c*layers+k, diagnostics[day][c,k][name]))
                for k in range(layers+1):
                    z = fields[day][c,k][0] if k < layers else fields[day][c,layers-1][1]
                    writers[interfaces_id].writerow((absolute, interfaces_id, c*(layers+1)+k, z))
            for name in scalars:
                record = daily[day] if name in daily[day] else budgets[day] if name in budgets[day] else momentum[day]
                writers[name].writerow((absolute, name, 0, record[name]))
    clock = {"count": len(days), "start": epoch+days[0], "end": epoch+days[-1], "units": "julian_day"}
    def entry(name, shape, dimensions, units):
        result = {"component": "degray", "item": name, "artifact": artifact(name), "format": "csv",
                "variable": name, "kind": "Float64", "shape": [len(days)]+shape,
                "dimensions": ["time"]+dimensions, "time": clock, "units": units}
        powers = {"m":[1,0,0,0,0], "m3":[3,0,0,0,0], "J":[2,1,-2,0,0],
                  "degC":[0,0,0,0,1], "m/s":[1,0,-1,0,0], "m2/s":[2,0,-1,0,0],
                  "m/s2":[1,0,-2,0,0], "kg*m/s":[1,1,-1,0,0], "1":[0,0,0,0,0]}
        if units in powers:
            result["quantity"]={"scale":1.,"offset":273.15 if units=="degC" else 0.,
                                "value_kind":1,"powers":powers[units]+[0]*5}
        return result
    vertical = entry(interfaces_id, [columns, layers+1], ["columns", "interfaces"], "m")
    # Nested vertical entry has no independent component or quantity.
    temperature = entry("temperature", [columns, layers], ["columns", "layers"], "degC")
    temperature.update(location="face", description="DeGray inferred W2 chain; cell-centre temperature. Plan footprints do not express hypsographic volumes.",
                       geometry={"artifact": "gui_geometry.json", "entity_axis": 1, "layer_axis": 2,
                                 "vertical_kind": 1, "vertical": vertical},
                       quantity={"offset": 273.15, "scale": 1., "value_kind": 1, "powers": [0.,0.,0.,0.,1.,0.,0.,0.,0.,0.]})
    results = [temperature] + [entry(name, [1], ["value"], units) for name, units in scalars.items()]
    for name, units in diagnostic_fields.items():
        diagnostic = entry(name,[columns,layers],["columns","layers"],units)
        diagnostic.update(location="face",geometry=temperature["geometry"],
            description=("Cumulative amount on moving sigma cells since initialization" if units == "J" else
                         "Last-step diagnostic snapshot; not a time average") + "; DeGray inferred W2 chain, field acceptance pending")
        results.append(diagnostic)
    completed = json.loads((directory/"provenance.json").read_text()) if (directory/"provenance.json").exists() else {}
    manifest = {"schema_version": "1.0", "run": {"id": str(uuid.uuid4()), "caption": caption,
                  "status": "completed", "wall_seconds": completed.get("wall_seconds", 0)},
                "components": [{"id": "degray", "caption": caption, "final_status": "Finished"}], "results": results,
                "provenance": {"source_csv": {p.name: sha(p) for p in directory.glob("degray_*.csv")},
                               "converter_sha256": sha(__file__), "geometry_scope": driver["geometry_scope"]}}
    save(directory/"gui_run.json", manifest)
    save(directory/"gui_review.composition.json", {"schema_version": "1.1", "metadata": {"id": "degray-recorded-review", "caption": caption},
         "components": [{"id": "degray", "caption": caption, "execution": {"mode": "open", "results": "gui_run.json"}}]})
    return len(days), columns, layers


def figures(root, cases, common):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 160})
    review = root/"review"
    baseline = read_profiles(root/"baseline/degray_profiles.csv")
    registration = json.loads((root/"registration.json").read_text())
    casts = {c["jday"]: c for c in registration["casts"]}
    fig, axes = plt.subplots(2, 3, figsize=(11, 8), sharey=True, constrained_layout=True)
    for ax, day in zip(axes.flat, (77.5,147.5,203.5,259.5,301.5,357.5)):
        c = casts[day]
        ax.plot(c["temperature_c"], c["depth_m"], "k.", markersize=4, label="Measured")
        for case in [c for c in ("baseline", "no_wind_stress", "outlet115") if c in cases]:
            profile = read_profiles(root/case/"degray_profiles.csv")
            data = profile.get(day, [])
            ax.plot([t for z,t in data], [z for z,t in data], label=case.replace("_", " "))
        ax.set_title(f"Day {day:g} · {'reserved review' if day >= 273.5 else 'development'}")
        ax.grid(alpha=.2); ax.set_xlabel("Temperature (°C)")
    axes[0,0].invert_yaxis(); axes[0,0].set_ylabel("Depth below surface (m)"); axes[1,0].set_ylabel("Depth below surface (m)")
    axes[0,0].legend(fontsize=8)
    fig.suptitle("DeGray dam profiles, 1980 · wind diagnostic is not a calibration")
    fig.savefig(review/"profiles.png"); plt.close(fig)

    points = rows(root/"baseline/comparison.csv")
    days = sorted(set(float(p["jday"]) for p in points))
    depths = sorted(set(float(p["depth_m"]) for p in points))
    lookup = {(float(p["jday"]),float(p["depth_m"])):float(p["error_c"]) for p in points}
    matrix = [[lookup.get((day,z),np.nan) for day in days] for z in depths]
    fig, ax = plt.subplots(figsize=(11,4.5), constrained_layout=True)
    limit = max(abs(e) for e in lookup.values())
    im = ax.pcolormesh(days, depths, matrix, cmap="RdBu_r", vmin=-limit, vmax=limit, shading="nearest")
    ax.axvline(273.5, color="black", linestyle="--", label="Reserved review begins")
    ax.invert_yaxis(); ax.set(xlabel="W2 Julian day, 1980", ylabel="Sampled depth (m)", title="Baseline minus measured temperature · unobserved depths stay blank")
    ax.legend(fontsize=8); fig.colorbar(im, ax=ax, label="Error (°C)")
    fig.savefig(review/"seasonal_errors.png"); plt.close(fig)

    fig, axes = plt.subplots(1,2, figsize=(11,4.5), constrained_layout=True)
    labels = list(cases)
    x = np.arange(len(labels)); width=.38
    for i, split in enumerate(("development", "reserved_review")):
        axes[0].bar(x+(i-.5)*width,[common["observation_errors"][c][split]["rmse_c"] for c in labels],width,label=split.replace("_"," "))
    axes[0].set_ylabel("RMSE on common observed depths (°C)"); axes[0].legend(fontsize=8)
    axes[1].bar(x,[common["differences_from_baseline"][c]["rmse_c"] for c in labels])
    axes[1].set_ylabel("RMS difference from baseline (°C)")
    for ax in axes:
        ax.set_xticks(x, labels, rotation=55, ha="right"); ax.grid(axis="y",alpha=.2)
    fig.suptitle("Numerical refinement and geometry/wind diagnostics · no parameter fitting")
    fig.savefig(review/"refinement.png"); plt.close(fig)

    daily = rows(root/"baseline/degray_daily.csv"); budgets = rows(root/"baseline/degray_budgets.csv")
    fig, axes = plt.subplots(3,1, figsize=(11,9), sharex=True, constrained_layout=True)
    for name in ("shortwave", "longwave_in", "back_radiation", "latent", "sensible"):
        axes[0].plot([float(r["jday"]) for r in budgets], [float(r["budget_heat_"+name])/1e15 for r in budgets], label=name.replace("_"," "))
    axes[0].set_ylabel("Cumulative surface energy (PJ)"); axes[0].legend(ncol=3,fontsize=8)
    axes[1].plot([float(r["jday"]) for r in daily], [float(r["source_heat_in_independent_J"])/1e15 for r in daily],label="Inflow heat")
    axes[1].plot([float(r["jday"]) for r in daily], [-float(r["source_heat_out_independent_J" if "source_heat_out_independent_J" in r else "source_heat_out_audit_J"])/1e15 for r in daily],label="Withdrawal heat")
    axes[1].set_ylabel("Independent source energy (PJ)"); axes[1].legend(fontsize=8)
    for name in ("water", "heat"):
        key="budget_"+name+"_residual"
        scale="volume_m3" if name=="water" else "heat_content_J"
        axes[2].plot([float(r["jday"]) for r in budgets],[abs(float(b[key]))/max(abs(float(d[scale])),1) for b,d in zip(budgets,daily)], label=name)
    axes[2].set_yscale("symlog",linthresh=1e-16); axes[2].set_ylabel("Relative inventory residual"); axes[2].legend(fontsize=8)
    axes[2].set_xlabel("W2 Julian day, 1980")
    for ax in axes: ax.grid(alpha=.2)
    fig.suptitle("DeGray water and heat audit · latent term is energy, not evaporation volume")
    fig.savefig(review/"budgets.png"); plt.close(fig)

    fig, axes = plt.subplots(1,2, figsize=(11,4.5), constrained_layout=True)
    for group, ax in zip((("baseline","dt150","dt75"),("layers15","baseline","layers60")), axes):
        for case in group:
            data=rows(root/case/"degray_daily.csv")
            ax.plot([float(r["jday"]) for r in data],[float(r["bottom_temperature_C"]) for r in data],label=case)
        ax.set(xlabel="W2 Julian day, 1980",ylabel="Deepest dam cell (°C)"); ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle("Timestep and vertical refinement · cell-centre depths differ by grid")
    fig.savefig(review/"deep_temperature.png"); plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seasonal",type=Path)
    parser.add_argument("--outlet-case", type=Path)
    parser.add_argument("--split-fields", action="store_true", help="Separate variables so field access does not scan unrelated annual diagnostics")
    parser.add_argument("--caption-prefix", default="DeGray 1980")
    args=parser.parse_args(); root=args.seasonal.resolve()
    qualification=json.loads((root/"qualification.json").read_text())
    cases=list(json.loads((root/"registration.json").read_text())["cases"])
    if args.outlet_case:
        target = root/"outlet115"
        if not target.exists(): target.symlink_to(args.outlet_case.resolve(), target_is_directory=True)
        cases.append("outlet115")
    common=common_comparison(root,cases)
    save(root/"review/status.json", {"experiments": cases,
         "numerical_qualification": {case: json.loads((root/case/"metrics.json").read_text())["numerically_qualified"] for case in cases},
         "field_acceptance": qualification["field_acceptance"],
         "geometry_scope": "Inferred W2 geometry, not surveyed 3D bathymetry. Footprint prisms do not represent hypsographic volumes."})
    figures(root,cases,common)
    for case in cases:
        make_run(root/case,args.caption_prefix+" · "+case.replace("_"," ")+" · inferred W2 geometry · field review pending",args.split_fields)
    print(json.dumps(common,indent=2))


if __name__=="__main__":
    main()
