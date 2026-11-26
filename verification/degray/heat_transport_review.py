#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independently close saved cell heat ledgers and build diagnostic review.

Amounts follow moving sigma cells. A Celsius heat reference is used consistently;
the volume term is explicit when attributing a temperature change. No fitting.
"""
import argparse
import csv
import json
import math
from pathlib import Path

from review_artifacts import rows, make_run
from seasonal_validation import save, sha

CAPACITY = 1000.0 * 4186.0
TERMS = ("horizontal_advection_J", "vertical_advection_J", "sources_J",
         "horizontal_diffusion_J", "vertical_diffusion_J", "surface_J", "other_J")
LABELS = {"baseline": "Normal wind", "no_wind_stress": "Zero wind stress",
          "reduced_shear_mixing": "Reduced shear mixing", "no_sources": "Closed basin with wind"}


def evaluate(directory):
    driver = json.loads((directory/"driver.json").read_text())
    data = [{key: float(value) for key, value in r.items()}
            for r in rows(directory/"degray_heat_transport.csv")]
    frames = {}
    for r in data:
        if not all(map(math.isfinite, r.values())):
            raise ValueError("nonfinite heat transport record")
        key = int(r["column"]),int(r["layer"])
        frame = frames.setdefault(r["jday"], {})
        if key in frame:
            raise ValueError("duplicate cell heat record")
        frame[key] = r
    initial = frames[min(frames)]
    scale_residual = 0.0
    interface_errors = {"advection": 0.0, "diffusion": 0.0}
    global_internal = {key: 0.0 for key in ("horizontal_advection_J","vertical_advection_J",
                                            "horizontal_diffusion_J","vertical_diffusion_J")}
    series = []
    for day, frame in sorted(frames.items()):
        if len(frame) != driver["columns"]*driver["layers"]:
            raise ValueError("incomplete heat frame")
        for key,r in frame.items():
            expected = r["initial_heat_J"] + math.fsum(r[t] for t in TERMS)
            scale = max(1,abs(r["heat_J"]),abs(r["initial_heat_J"]))
            scale_residual = max(scale_residual,abs(r["heat_J"]-expected)/scale)
            for process in interface_errors:
                amount = r[f"vertical_{process}_J"]
                face_amount = r[f"upward_{process}_bottom_J"]-r[f"upward_{process}_top_J"]
                interface_errors[process] = max(interface_errors[process],abs(amount-face_amount)/scale)
        inventory = math.fsum(r["heat_J"] for r in frame.values())
        for term in global_internal:
            global_internal[term] = max(global_internal[term],
                abs(math.fsum(r[term] for r in frame.values()))/max(1,abs(inventory)))
        bottom_key = driver["dam_column"],driver["layers"]-1
        bottom, first = frame[bottom_key],initial[bottom_key]
        current_capacity = CAPACITY*bottom["volume_m3"]
        initial_t = first["heat_J"]/(CAPACITY*first["volume_m3"])
        geometry_c = initial_t * (first["volume_m3"]/bottom["volume_m3"]-1)
        contribution = {t: bottom[t]/current_capacity for t in TERMS}
        series.append({"jday":day,"bottom_temperature_C":bottom["heat_J"]/current_capacity,
                       "initial_temperature_C":initial_t,"volume_contribution_C":geometry_c,
                       "advection_C":contribution["horizontal_advection_J"]+contribution["vertical_advection_J"],
                       "advection_and_sources_C":contribution["horizontal_advection_J"]+contribution["vertical_advection_J"]+contribution["sources_J"],
                       "diffusion_C":contribution["horizontal_diffusion_J"]+contribution["vertical_diffusion_J"],
                       "sources_C":contribution["sources_J"],"surface_C":contribution["surface_J"],
                       "other_C":contribution["other_J"],
                       "horizontal_advection_C":contribution["horizontal_advection_J"],
                       "vertical_advection_C":contribution["vertical_advection_J"]})
    for r in series:
        closed = r["initial_temperature_C"]+sum(r[k] for k in (
            "volume_contribution_C","advection_C","diffusion_C","sources_C","surface_C","other_C"))
        if abs(closed-r["bottom_temperature_C"]) > 1e-6:
            raise ValueError("temperature attribution does not close")
    with (directory/"dam_bottom_attribution.csv").open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=series[0]);writer.writeheader();writer.writerows(series)
    last = frames[max(frames)]
    dam = [last[driver["dam_column"],k] for k in range(driver["layers"])]
    review = {"maximum_relative_cell_residual":scale_residual,
              "maximum_relative_interface_divergence_residual":interface_errors,
              "maximum_relative_global_internal_residual":global_internal,
              "dam_bottom_final":series[-1],
              "last_snapshot_max_speed_m_s":max(math.hypot(r["u_m_s"],r["v_m_s"]) for r in last.values()),
              "last_snapshot_max_baroclinic_acceleration_m_s2":max(math.hypot(r["baroclinic_x_m_s2"],r["baroclinic_y_m_s2"]) for r in last.values()),
              "last_snapshot_dam_max_kz_m2_s":max(r["kz_bottom_m2_s"] for r in dam),
              "source_sha256":sha(directory/"degray_heat_transport.csv")}
    review["audit_passed"] = max([scale_residual]+list(interface_errors.values())+list(global_internal.values()))<=1e-8
    save(directory/"heat_audit_review.json",review)
    return series,dam,review


def figures(root, cases, evaluated):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                         "axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":160})
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for case in cases:
        series,dam,_=evaluated[case]
        x=[r["jday"] for r in series]
        axes[0,0].plot(x,[r["bottom_temperature_C"] for r in series],label=LABELS.get(case,case))
        z=[(dam[0]["top_m"]-.5*(r["top_m"]+r["bottom_m"])) for r in dam]
        axes[0,1].plot([r["heat_J"]/(CAPACITY*r["volume_m3"]) for r in dam],z,label=LABELS.get(case,case))
    baseline=evaluated[cases[0]][0]; x=[r["jday"] for r in baseline]
    for key,label in (("advection_and_sources_C","Advection + sources"),("diffusion_C","Diffusion"),
                      ("surface_C","Surface exchange"),("volume_contribution_C","Changing volume")):
        axes[1,0].plot(x,[r[key] for r in baseline],label=label)
    for case in cases:
        _,dam,_=evaluated[case]
        z=[dam[0]["top_m"]-.5*(r["top_m"]+r["bottom_m"]) for r in dam]
        axes[1,1].plot([math.hypot(r["u_m_s"],r["v_m_s"]) for r in dam],z,label=LABELS.get(case,case))
    axes[0,0].set(title="Dam bottom temperature",xlabel="1980 W2 Julian day",ylabel="Temperature (°C)")
    axes[0,1].set(title="Final dam profile",xlabel="Temperature (°C)",ylabel="Depth below surface (m)")
    axes[1,0].set(title="Normal wind: bottom-cell change attribution",xlabel="1980 W2 Julian day",ylabel="Contribution to change (°C)")
    axes[1,1].set(title="Final dam speed · instantaneous snapshot",xlabel="Horizontal speed (m/s)",ylabel="Depth below surface (m)")
    for ax in (axes[0,1],axes[1,1]):ax.invert_yaxis()
    for ax in axes.flat:ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8);axes[1,0].legend(fontsize=8)
    fig.suptitle("DeGray spring heat-transport diagnosis · no parameter fitting")
    fig.savefig(root/"review/spring_diagnosis.png");plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("--closed-basin",type=Path)
    args=parser.parse_args();root=args.root.resolve();(root/"review").mkdir(exist_ok=True)
    registration=json.loads((root/"registration.json").read_text())
    cases=list(registration["cases"])
    qualification=json.loads((root/"qualification.json").read_text())
    evaluated={case:evaluate(root/case) for case in cases}
    directories={case:root/case for case in cases}
    external=None
    if args.closed_basin:
        external=args.closed_basin.resolve()
        other=json.loads((external/"registration.json").read_text())
        if any(registration[k] != other[k] for k in ("driver_sha256","core_sha256","inputs","model_interval_jday")):
            raise ValueError("closed-basin experiment inputs/software/interval differ")
        case="no_sources";cases.append(case);directories[case]=external/case
        evaluated[case]=evaluate(directories[case])
        qualification["cases"][case]=json.loads((external/"qualification.json").read_text())["cases"][case]
    for case in cases:
        if not qualification["cases"][case]["numerically_qualified"] or not evaluated[case][2]["audit_passed"]:
            raise ValueError(f"{case} failed numerical/audit qualification")
        make_run(directories[case],f"DeGray spring · {LABELS.get(case,case)} · diagnostic review")
    figures(root,cases,evaluated)
    save(root/"review/status.json",{"cases":{case:evaluated[case][2] for case in cases},
         "numerical_and_audit_passed":True,"field_acceptance":"pending",
         "closed_basin_experiment":str(external) if external else None,
         "scope":"Spring diagnostic only. No seasonal replacement, calibration, or pressure/conveyance causality claim.",
         "interpretation":"Cumulative amounts follow moving sigma cells. Bottom temperature contributions include changing volume. Advection includes physical circulation and numerical mixing; separating those requires refinement/operators. Velocity/Kz/acceleration are last-step snapshots.",
         "script_sha256":sha(__file__)})


if __name__=="__main__":main()
