#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare phase5c to its frozen original baseline on common observations.

This diagnoses numerical/scientific discrepancies; it does not fit parameters
or declare field acceptance. Figure generation additionally uses matplotlib.
"""
import argparse
import csv
import json
from pathlib import Path
from seasonal_validation import metrics,save,sha,read_profiles


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase5c",type=Path);args=parser.parse_args();root=args.phase5c.resolve()
    corrected=root/"qualified_seasonal";original=root/"original_seasonal"
    a=json.loads((corrected/"registration.json").read_text());b=json.loads((original/"registration.json").read_text())
    if ({k:v["sha256"] for k,v in a["inputs"].items()}!={k:v["sha256"] for k,v in b["inputs"].items()} or
        a["observations"]["sha256"]!=b["observations"]["sha256"] or a["model_interval_jday"]!=b["model_interval_jday"]):
        raise ValueError("Input, observation or interval fingerprints differ")
    directories={"Original solver":original/"baseline","Corrected solver":corrected/"baseline",
        "Stable scalar mixing off":corrected/"scalar_mixing_off","Stable viscosity off":corrected/"viscosity_off",
        "Wind stress off":corrected/"no_wind_stress"}
    points={}
    for name,directory in directories.items():
        with (directory/"comparison.csv").open() as stream:
            points[name]={(float(row["jday"]),float(row["depth_m"])):row for row in csv.DictReader(stream)}
    common=set.intersection(*(set(table) for table in points.values()))
    for key in common:
        if len({table[key]["observed_c"] for table in points.values()})!=1:
            raise ValueError("Observed values differ across models")
    summary={name:{split:metrics([float(table[key]["error_c"]) for key in common
                   if split=="all" or table[key]["split"]==split])
                   for split in ("all","development","reserved_review")} for name,table in points.items()}
    review=root/"review";review.mkdir(parents=True,exist_ok=True)
    data={"common_points":len(common),"errors":summary,
        "numerical_qualification":{name:json.loads((directory/"metrics.json").read_text())["numerically_qualified"] for name,directory in directories.items()},
        "registrations_sha256":{"original":sha(original/"registration.json"),"corrected":sha(corrected/"registration.json")},
        "interpretation":"Verified operator corrections do not remove the seasonal warm bias. Off cases are diagnostic, not calibrated models. Autumn reservation is retrospective, not blind validation.",
        "next_scientific_work":["Audit wind projection/exposure against source conventions and mesh axes",
            "Qualify storage/conveyance momentum and sigma transport on sloping geometry",
            "Use independently verified solver derivatives to diagnose separate mixing and forcing sensitivities",
            "Carry structural discrepancy separately from conditional parameter uncertainty"]}
    save(review/"diagnosis.json",data)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    styles=(("#777777","--"),("#0072b2","-"),("#e69f00","-."),("#cc79a7",":"),("#009e73","--"))
    profiles={name:read_profiles(directory/"degray_profiles.csv") for name,directory in directories.items()}
    casts={c["jday"]:c for c in a["casts"]}
    fig,axes=plt.subplots(2,3,figsize=(13,9),layout="constrained",sharey=True)
    for ax,day in zip(axes.flat,(77.5,147.5,203.5,259.5,301.5,357.5)):
        cast=casts[day];ax.plot(cast["temperature_c"],cast["depth_m"],"k.",label="Observed")
        for (name,profile),(color,style) in zip(profiles.items(),styles):
            rows=profile[day];ax.plot([t for z,t in rows],[z for z,t in rows],color=color,linestyle=style,label=name)
        ax.set(xlabel="Temperature (°C)",title=f"Day {day:g} · {'reserved review' if day>=273.5 else 'development'}")
        ax.grid(alpha=.2)
    axes[0,0].invert_yaxis();axes[0,0].set_ylabel("Depth (m)");axes[1,0].set_ylabel("Depth (m)")
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=3,fontsize=9)
    fig.suptitle("DeGray 1980 · original and corrected operators\nWind and mixing interventions are diagnostic; substantial model discrepancy remains",fontsize=14)
    fig.savefig(review/"before_after_profiles.png",dpi=160);fig.savefig(review/"before_after_profiles.pdf");plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.5),layout="constrained")
    names=list(summary);values=[summary[name]["all"]["rmse_c"] for name in names]
    ax.barh(names,values,color=[style[0] for style in styles]);ax.invert_yaxis()
    for i,value in enumerate(values): ax.text(value+.07,i,f"{value:.2f}",va="center")
    ax.set(xlabel=f"Temperature RMSE (°C) · {len(common)} common supported observations",xlim=(0,max(values)*1.12),
           title="Numerical conservation passes; field accuracy remains unresolved")
    fig.savefig(review/"temperature_errors.png",dpi=160);plt.close(fig)
    print(json.dumps({name:values["all"] for name,values in summary.items()},indent=2))


if __name__=="__main__": main()
