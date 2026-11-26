#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Convert closed-basin diagnostics into scientific figures and Composer data.

Conversion uses the standard library; --figures additionally needs matplotlib
and numpy. Grid comparisons use common physical positions, not layer indices.
"""
import argparse
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import uuid


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_run(case, root, source_sha):
    directory=root/case["id"]
    directory.mkdir(parents=True, exist_ok=True)
    geometry=dict(case["geometry"], schema="hc.recorded-mesh/1.0")
    offsets, nodes=geometry["face_offsets"], geometry["face_nodes"]
    edges=set()
    for a,b in zip(offsets,offsets[1:]):
        ring=nodes[a:b]
        edges.update(tuple(sorted(pair)) for pair in zip(ring,ring[1:]+ring[:1]))
    geometry["edge_nodes"]=[n for pair in sorted(edges) for n in pair]
    geometry["reference"]={"authority":"", "id":0, "text":"", "units":0,
        "vertical_authority":"", "vertical_id":0, "vertical_text":"", "vertical_units":0}
    save(directory/"geometry.json",geometry)
    epoch=datetime.date(2026,10,13).toordinal()+1721424.5
    frames=case["frames"]; n,k=case["columns"],case["layers"]
    if not frames or n<1 or k<1 or len(offsets)!=n+1:
        raise ValueError("Empty or inconsistent case geometry")
    times=[f["time_s"] for f in frames]
    if times[0]!=0 or any(not math.isfinite(t) for t in times) or any(a>=b for a,b in zip(times,times[1:])):
        raise ValueError("Recorded times must increase strictly from zero")
    fields={"temperature":"degC", "tracer":"kg/m3", "velocity_x":"m/s",
            "velocity_y":"m/s", "pressure_x":"m/s2"}
    scalars={key:("J" if key.endswith("_J") else "kg*m/s" if key.endswith("_kg_m_s") else
                 "m/s" if key.endswith("_m_s") or key=="max_speed" else
                 "kg" if key=="tracer_inventory_kg" else "degC*m3" if key=="heat_inventory_C_m3" else "1")
             for key in frames[0]["scalars"]}
    with (directory/"fields.csv").open("w") as stream:
        writer=csv.writer(stream);writer.writerow(("time","item","index","value"))
        for frame in frames:
            if any(len(frame[item])!=n*(k+1) for item in ("interfaces","viscosity","diffusivity")):
                raise ValueError("Wrong interface shape")
            if set(frame["scalars"])!=set(scalars) or any(not math.isfinite(v) for v in frame["scalars"].values()):
                raise ValueError("Inconsistent or nonfinite scalar record")
            if any(not math.isfinite(v) for item in ("interfaces","viscosity","diffusivity") for v in frame[item]):
                raise ValueError("Nonfinite interface data")
            t=epoch+frame["time_s"]/86400
            for item in fields:
                if len(frame[item])!=n*k: raise ValueError("Wrong field shape")
                for i,value in enumerate(frame[item]):
                    if not math.isfinite(value): raise ValueError("Nonfinite field")
                    writer.writerow((t,item,i,value))
            for i,z in enumerate(frame["interfaces"]): writer.writerow((t,"interfaces",i,z))
            for item in ("viscosity","diffusivity"):
                for c in range(n):
                    for layer in range(k): writer.writerow((t,item,c*k+layer,frame[item][c*(k+1)+layer]))
            for item,value in frame["scalars"].items(): writer.writerow((t,item,0,value))
    clock={"count":len(frames),"start":epoch,"end":epoch+frames[-1]["time_s"]/86400,"units":"julian_day"}
    powers={"degC":[0,0,0,0,1],"kg/m3":[-3,1,0,0,0],"m/s":[1,0,-1,0,0],
            "m/s2":[1,0,-2,0,0],"m2/s":[2,0,-1,0,0],"J":[2,1,-2,0,0],
            "kg*m/s":[1,1,-1,0,0],"kg":[0,1,0,0,0],"1":[0,0,0,0,0]}
    def entry(item,shape,dimensions,units):
        result={"component":"reservoir","item":item,"artifact":"fields.csv","format":"csv","variable":item,
            "kind":"Float64","shape":[len(frames)]+shape,"dimensions":["time"]+dimensions,"time":clock,"units":units}
        if units in powers: result["quantity"]={"scale":1,"offset":273.15 if units=="degC" else 0,
            "value_kind":1,"powers":powers[units]+[0]*5}
        return result
    vertical=entry("interfaces",[n,k+1],["columns","interfaces"],"m")
    results=[]
    fields.update(viscosity="m2/s",diffusivity="m2/s")
    for item,units in fields.items():
        result=entry(item,[n,k],["columns","layers"],units)
        result.update(location="face",geometry={"artifact":"geometry.json","entity_axis":1,"layer_axis":2,
            "vertical_kind":1,"vertical":vertical},description="Synthetic closed basin; prescribed storage and rectangular conveyance. Footprints do not express tapered volumes. Mixing is the upper-interface last-step coefficient.")
        results.append(result)
    results.extend(entry(item,[1],["value"],units) for item,units in scalars.items())
    caption="Closed basin — "+case["id"].replace("_"," ")
    manifest={"schema_version":"1.0","run":{"id":str(uuid.uuid4()),"caption":caption,"status":"completed"},
        "components":[{"id":"reservoir","caption":caption,"final_status":"Finished"}],"results":results,
        "provenance":{"source_sha256":source_sha,"converter_sha256":sha(Path(__file__)),
            "fields_sha256":sha(directory/"fields.csv"),"settings":{key:value for key,value in case.items() if key not in ("geometry","frames")},
            "qualification":"Numerical closed-basin diagnostic; not field validation or full mechanical energy closure"}}
    save(directory/"run.json",manifest)
    return directory/"run.json"


def make_composition(root,manifests):
    baseline=manifests["baseline"]
    composition={"schema_version":"1.1","metadata":{"id":"phase5c-basin-review","caption":"Phase 5c — controlled circulation review"},
        "components":[{"id":"reservoir","caption":"Synthetic closed basin","execution":{"mode":"open","results":str(baseline)}}]}
    save(root/"basin_review.composition.json",composition)
    save(root/"basin_native_review.composition.json",composition)
    layer={"type":"recorded-result","manifest":str(baseline),"component":"reservoir","item":"temperature",
        "name":"Closed basin baseline — temperature","units":"degC","visible":True,"opacity":1,"verticalLayer":1,
        "scene":{"shownIn3D":True,"terrainEnabled":True,"zMode":2,"offset":0}}
    run=json.loads(baseline.read_text())["run"]
    layer.update(runId=run["id"],runTitle=run["caption"],runStatus=run["status"])
    save(root/"basin_review.composition.composer.json",{"sidecar_version":1,"layers":[layer],
        "results":{"version":1,"runs":[str(p) for p in manifests.values()],"alignment":0}})
    save(root/"basin_native_review.composition.composer.json",{"sidecar_version":1,"layers":[layer],
        "results":{"version":1,"runs":[str(baseline)],"alignment":0}})


def figures(data,root):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cases={c["id"]:c for c in data["cases"]}
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout="constrained")
    palette=("#0072b2","#e69f00","#009e73","#000000","#cc79a7")
    styles=("-","--","-.",":","-")
    for name,color,style in zip(("baseline","scalar_mixing_off","viscosity_off","no_wind","narrow_conveyance"),palette,styles):
        c=cases[name];frames=c["frames"]
        t=[f["time_s"]/3600 for f in frames]
        axes[0,0].plot(t,[f["scalars"]["max_speed"] for f in frames],label=name.replace("_"," "),color=color,linestyle=style)
        f=frames[-1];n,k=c["columns"],c["layers"];i=n//2
        z=np.array(f["interfaces"]).reshape(n,k+1)[i]
        axes[0,1].plot(np.array(f["temperature"]).reshape(n,k)[i],-.5*(z[:-1]+z[1:]),label=name.replace("_"," "),color=color,linestyle=style)
    axes[0,0].set(xlabel="Time (hours)",ylabel="Maximum speed (m/s)",title="Independent mixing interventions")
    axes[0,0].legend(fontsize=8)
    axes[0,1].set(xlabel="Temperature (°C)",ylabel="Depth (m)",title="Mid-basin temperature at two hours");axes[0,1].invert_yaxis()
    for profile in ("linear","exponential"):
        values=[v for v in data["pressure"] if v["profile"]==profile]
        axes[1,0].loglog([v["layers"] for v in values],[max(v["max_acceleration_m_s2"],1e-20) for v in values],"o-",label=profile)
    axes[1,0].set(xlabel="Layer count",ylabel="Rest-state residual acceleration (m/s²)",title="Steep-bed pressure balance");axes[1,0].legend()
    f=cases["baseline"]["frames"][-1]["scalars"]
    labels=["Wind work","Viscosity loss","Bed loss","Implicit loss","Mode correction"]
    keys=["wind_work_J","viscosity_dissipation_J","bed_dissipation_J","implicit_dissipation_J","mode_energy_J"]
    axes[1,1].barh(labels,[f[key]/1000 for key in keys],color=["#0072b2"]*4+["#d55e00"])
    axes[1,1].axvline(0,color="black",linewidth=.7)
    axes[1,1].set(xlabel="Cumulative energy (kJ)",title="Vertical operator and separate mode correction")
    fig.suptitle("Phase 5c · controlled closed-basin diagnostics\nSynthetic geometry; mode correction is not a physical dissipation term",fontsize=14)
    fig.savefig(root/"circulation_review.png",dpi=160);fig.savefig(root/"circulation_review.pdf");plt.close(fig)
    def samples(c,item):
        f=c["frames"][-1];n,k=c["columns"],c["layers"]
        z=np.array(f["interfaces"]).reshape(n,k+1)
        centres=.5*(z[:,:-1]+z[:,1:]);values=np.array(f[item]).reshape(n,k)
        x=(np.arange(n)+.5)*1000/n
        result=[]
        for depth in (-3,-6,-10,-14,-17):
            by_column=[np.interp(depth,centres[i][::-1],values[i][::-1]) for i in range(n)]
            result.extend(np.interp([100,300,500,700,900],x,by_column))
        return np.array(result)
    refinement={}
    # Three separate families: time, layers and horizontal resolution.
    for family,names in {"time":["baseline","dt7_5","dt3_75"],
                         "layers":["layers10","dt3_75","layers40","layers80"],
                         "columns":["columns10","dt3_75","columns40"]}.items():
        refinement[family]={}
        for item in ("temperature","velocity_x","tracer"):
            values=[samples(cases[name],item) for name in names]
            differences=[float(np.sqrt(np.mean((a-b)**2))) for a,b in zip(values,values[1:])]
            refinement[family][item]={"cases":names,"rms_differences":differences,
                "ratios":[a/b if b else None for a,b in zip(differences,differences[1:])]}
    save(root/"refinement.json",{"families":refinement,"sampling":"25 common physical x/z positions; linear centre interpolation within support",
        "interpretation":"Measured finite-grid differences. A ratio below one is unresolved convergence, not a passing accuracy gate."})
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout="constrained")
    for ax,(family,values) in zip(axes,refinement.items()):
        for item in ("temperature","velocity_x","tracer"):
            # Normalize each variable to its first difference; never mix units.
            differences=values[item]["rms_differences"];a=differences[0]
            ax.plot(range(1,len(differences)+1),[b/a if a else 0 for b in differences],"o-",label=item)
        count=len(next(iter(values.values()))["rms_differences"])
        ax.set(xticks=list(range(1,count+1)),xlabel="Successive refinement pairs",ylabel="RMS difference / first difference",title=family)
    axes[0].legend(fontsize=8);fig.suptitle("Refinement at common physical positions · smaller differences need further accuracy review")
    fig.savefig(root/"refinement.png",dpi=160);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",type=Path);parser.add_argument("output",type=Path);parser.add_argument("--figures",action="store_true")
    args=parser.parse_args();data=json.loads(args.source.read_text());root=args.output.resolve();root.mkdir(parents=True,exist_ok=True)
    manifests={c["id"]:make_run(c,root,sha(args.source)) for c in data["cases"]}
    make_composition(root,manifests)
    save(root/"qualification.json",{"cases":{c["id"]:c["frames"][-1]["scalars"] for c in data["cases"]},
        "source_sha256":sha(args.source),"scope":data["scope"],"pressure":data["pressure"]})
    if args.figures: figures(data,root)


if __name__=="__main__": main()
