#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Qualify actual passive-scalar gradients and build a Bellinge GUI review.

The receiver geometry and tracer are synthetic. Recorded Bellinge water/heat
are fixed; this does not differentiate SWMM parameters or qualify DeGray heat.
Standard-library workflow; --figures additionally requires matplotlib/numpy.
"""
import argparse
import concurrent.futures
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import time
import uuid

from fvqual import Library,Simulation
from fvqual.recorded import open_mesh,SCHEMA,UNITS
from fvqual.scalar_uncertainty import FixedObservation,local_covariance


def save(path,data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def run(library,spec,base,initial=.0002,load=1,record=False):
    if initial<0 or load<0: raise ValueError("Negative concentration or load in physical ensemble")
    with open_mesh(library,spec,base) as mesh:
        with Simulation(mesh,layers=spec["layers"],temperature=10,
                        species={"tracer":{"units":"kg/m3","initial":initial}},options=spec["options"]) as sim:
            if record:sim.begin_scalar_trace("tracer")
            initial_mass=sim.mass("tracer");accepted=0
            snapshots=[sim.snapshot("tracer")]
            for original in spec["frames"]:
                frame=copy.deepcopy(original)
                for data in frame["inflows"]:data["species_mass_kg"]["tracer"]*=load
                receipt=sim.apply(frame);accepted+=receipt["species_mass_kg"]["tracer"]
                snapshots.append(sim.snapshot("tracer"))
            residual=sim.mass("tracer")-initial_mass-accepted
            if abs(residual)>1e-8:raise RuntimeError("Tracer mass ledger does not close")
            return {"snapshots":snapshots,"mass_residual_kg":residual,
                    "volume_m3":sim.volume,"temperature_c":sim.mean_temperature,
                    "trace":sim.capture_scalar_trace() if record else None}


def percentile(values,p):
    values=sorted(values);index=(len(values)-1)*p;lo=int(index);hi=min(lo+1,len(values)-1)
    return values[lo]+(values[hi]-values[lo])*(index-lo)


def build_review(library,source,output,samples=128,seed=60213):
    if samples<32:raise ValueError("At least 32 samples required for convergence review")
    output.mkdir(parents=True,exist_ok=False)
    donor=source/"replay_inflows.json";config=source/"configuration.json"
    original=json.loads(donor.read_text());geometry=json.loads(config.read_text())["receiver_geometry"]
    spec={"schema":SCHEMA,"units":UNITS,"epoch_julian_day":original["epoch_julian_day"],
          "mesh":dict(geometry,units="m",elevation_kind="bed-elevation",vertical_datum="synthetic flat bed"),
          "layers":3,"temperature_c":10,"species":{"tracer":{"units":"kg/m3","initial":.0002}},
          "options":{"advectionScheme":0,"maxTimeStep":300,"backgroundDiffusivity":1e-4},
          "frames":copy.deepcopy(original["frames"])}
    for frame in spec["frames"]:
        for data in frame["inflows"]:data["species_mass_kg"]={"tracer":.001*data["volume_m3"]}
    save(output/"recorded_tracer_case.json",spec)
    started=time.perf_counter();baseline=run(library,spec,output,record=True);trace=baseline.pop("trace")
    with trace:
        n=trace.cells;frame_indices=[]
        for snap in baseline["snapshots"]:
            matches=[i for i,t in enumerate(trace.times) if abs(t-snap["time_s"])<1e-8]
            if len(matches)!=1:raise RuntimeError("Observation time has ambiguous or absent step support")
            frame_indices.append(matches[0])
        observations=FixedObservation((trace.steps+1,n),[[(i,c,1)] for i in frame_indices for c in range(n)])
        rho=-.35;load_sd=.2;initial_sd=.00002
        modes=[{"initial":[initial_sd]*n,"sources":[rho*load_sd*x for x in trace.sources]},
               {"sources":[math.sqrt(1-rho*rho)*load_sd*x for x in trace.sources]}]
        local=local_covariance(trace,modes,observations)
        # Independent full solver runs, not finite differences of a trace replay.
        direction={"initial":[initial_sd]*n,"sources":[load_sd*x for x in trace.sources]}
        tangent=observations.jvp(trace.jvp(**direction));fd=[]
        for epsilon in (.01,.001,.0001,.00001,.000001):
            plus=run(library,spec,output,.0002+epsilon*initial_sd,1+epsilon*load_sd)
            minus=run(library,spec,output,.0002-epsilon*initial_sd,1-epsilon*load_sd)
            a=[x for s in plus["snapshots"] for x in s["values"]]
            b=[x for s in minus["snapshots"] for x in s["values"]]
            error=max(abs((x-y)/(2*epsilon)-d) for x,y,d in zip(a,b,tangent))
            fd.append({"epsilon":epsilon,"max_absolute_error_kg_m3":error})
            if error>2e-11:raise RuntimeError("Whole-run finite difference gate failed")
        rng=random.Random(seed);w=[rng.uniform(-1,1) for _ in tangent]
        bar=trace.vjp(observations.vjp(w))
        lhs=math.fsum(x*y for x,y in zip(w,tangent))
        rhs=math.fsum(x*y for x,y in zip(direction["initial"],bar["initial"]))+math.fsum(x*y for x,y in zip(direction["sources"],bar["sources"]))
        dot_error=abs(lhs-rhs)
        if dot_error>1e-13:raise RuntimeError("Mapped adjoint identity failed")
        replay=trace.evaluate();actual_error=max(abs(a-b) for i in frame_indices for a,b in zip(replay[i],trace.states[i]))
        if actual_error>2e-14:raise RuntimeError("Recorded forward replay differs from production solver")
        # Same draws for all prefix comparisons. Correlation is implemented via L.
        rng=random.Random(seed);draws=[(rng.gauss(0,1),rng.gauss(0,1)) for _ in range(samples)]
        def member(pair):
            index,(z0,z1)=pair
            initial=.0002+initial_sd*z0
            load=1+load_sd*(rho*z0+math.sqrt(1-rho*rho)*z1)
            try:
                result=run(library,spec,output,initial,load)
                values=[x for snap in result["snapshots"] for x in snap["values"]]
                predicted=[m+a*z0+b*z1 for m,a,b in zip(local["mean"],*local["mode_sensitivities"])]
                difference=max(abs(a-b) for a,b in zip(values,predicted))
                result.pop("trace");result.update(index=index,z=[z0,z1],initial_kg_m3=initial,load_multiplier=load,
                    max_tangent_prediction_error_kg_m3=difference,status="completed")
                return result
            except Exception as error:return {"index":index,"status":"failed","error":str(error),"z":[z0,z1]}
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            members=list(executor.map(member,enumerate(draws)))
        save(output/"ensemble_runs.json",members)
        failures=[r for r in members if r["status"]!="completed"]
        if failures:raise RuntimeError("Ensemble failures retained in ensemble_runs.json; review incomplete")
        max_prediction=max(r["max_tangent_prediction_error_kg_m3"] for r in members)
        if max_prediction>2e-14:raise RuntimeError("Actual whole-run ensemble differs from tangent prediction")
        member_values=[[x for s in r["snapshots"] for x in s["values"]] for r in members]
        lower=[percentile([v[i] for v in member_values],.025) for i in range(len(tangent))]
        upper=[percentile([v[i] for v in member_values],.975) for i in range(len(tangent))]
        convergence=[]
        for count in sorted(set([32,64,samples])):
            if count>samples:continue
            errors=[]
            for i,sd in enumerate(local["standard_deviation"]):
                a=[v[i] for v in member_values[:count]];mean=math.fsum(a)/count
                empirical=math.sqrt(math.fsum((v-mean)**2 for v in a)/(count-1))
                if sd>1e-14:errors.append(abs(empirical/sd-1))
            convergence.append({"samples":count,"max_relative_standard_deviation_difference":max(errors)})
        # Field responses to one SD of each marginal parameter, before correlation.
        initial_response=observations.jvp(trace.jvp(initial=[initial_sd]*n))
        load_response=observations.jvp(trace.jvp(sources=[load_sd*x for x in trace.sources]))
        report={"status":"qualified restricted passive-scalar gate", "coverage":"actual upwind advection and implicit vertical diffusion; multi-step initialization and gross source concentration derivatives; fixed observation transpose",
            "frozen":["hydraulics","geometry","placement","Kz trajectory","substep schedule","recorded SWMM water and heat"],
            "unsupported":["temperature/wind/geometry derivatives","reaction/settling derivatives","MUSCL/FCT derivatives","SWMM parameter gradients","feedback","posterior calibration"],
            "science":"Bellinge water/heat forcing; synthetic flat mixed-mesh receiver and synthetic conservative tracer. Parameter-only envelopes, not predictive uncertainty or observed Bellinge quality.",
            "parameters":{"initial_kg_m3":{"mean":.0002,"sd":initial_sd},"load_multiplier":{"mean":1,"sd":load_sd},"correlation":rho,"distribution":"joint Gaussian; negative physical draws fail and are counted"},
            "seed":seed,"samples":samples,"successful":len(members),"failed":len(failures),"convergence":convergence,
            "whole_run_finite_differences":fd,"observation_adjoint_dot_error":dot_error,"recorded_forward_max_difference":actual_error,
            "whole_run_ensemble_max_tangent_prediction_error":max_prediction,
            "max_mass_residual_kg":max(abs(r["mass_residual_kg"]) for r in members),
            "hydraulic_volume_max_difference_m3":max(abs(r["volume_m3"]-baseline["volume_m3"]) for r in members),
            "temperature_max_difference_c":max(abs(r["temperature_c"]-baseline["temperature_c"]) for r in members),
            "internal_steps":trace.steps,"trace_payload_bytes":trace.storage_bytes,"wall_seconds":time.perf_counter()-started,
            "provenance":{"library_sha256":sha(library.path),"driver_sha256":sha(Path(__file__)),"donor_sha256":sha(donor),
                          "configuration_sha256":sha(config),"case_sha256":sha(output/"recorded_tracer_case.json")}}
        save(output/"qualification.json",report)
        fields={"tracer":local["mean"],"parameter_sd":local["standard_deviation"],
                "parameter_lower_95":[m-1.96*s for m,s in zip(local["mean"],local["standard_deviation"])],
                "parameter_upper_95":[m+1.96*s for m,s in zip(local["mean"],local["standard_deviation"])],
                "ensemble_lower_95":lower,"ensemble_upper_95":upper,
                "initial_sensitivity":initial_response,"load_sensitivity":load_response}
        save(output/"review_fields.json",{"frames":baseline["snapshots"],"fields":fields,"qualification":report})
        make_gui(output,source,spec,baseline["snapshots"],fields,report)
        return report


def make_gui(root,source,spec,snapshots,fields,report):
    n,k=snapshots[0]["columns"],snapshots[0]["layers"];cells=n*k
    # Reuse the registered geometric identity/CRS, not invented coordinate units.
    manifest=json.loads((source/"live_output/run.json").read_text())
    result=next(r for r in manifest["results"] if r["component"]=="reservoir" and r["item"]=="temperature")
    geo=source/"live_output"/result["geometry"]["artifact"]
    geometry=json.loads(geo.read_text());save(root/"geometry.json",geometry)
    epoch=spec["epoch_julian_day"]
    with (root/"fields.csv").open("w") as stream:
        writer=csv.writer(stream);writer.writerow(("time","item","index","value"))
        for t,snap in enumerate(snapshots):
            clock=epoch+snap["time_s"]/86400
            for item,values in fields.items():
                for i in range(cells):writer.writerow((clock,item,i,values[t*cells+i]))
            for i,z in enumerate(snap["interfaces"]):writer.writerow((clock,"interfaces",i,z))
    clock={"count":len(snapshots),"start":epoch,"end":epoch+snapshots[-1]["time_s"]/86400,"units":"julian_day"}
    def entry(item,shape,dimensions,units):
        return {"component":"reservoir","item":item,"artifact":"fields.csv","format":"csv","variable":item,
                "kind":"Float64","shape":[len(snapshots)]+shape,"dimensions":["time"]+dimensions,"time":clock,"units":units}
    vertical=entry("interfaces",[n,k+1],["columns","interfaces"],"m");results=[]
    for item in fields:
        r=entry(item,[n,k],["columns","layers"],"kg/m3")
        r.update(location="face",geometry={"artifact":"geometry.json","entity_axis":1,"layer_axis":2,"vertical_kind":1,"vertical":vertical},
                 quantity={"scale":1,"offset":0,"value_kind":1,"powers":[-3,1,0,0,0,0,0,0,0,0]},
                 description="Synthetic passive tracer; fixed Bellinge water/heat. Parameter-only uncertainty, no structural or observation error.")
        results.append(r)
    run_id=str(uuid.uuid4());caption="Bellinge — passive tracer sensitivity review"
    save(root/"run.json",{"schema_version":"1.0","run":{"id":run_id,"caption":caption,"status":"completed"},
        "components":[{"id":"reservoir","caption":caption,"final_status":"Finished"}],"results":results,"provenance":report})
    composition=root/"Bellinge_sensitivity_review.composition.json"
    save(composition,{"schema_version":"1.1","metadata":{"id":"phase6a-bellinge","caption":caption},
        "components":[{"id":"reservoir","caption":caption,"execution":{"mode":"open","results":str(root/"run.json")}}]})
    layer={"type":"recorded-result","manifest":str(root/"run.json"),"component":"reservoir","item":"tracer",
        "name":"Bellinge — synthetic tracer","units":"kg/m3","visible":True,"opacity":1,"verticalLayer":1,
        "scene":{"shownIn3D":True,"terrainEnabled":True,"zMode":2,"offset":0},"runId":run_id,"runTitle":caption,"runStatus":"completed"}
    requests=[]
    for item,label,color in [("tracer","Baseline","#0072b2"),("parameter_lower_95","Parameter lower 95%","#009e73"),("parameter_upper_95","Parameter upper 95%","#d55e00")]:
        # Envelopes describe the same physical variable. Distinct recorded
        # views preserve identity while satisfying the plot's variable/unit
        # guard; never relabel a sensitivity as a concentration to bypass it.
        view=root/"run.json";view_id=run_id;view_title=caption
        if item!="tracer":
            view=root/("run_"+item+".json");view_id=str(uuid.uuid4());view_title=caption+" — "+label
            envelope=copy.deepcopy(results[0]);envelope["variable"]=item
            envelope["description"]="Derived marginal parameter envelope of tracer concentration; not an independently simulated trajectory. Structural/observation error excluded."
            save(view,{"schema_version":"1.0","run":{"id":view_id,"caption":view_title,"status":"completed"},
                "components":[{"id":"reservoir","caption":view_title,"final_status":"Finished"}],"results":[envelope],
                "provenance":{"derived_from":str(root/"run.json"),"parent_run_id":run_id,"envelope":item}})
        requests.append(dict(layer,type="recorded-plot",item="tracer",manifest=str(view),runId=view_id,runTitle=view_title,feature=0,layerCount=k,legendLabel=label,
            plotName=label+" — face 0, layer 1",plotColor=color,plotVisible=True,plotPoints=True,
            unitSignature={"caption":"kg/m3","scale":1,"offset":0,"powers":[-3,1,0,0,0,0,0,0,0,0]}))
    save(root/"Bellinge_sensitivity_review.composition.composer.json",{"sidecar_version":1,"layers":[layer],
        "results":{"version":1,"runs":[request["manifest"] for request in requests],"alignment":0,"active":requests[0],
                   "plot":{"version":1,"pinned":True,"requests":requests}},"navigation":{"version":1,"workspace":"results"}})
    # Native acceptance starts with an unpinned single map layer.
    save(root/"Bellinge_sensitivity_native.composition.json",json.loads(composition.read_text()))
    save(root/"Bellinge_sensitivity_native.composition.composer.json",{"sidecar_version":1,"layers":[layer],
        "results":{"version":1,"runs":[str(root/"run.json")],"alignment":0}})


def figures(root):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data=json.loads((root/"review_fields.json").read_text());fields=data["fields"];frames=data["frames"]
    times=np.array([f["time_s"] for f in frames])/3600;shape=(len(frames),frames[0]["columns"],frames[0]["layers"])
    f={key:np.array(value).reshape(shape)*1000 for key,value in fields.items()}
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    colors=("#0072b2","#d55e00","#009e73")
    for c,color in enumerate(colors):axes[0,0].plot(times,f["tracer"][:,c,1],color=color,label=f"Face {c}, layer 1")
    axes[0,0].legend();axes[0,0].set(title="Synthetic tracer on recorded Bellinge forcing",xlabel="Elapsed time (hours)",ylabel="Concentration (mg/L)")
    ax=axes[0,1];ax.fill_between(times,f["parameter_lower_95"][:,0,1],f["parameter_upper_95"][:,0,1],color=colors[0],alpha=.18,label="Local Gaussian parameter 95%")
    ax.plot(times,f["tracer"][:,0,1],color=colors[0],label="Baseline")
    ax.plot(times,f["ensemble_lower_95"][:,0,1],"--",color=colors[1],label="Whole-run sample 2.5 / 97.5%")
    ax.plot(times,f["ensemble_upper_95"][:,0,1],"--",color=colors[1]);ax.legend(fontsize=8)
    ax.set(title="Face 0, layer 1 · parameter envelopes",xlabel="Elapsed time (hours)",ylabel="Concentration (mg/L)")
    for c,color in enumerate(colors):
        axes[1,0].plot(times,f["initial_sensitivity"][:,c,1],color=color,label=f"Initial · face {c}")
        axes[1,0].plot(times,f["load_sensitivity"][:,c,1],"--",color=color,label=f"Load · face {c}")
    axes[1,0].legend(fontsize=8,ncol=2);axes[1,0].set(title="Responses to one parameter standard deviation",xlabel="Elapsed time (hours)",ylabel="Concentration response (mg/L)")
    report=data["qualification"];values=report["convergence"]
    axes[1,1].plot([v["samples"] for v in values],[100*v["max_relative_standard_deviation_difference"] for v in values],"o-",color=colors[0])
    axes[1,1].set(title="Sample covariance convergence diagnostic",xlabel="Whole-run samples",ylabel="Maximum SD difference from local Gaussian (%)")
    fig.suptitle("Phase 6a · actual passive-scalar derivatives\nSynthetic receiver/tracer; fixed hydraulics. Bands omit structural and observation uncertainty.",fontsize=13)
    fig.savefig(root/"sensitivity_review.png",dpi=160);fig.savefig(root/"sensitivity_review.pdf");plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("source",type=Path);parser.add_argument("output",type=Path)
    parser.add_argument("--library");parser.add_argument("--samples",type=int,default=128);parser.add_argument("--figures",action="store_true")
    args=parser.parse_args();report=build_review(Library(args.library),args.source.resolve(),args.output.resolve(),args.samples)
    if args.figures:figures(args.output)
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()
