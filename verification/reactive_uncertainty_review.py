#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bellinge native reactive derivatives, SDK replay and conditional uncertainty.

Recorded water/heat are fixed. Receiver geometry, A->B chemistry, priors and
observations are synthetic. This is not SWMM parameter differentiation, DeGray
calibration or field receiving-water validation.
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
import subprocess
import time
import uuid
from fvqual import Library,Simulation
from fvqual.recorded import open_mesh
from fvqual.scalar_uncertainty import FixedObservation
from fvqual.reactive_uncertainty import local_covariance,gaussian_diagnostics


def save(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def flat(rows):return [x for row in rows for x in row]
def quantile(values,p):
    v=sorted(values);q=(len(v)-1)*p;i=int(q);return v[i]+(v[min(i+1,len(v)-1)]-v[i])*(q-i)
def deck(initial_a=.0002,initial_b=.00005,k=2e-5):
    return f"[SPECIES]\nA kg/m3 {initial_a:.17g}\nB kg/m3 {initial_b:.17g}\n[CONSTANTS]\nk {k:.17g}\ncref .001\n[RATES]\nA WATER -k*A*(1+B/cref)\nB WATER k*A*(1+B/cref)\n"


def run(library,spec,base,initial_a=.0002,initial_b=.00005,k=2e-5,record=False):
    if min(initial_a,initial_b,k)<0:raise ValueError("Negative physical parameter draw")
    with open_mesh(library,{"mesh":spec["mesh"]},base) as mesh:
        with Simulation(mesh,layers=3,temperature=10,options={"advectionScheme":0,"maxTimeStep":300},
                        operations={"reaction_deck":deck(initial_a,initial_b,k)}) as sim:
            if record:sim.begin_reactive_trace()
            initial=sim.mass("A")+sim.mass("B");accepted=0
            frames=[sim.snapshot("A")];a=[frames[0]["values"]];b=[sim.snapshot("B")["values"]]
            for frame in spec["frames"]:
                receipt=sim.apply(frame);accepted+=sum(receipt["species_mass_kg"].values())
                frames.append(sim.snapshot("A"));a.append(frames[-1]["values"]);b.append(sim.snapshot("B")["values"])
            residual=sim.mass("A")+sim.mass("B")-initial-accepted
            ledgers={s:sim.budget("budget_"+s+"_residual") for s in ("A","B")}
            if abs(residual)>1e-8 or max(map(abs,ledgers.values()))>1e-8:raise RuntimeError("Native reactive inventories do not close")
            return {"frames":frames,"A":flat(a),"B":flat(b),"combined_mass_residual_kg":residual,
                    "species_residual_kg":ledgers,"temperature_c":sim.mean_temperature,"volume_m3":sim.volume,
                    "trace":sim.capture_reactive_trace() if record else None}


def build(library,source,output,sdk_driver,samples=256,diagnostics=False):
    output.mkdir(parents=True,exist_ok=False);started=time.perf_counter()
    donor=source/"replay_inflows.json";original=json.loads(donor.read_text());geometry=json.loads((source/"configuration.json").read_text())["receiver_geometry"]
    spec={"mesh":dict(geometry,units="m",elevation_kind="bed-elevation",vertical_datum="synthetic flat receiving bed"),"frames":copy.deepcopy(original["frames"])}
    for frame in spec["frames"]:
        for item in frame["inflows"]:item["species_mass_kg"]={"A":.001*item["volume_m3"],"B":0.}
    save(output/"case.json",spec);save(output/"native_inflows.json",dict(original,schema="fvqual.inflow-series/1.0",species={},frames=spec["frames"]))
    (output/"reaction.rxn").write_text(deck())
    save(output/"sdk_case.json",{"geometry":geometry,"epoch_julian_day":original["epoch_julian_day"],"duration_s":spec["frames"][-1]["end_s"],"frames":spec["frames"],"rate_s_1":2e-5})
    baseline=run(library,spec,output,record=True);trace=baseline.pop("trace")
    with trace:
        indices=[]
        for frame in baseline["frames"]:
            matches=[i for i,t in enumerate(trace.times) if abs(t-frame["time_s"])<1e-8]
            if len(matches)!=1:raise RuntimeError("Ambiguous or absent frame support in reactive history")
            indices.append(matches[0])
        cells=trace.cells;nt=len(indices)
        observations=FixedObservation((trace.steps+1,trace.state_size),[[(t,s*cells+c,1)] for t in indices for s in range(2) for c in range(cells)])
        initial_sd=2e-5;rate_sd=3e-6;b_sd=1e-5;rho=-.35
        modes=[{"initial":{"A":[initial_sd]*cells},"constants":{"k":rho*rate_sd}},
               {"constants":{"k":math.sqrt(1-rho*rho)*rate_sd}},
               {"initial":{"B":[b_sd]*cells}}]
        local=local_covariance(trace,modes,observations)
        def unpack(x,s):return [x[t*2*cells+s*cells+c] for t in range(nt) for c in range(cells)]
        actual=[v for t in range(nt) for s in ("A","B") for v in baseline[s][t*cells:(t+1)*cells]]
        parity=max(abs(a-b) for a,b in zip(actual,local["mean"]))
        if parity>1e-12:raise RuntimeError("Native recorded reactive states differ")
        replay=observations.apply(trace.evaluate())
        if max(abs(a-b) for a,b in zip(actual,replay))>1e-12:raise RuntimeError("Native reactive replay differs")
        direction={"initial":{"A":[initial_sd]*cells,"B":[b_sd]*cells},"constants":{"k":rate_sd}}
        tangent=observations.jvp(trace.jvp(**direction));fd=[]
        for epsilon in (.01,.001,.0001,.00001):
            runs=[run(library,spec,output,.0002+sign*epsilon*initial_sd,.00005+sign*epsilon*b_sd,2e-5+sign*epsilon*rate_sd) for sign in (1,-1)]
            def mapped(r):return [v for t in range(nt) for s in ("A","B") for v in r[s][t*cells:(t+1)*cells]]
            error=max(abs((a-b)/(2*epsilon)-d) for a,b,d in zip(mapped(runs[0]),mapped(runs[1]),tangent))
            if error>1e-9:raise RuntimeError("Independent native reactive directional derivative failed")
            fd.append({"epsilon":epsilon,"max_absolute_concentration_error_kg_m3":error})
        weights=[math.sin(.7*i) for i in range(len(actual))];bar=trace.vjp(observations.vjp(weights))
        left=sum(a*b for a,b in zip(weights,tangent));right=initial_sd*sum(bar["initial"]["A"])+b_sd*sum(bar["initial"]["B"])+rate_sd*bar["constants"]["k"]
        if abs(left-right)>1e-12:raise RuntimeError("Mapped reactive adjoint identity failed")
        rng=random.Random(61017);draws=[[rng.gauss(0,1) for _ in range(3)] for _ in range(samples)]
        def parameters(z):return {"initial_a":.0002+initial_sd*z[0],"initial_b":.00005+b_sd*z[2],"k":2e-5+rate_sd*(rho*z[0]+math.sqrt(1-rho*rho)*z[1])}
        def member(pair):
            index,z=pair;p=parameters(z)
            try:
                r=run(library,spec,output,**p);r.pop("trace");r.pop("frames")
                return dict(r,index=index,z=z,parameters=p,status="completed")
            except Exception as e:return {"index":index,"z":z,"parameters":p,"status":"failed","error":str(e)}
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:members=list(pool.map(member,enumerate(draws)))
        save(output/"ensemble.json",members)
        if any(r["status"]!="completed" for r in members):raise RuntimeError("Ensemble failures retained; review incomplete")
        fields={}
        for s,name in enumerate(("A","B")):
            fields[name]=baseline[name];fields["parameter_sd_"+name]=unpack(local["standard_deviation"],s)
            fields["ensemble_lower_95_"+name]=[quantile([r[name][i] for r in members],.025) for i in range(nt*cells)]
            fields["ensemble_upper_95_"+name]=[quantile([r[name][i] for r in members],.975) for i in range(nt*cells)]
            for sign,label in ((-1,"local_lower_95_"),(1,"local_upper_95_")):
                fields[label+name]=[m+sign*1.96*sd for m,sd in zip(fields[name],fields["parameter_sd_"+name])]
            fields["kinetic_response_"+name]=unpack(observations.jvp(trace.jvp(constants={"k":rate_sd})),s)
            fields["initial_response_"+name]=unpack(observations.jvp(trace.jvp(initial={"A":[initial_sd]*cells})),s)
        subprocess.run([str(sdk_driver),str(output)],check=True)
        sdk=json.loads((output/"sdk_qualification.json").read_text())
        for i,row in enumerate(sdk["every_step"]["history"]):
            for name in ("A","B"):
                if max(abs(a-b) for a,b in zip(row[name],baseline[name][(i+1)*cells:(i+2)*cells]))>1e-12:raise RuntimeError("C ABI and SDK native trajectories differ")
        endpoint_weights=[]
        for i in range(nt):endpoint_weights.extend(([1.]*cells+[.3]*cells) if i else [0.]*(2*cells))
        sdk_bar=trace.vjp(observations.vjp(endpoint_weights));sdk_difference=abs(sdk_bar["constants"]["k"]-sdk["every_step"]["constant_gradient"]["k"])
        if sdk_difference>1e-8:raise RuntimeError("C ABI and SDK gradient routing disagree")
        report={"status":"qualified conditional native reactive transport and SDK integration","samples":samples,"successful":samples,"failed":0,"seed":61017,
            "scope":"native mixed-mesh upwind transport, implicit vertical diffusion and nonlinear multispecies adaptive kinetics; initial species and kinetic constants",
            "prescribed":["SWMM water/heat exports","hydraulic/thermal trajectory","source composition and placement","geometry","mixing coefficients"],
            "unsupported":["horizontal thermal/momentum derivatives","general geometry/withdrawal placement derivatives","settling derivatives","MUSCL/FCT","self-shading","SWMM parameter gradients","wet/dry"],
            "uncertainty_scope":"Illustrative parameter-only marginal bands; receiver, chemistry and priors synthetic. Observation and structural discrepancy excluded; not field validation.",
            "priors":{"distribution":"joint Gaussian; negative physical draws fail and are counted","initial_A_mean_kg_m3":.0002,"initial_A_sd":initial_sd,"initial_B_mean_kg_m3":.00005,"initial_B_sd":b_sd,"k_mean_s_1":2e-5,"k_sd":rate_sd,"initial_A_k_correlation":rho},
            "whole_native_directional_sweep":fd,"mapped_adjoint_residual":abs(left-right),"recorded_forward_max_difference":parity,
            "max_mass_residual_kg":max(abs(r["combined_mass_residual_kg"]) for r in members),"max_species_budget_residual_kg":max(abs(v) for r in members for v in r["species_residual_kg"].values()),
            "temperature_max_difference_c":max(abs(r["temperature_c"]-baseline["temperature_c"]) for r in members),"volume_max_difference_m3":max(abs(r["volume_m3"]-baseline["volume_m3"]) for r in members),
            "sdk_parameter_gradient_difference":sdk_difference,"sdk_checkpoint_comparison":{p:{key:sdk[p][key] for key in ("stored_checkpoints","peak_checkpoints","replays")} for p in ("every_step","logarithmic")},
            "peak_cell_tape_nodes":bar["peak_cell_tape_nodes"],"checkpoint_scope":"composition checkpoint counts; latest-interval operator history remains bounded separately, not a process peak-memory benchmark",
            "provenance":{"library_sha256":sha(library.path),"driver_sha256":sha(__file__),"sdk_driver_sha256":sha(sdk_driver),"donor_sha256":sha(donor),"configuration_sha256":sha(output/"case.json")}}
        if diagnostics:
            requests=[(indices[t],s*cells+c,1) for t,s,c in ((3,0,1),(8,0,4),(14,0,7),(20,0,1),(8,1,1),(14,1,4),(20,1,7))]
            survey=FixedObservation((trace.steps+1,trace.state_size),[[r] for r in requests]);survey_local=local_covariance(trace,modes,survey)
            truth=run(library,spec,output,.00021,.000055,2.2e-5)
            frame_map={t:i for i,t in enumerate(indices)}
            def sample(r):return [r["A" if c<cells else "B"][frame_map[t]*cells+c%cells] for t,c,_ in requests]
            noise=random.Random(61018);observed=[x+noise.gauss(0,1e-5) for x in sample(truth)]
            design=gaussian_diagnostics(survey_local,[1e-5]*len(requests),observed)
            posterior=run(library,spec,output,**parameters(design["posterior_mode_mean"]))
            heldout=[(t,c) for t in (5,11,17,20) for c in (0,5,8)]
            def rmse(r):return math.sqrt(sum((r["A"][t*cells+c]-truth["A"][t*cells+c])**2 for t,c in heldout)/len(heldout))
            design.update(assumptions="Synthetic independent observation SD .00001 kg/m3; known synthetic truth and no structural error. Local Gaussian demonstration, not DeGray calibration.",
                requests=requests,observed_values=observed,heldout_requests=heldout,baseline_heldout_rmse_kg_m3=rmse(baseline),posterior_native_heldout_rmse_kg_m3=rmse(posterior),
                heldout_improved=rmse(posterior)<rmse(baseline))
            save(output/"observation_design.json",design);report["observation_design"]={k:design[k] for k in ("numerical_rank","mode_count","condition_number","baseline_heldout_rmse_kg_m3","posterior_native_heldout_rmse_kg_m3")}
        report["wall_seconds"]=time.perf_counter()-started
        save(output/"qualification.json",report);save(output/"review_fields.json",{"frames":baseline["frames"],"fields":fields,"qualification":report})
        make_gui(output,source,baseline["frames"],fields,original["epoch_julian_day"],report)
        return report


def make_gui(root,source,frames,fields,epoch,report):
    old=json.loads((source/"live_output/run.json").read_text());reference=next(r for r in old["results"] if r["component"]=="reservoir" and r["item"]=="temperature")
    save(root/"geometry.json",json.loads((source/"live_output"/reference["geometry"]["artifact"]).read_text()))
    with (root/"fields.csv").open("w") as stream:
        writer=csv.writer(stream);writer.writerow(("time","item","index","value"))
        for t,frame in enumerate(frames):
            clock=epoch+frame["time_s"]/86400
            for name,values in fields.items():
                for cell in range(9):writer.writerow((clock,name,cell,values[t*9+cell]))
            for i,z in enumerate(frame["interfaces"]):writer.writerow((clock,"interfaces",i,z))
    clock={"count":len(frames),"start":epoch,"end":epoch+frames[-1]["time_s"]/86400,"units":"julian_day"}
    vertical={"component":"reservoir","item":"interfaces","artifact":"fields.csv","format":"csv","variable":"interfaces","kind":"Float64","shape":[len(frames),3,4],"dimensions":["time","columns","interfaces"],"time":clock,"units":"m"}
    quantity={"scale":1,"offset":0,"value_kind":1,"powers":[-3,1,0,0,0,0,0,0,0,0]}
    entries=[]
    for name in fields:
        entries.append({"component":"reservoir","item":name,"artifact":"fields.csv","format":"csv","variable":name,"kind":"Float64","shape":[len(frames),3,3],"dimensions":["time","columns","layers"],"time":clock,"units":"kg/m3",
            "location":"face","geometry":{"artifact":"geometry.json","entity_axis":1,"layer_axis":2,"vertical_kind":1,"vertical":vertical},"quantity":quantity,"description":report["uncertainty_scope"]})
    caption="Bellinge — synthetic A to B reaction and conditional uncertainty"
    def manifest(path,entries,caption,provenance):
        run={"id":str(uuid.uuid4()),"caption":caption,"status":"completed"};save(path,{"schema_version":"1.0","run":run,"components":[{"id":"reservoir","caption":caption,"final_status":"Finished"}],"results":entries,"provenance":provenance});return run
    path=root/"run.json";run=manifest(path,entries,caption,report)
    layer={"type":"recorded-result","manifest":str(path),"component":"reservoir","item":"A","name":caption,"units":"kg/m3","opacity":1,"visible":True,"verticalLayer":1,"scene":{"shownIn3D":True,"terrainEnabled":True,"zMode":2,"offset":0},"runId":run["id"],"runTitle":caption,"runStatus":"completed"}
    requests=[]
    for item,label,color in (("A","A baseline","#0072b2"),("ensemble_lower_95_A","A sample lower 95%","#009e73"),("ensemble_upper_95_A","A sample upper 95%","#d55e00")):
        view=path;record=run
        if item!="A":
            view=root/(item+".json");entry=copy.deepcopy(entries[0]);entry["variable"]=item
            record=manifest(view,[entry],caption+" — "+label,{"derived_from":str(path),"scope":"Marginal sample envelope of A, not an independent trajectory"})
        requests.append(dict(layer,type="recorded-plot",manifest=str(view),runId=record["id"],runTitle=record["caption"],feature=0,layerCount=3,legendLabel=label,plotName=label,plotColor=color,plotVisible=True,plotPoints=True,unitSignature=dict(quantity,caption="kg/m3")))
    composition={"schema_version":"1.1","metadata":{"id":"reactive-uncertainty-review","caption":caption},"components":[{"id":"reservoir","caption":caption,"execution":{"mode":"open","results":str(path)}}]}
    save(root/"Bellinge_reactive_review.composition.json",composition)
    save(root/"Bellinge_reactive_review.composition.composer.json",{"sidecar_version":1,"layers":[layer],"results":{"version":1,"runs":[r["manifest"] for r in requests],"alignment":0,"active":requests[0],"plot":{"version":1,"pinned":True,"requests":requests}},"navigation":{"version":1,"workspace":"results"}})
    save(root/"Bellinge_reactive_native.composition.json",composition)
    save(root/"Bellinge_reactive_native.composition.composer.json",{"sidecar_version":1,"layers":[layer],"results":{"version":1,"runs":[str(path)],"alignment":0}})


def figures(root):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data=json.loads((root/"review_fields.json").read_text());f={k:np.asarray(v).reshape(-1,3,3)*1000 for k,v in data["fields"].items()};times=np.array([r["time_s"] for r in data["frames"]])/3600
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout="constrained");colors=("#0072b2","#d55e00","#009e73")
    for s,color in zip(("A","B"),colors):axes[0,0].plot(times,f[s][:,0,1],color=color,label=s)
    axes[0,0].set(title="Synthetic conversion A → B · face 0, middle layer",xlabel="Elapsed hours",ylabel="Concentration (mg/L)");axes[0,0].legend()
    ax=axes[0,1];ax.fill_between(times,f["ensemble_lower_95_A"][:,0,1],f["ensemble_upper_95_A"][:,0,1],color=colors[0],alpha=.2,label="Native sample 2.5 / 97.5%")
    ax.plot(times,f["A"][:,0,1],color=colors[0],label="Baseline")
    for item in ("local_lower_95_A","local_upper_95_A"):ax.plot(times,f[item][:,0,1],"--",color=colors[1],label="Local Gaussian 95%" if "lower" in item else None)
    ax.set(title="Conditional parameter envelopes for A",xlabel="Elapsed hours",ylabel="Concentration (mg/L)");ax.legend(fontsize=8)
    for name,color in zip(("kinetic_response_A","initial_response_A"),colors):axes[1,0].plot(times,f[name][:,0,1],color=color,label=name.replace("_"," "))
    axes[1,0].set(title="Analytic responses to one parameter SD",xlabel="Elapsed hours",ylabel="Response (mg/L)");axes[1,0].legend(fontsize=8)
    design=json.loads((root/"observation_design.json").read_text());axes[1,1].bar(range(1,len(design["singular_values"])+1),design["singular_values"],color=colors[0])
    axes[1,1].set(title="Synthetic observation design · whitened sensitivities",xlabel="Sensitivity direction (rank order)",ylabel="Singular value");axes[1,1].set_xticks(range(1,len(design["singular_values"])+1))
    fig.suptitle("Bellinge reactive transport and SDK uncertainty review\nSynthetic basin, chemistry, priors and observations. Recorded water/heat fixed; structural discrepancy excluded.",fontsize=12)
    fig.savefig(root/"reactive_review.png",dpi=160);fig.savefig(root/"reactive_review.pdf");plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("source",type=Path);parser.add_argument("output",type=Path);parser.add_argument("--sdk-driver",type=Path,required=True);parser.add_argument("--library");parser.add_argument("--samples",type=int,default=256);parser.add_argument("--figures",action="store_true")
    args=parser.parse_args()
    if args.samples<32:raise ValueError("At least 32 samples required")
    report=build(Library(args.library),args.source.resolve(),args.output.resolve(),args.sdk_driver.resolve(),args.samples,args.figures)
    if args.figures:figures(args.output.resolve())
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()
