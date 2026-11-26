#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Combined thermal/geometry/replay review and Bellinge full-run uncertainty.

Controlled-column derivatives are analytic. Bellinge sensitivities are named
finite-difference references; its parameter bands come from full native runs.
No DeGray calibration, general reservoir adjoint or SWMM gradient is claimed.
"""
import argparse
import concurrent.futures
import copy
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import random
import time
import uuid
from fvqual import Library,ThermalColumn,SigmaDepthObservation,Simulation
from fvqual.recorded import open_mesh


def save(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def flat(rows):return [x for row in rows for x in row]
def percentile(values,p):
    values=sorted(values);q=(len(values)-1)*p;i=int(q);j=min(i+1,len(values)-1)
    return values[i]+(values[j]-values[i])*(q-i)
def bands(rows):
    return ([percentile([r[i] for r in rows],.025) for i in range(len(rows[0]))],
            [percentile([r[i] for r in rows],.975) for i in range(len(rows[0]))])
def field(rows,layers,block=0):return [v for row in rows for v in row[block*layers:(block+1)*layers]]


def column_review(library,output,samples,seed):
    config={"layers":8,"depth_m":10,"surface_area_m2":10000,"bed_area_fraction":.35,"duration_s":21600,"dt_s":60}
    started=time.perf_counter()
    with ThermalColumn(library,config) as model:
        n=model.layers;shape=(model.time_count,model.state_size)
        observation=SigmaDepthObservation(shape,n,10,[(model.time_count-1,2.6),(model.time_count-1,7.4),(model.time_count//2,2.6)])
        direction={"wind_speed_m_s":.3,"depth_scale":.025,"air_temperature_c":.5,"solar_attenuation_m_1":.03}
        dx=[.2]*n+[.001*math.sin(i) for i in range(2*n)]
        tangent=model.jvp(dx,direction);mapped=observation.jvp(model.states,tangent,.025)
        fd=[]
        for epsilon in (.01,.001,.0001,.00001):
            values=[]
            for sign in (1,-1):
                initial=[x+sign*epsilon*v for x,v in zip(model.initial,dx)]
                p={key:model.parameters[key]+sign*epsilon*v for key,v in direction.items()}
                values.append(observation.apply(model.evaluate(initial,p),1+sign*epsilon*.025))
            error=max(abs((a-b)/(2*epsilon)-d) for a,b,d in zip(*values,mapped))
            if error>3e-6:raise RuntimeError("Native thermal/geometry directional derivative failed")
            fd.append({"epsilon":epsilon,"max_absolute_temperature_error_c":error})
        weights=[.3,-.7,.2];obsbar=observation.vjp(model.states,weights);bar=model.vjp(obsbar["states"])
        left=sum(a*b for a,b in zip(weights,mapped))
        right=sum(a*b for a,b in zip(dx,bar["initial"]))+sum(v*bar["parameters"][key] for key,v in direction.items())+.025*obsbar["depth_scale"]
        if abs(left-right)>1e-10:raise RuntimeError("Moving-depth adjoint identity failed")
        # Independent actual solver runs for all covariance-factor samples.
        modes=[([.3]*n+[0.]*(2*n),{}),([0.]*model.state_size,{"wind_speed_m_s":.3,"air_temperature_c":.35}),
               ([0.]*model.state_size,{"air_temperature_c":math.sqrt(1-.35**2)}),
               ([0.]*model.state_size,{"depth_scale":.025}),([0.]*model.state_size,{"solar_attenuation_m_1":.03})]
        tangents=[model.jvp(x,p) for x,p in modes]
        mean=field(model.states,n);columns=[field(t,n) for t in tangents]
        sd=[math.sqrt(sum(col[i]**2 for col in columns)) for i in range(len(mean))]
        rng=random.Random(seed);draws=[[rng.gauss(0,1) for _ in modes] for _ in range(samples)]
        def member(pair):
            index,z=pair
            try:
                initial=[x+sum(v*mode[0][i] for v,mode in zip(z,modes)) for i,x in enumerate(model.initial)]
                p={name:model.parameters[name]+sum(v*mode[1].get(name,0) for v,mode in zip(z,modes)) for name in model.parameter_names}
                native=model.evaluate(initial,p);values=field(native,n)
                return {"index":index,"z":z,"parameters":p,"temperature":values,
                    "instrument_temperature_c":observation.apply(native,p["depth_scale"]),"status":"completed"}
            except Exception as error:return {"index":index,"z":z,"status":"failed","error":str(error)}
        # Calls on one handle are serialized; each native rerun owns its solver.
        members=[member(pair) for pair in enumerate(draws)]
        save(output/"column_ensemble.json",members)
        failures=[m for m in members if m["status"]!="completed"]
        if failures:raise RuntimeError("Column ensemble failures retained; qualification incomplete")
        lower,upper=bands([m["temperature"] for m in members])
        all_weights=[[0.]*model.state_size for _ in model.times];all_weights[-1][n//2]=1
        default_bar=model.vjp(all_weights)
        with ThermalColumn(library,dict(config,checkpoint_stride=1)) as full:
            fullbar=full.vjp(all_weights)
            difference=max(abs(default_bar["parameters"][k]-fullbar["parameters"][k]) for k in model.parameter_names)
            if difference>1e-9:raise RuntimeError("Checkpoint policy gradients differ")
            policies={"every_step_checkpoints":full.checkpoint_count,"blocked_checkpoints":model.checkpoint_count,
                "max_parameter_gradient_difference":difference,"blocked":default_bar["stats"],"every_step":fullbar["stats"],
                "scope":"checkpoint and tape payload counts, not process peak memory or logarithmic scheduling"}
        instrument_columns=[observation.jvp(model.states,t,mode[1].get("depth_scale",0)) for t,mode in zip(tangents,modes)]
        instrument_sd=[math.sqrt(sum(col[i]**2 for col in instrument_columns)) for i in range(3)]
        report={"status":"qualified restricted thermal/geometry/replay gate","configuration":config,"samples":samples,"successful":len(members),"failed":0,
            "seed":seed,"mode_definition":[{"initial":x,"parameters":p} for x,p in modes],
            "parameter_assumptions":"Gaussian: initial sigma-profile offset SD .3 C; wind SD .3 m/s; air SD 1 C (wind/air correlation .35); depth scale SD .025; attenuation SD .03/m",
            "maximum_native_operator_difference":model.forward_difference,"minimum_buoyancy_margin":model.minimum_buoyancy_margin,
            "whole_native_directional_sweep":fd,"moving_observation_adjoint_residual":abs(left-right),"checkpoint_policies":policies,
            "instrument_requests":observation.requests,"instrument_local_sd_c":instrument_sd,
            "instrument_ensemble_sd_c":[math.sqrt(sum((m['instrument_temperature_c'][i]-sum(r['instrument_temperature_c'][i] for r in members)/samples)**2 for m in members)/(samples-1)) for i in range(3)],
            "limits":"one closed column; zero Manning/Coriolis; no sources, horizontal circulation, kinetics or SWMM gradients; initial temperatures fixed by sigma index; local derivatives refuse branch boundaries",
            "uncertainty_scope":"Illustrative parameter uncertainty; observation, structural, mesh and forcing-model discrepancy excluded; not DeGray field validation",
            "wall_seconds":time.perf_counter()-started}
        save(output/"column_qualification.json",report)
        fields={"temperature":mean,"parameter_sd":sd,"ensemble_lower_95":lower,"ensemble_upper_95":upper,
                "local_lower_95":[m-1.96*s for m,s in zip(mean,sd)],"local_upper_95":[m+1.96*s for m,s in zip(mean,sd)],
                "wind_response":field(model.jvp(parameters={"wind_speed_m_s":.3}),n),
                "depth_response":field(model.jvp(parameters={"depth_scale":.025}),n),
                "shear_u":field(model.states,n,1),"shear_v":field(model.states,n,2)}
        geometry={"schema":"hc.recorded-mesh/1.0","node_x":[0,100,100,0],"node_y":[0,0,100,100],"node_z":[-10]*4,
                  "face_offsets":[0,4],"face_nodes":[0,1,2,3],"edge_nodes":[0,1,1,2,2,3,3,0],
                  "reference":{"authority":"","id":0,"text":"","units":0,"vertical_authority":"","vertical_id":0,"vertical_text":"","vertical_units":0}}
        frames=[{"time_s":t,"interfaces":[-10*k/n for k in range(n+1)]} for t in model.times]
        return {"id":"column","caption":"Controlled thermal column — analytic sensitivities","columns":1,"layers":n,"frames":frames,"geometry":geometry,"fields":fields,"report":report}


def bellinge_review(library,source,output,samples,seed):
    original=json.loads((source/"replay_inflows.json").read_text());geometry=json.loads((source/"configuration.json").read_text())["receiver_geometry"]
    mesh=dict(geometry,units="m",elevation_kind="bed-elevation",vertical_datum="synthetic flat bed")
    def forward(delta=0,heat=100,wind=.01):
        with open_mesh(library,{"mesh":mesh},source) as horizontal:
            with Simulation(horizontal,layers=3,temperature=10,options={"advectionScheme":0,"maxTimeStep":300,"surfaceHeatFlux":heat,"windStressX":wind}) as sim:
                initial=sim.mass("temperature")*1000*4186;accepted=0;frames=[sim.snapshot("temperature")]
                for original_frame in original["frames"]:
                    frame=copy.deepcopy(original_frame)
                    for data in frame["inflows"]:data["heat_j"]+=1000*4186*data["volume_m3"]*delta
                    accepted+=sim.apply(frame)["acceptedHeatJ"];frames.append(sim.snapshot("temperature"))
                residual=sim.mass("temperature")*1000*4186-initial-accepted-heat*sim.area*sim.time
                if abs(residual)>2:raise RuntimeError("Bellinge heat conservation gate failed")
                return frames,residual
    started=time.perf_counter();baseline,residual=forward();mean=flat([f["values"] for f in baseline]);directions=[(1.,0,0),(0,20.,0),(0,0,.002)]
    reference=[];sensitivities=[]
    for d in directions:
        sweeps=[];previous=None
        for epsilon in (.01,.001,.0001):
            plus,_=forward(epsilon*d[0],100+epsilon*d[1],.01+epsilon*d[2]);minus,_=forward(-epsilon*d[0],100-epsilon*d[1],.01-epsilon*d[2])
            derivative=[(a-b)/(2*epsilon) for a,b in zip(flat([f["values"] for f in plus]),flat([f["values"] for f in minus]))]
            sweeps.append({"epsilon":epsilon,"max_change_from_previous_c":None if previous is None else max(abs(a-b) for a,b in zip(derivative,previous))});previous=derivative
        sensitivities.append(previous);reference.append({"direction":d,"steps":sweeps})
    rng=random.Random(seed);draws=[(rng.gauss(0,1),rng.gauss(0,1),rng.gauss(0,1)) for _ in range(samples)]
    def member(pair):
        index,z=pair
        try:
            frames,residual=forward(z[0],100+20*z[1],.01+.002*z[2])
            return {"index":index,"z":z,"temperature":flat([f["values"] for f in frames]),"heat_residual_j":residual,"status":"completed"}
        except Exception as error:return {"index":index,"z":z,"status":"failed","error":str(error)}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:members=list(executor.map(member,enumerate(draws)))
    save(output/"bellinge_ensemble.json",members)
    if any(m["status"]!="completed" for m in members):raise RuntimeError("Bellinge ensemble failures retained; review incomplete")
    lower,upper=bands([m["temperature"] for m in members]);sd=[math.sqrt(sum(col[i]**2 for col in sensitivities)) for i in range(len(mean))]
    report={"status":"full native ensemble with finite-difference sensitivity references","samples":samples,"successful":samples,"failed":0,"seed":seed,
        "assumptions":"Synthetic flat receiver, initial 10 C; inlet temperature offset SD 1 C; prescribed surface heat mean 100 W/m2 SD 20; x wind stress mean .01 Pa SD .002; independent Gaussian parameters",
        "frozen":"recorded Bellinge water discharge, mapping, synthetic receiver geometry and solver choices",
        "maximum_heat_residual_j":max(abs(m["heat_residual_j"]) for m in members),"finite_difference_reference":reference,
        "derivative_scope":"No analytic full Bellinge thermal adjoint; reference differences and whole-run ensemble retain horizontal circulation/thermal feedback",
        "uncertainty_scope":"Parameter-only marginal intervals; observations, structural discrepancy and geometry uncertainty excluded. Not observed Bellinge receiving-water validation.",
        "wall_seconds":time.perf_counter()-started}
    save(output/"bellinge_qualification.json",report)
    old=json.loads((source/"live_output/run.json").read_text());item=next(r for r in old["results"] if r["component"]=="reservoir" and r["item"]=="temperature")
    geom=json.loads((source/"live_output"/item["geometry"]["artifact"]).read_text())
    return {"id":"bellinge","caption":"Bellinge — full-run thermal parameter uncertainty","columns":3,"layers":3,"frames":baseline,"geometry":geom,
        "fields":{"temperature":mean,"ensemble_lower_95":lower,"ensemble_upper_95":upper,"reference_sd":sd,
                  "inlet_temperature_response":sensitivities[0],"surface_heat_response":sensitivities[1],"wind_stress_response":sensitivities[2]},"report":report}


def gui_case(root,case,epoch):
    directory=root/case["id"];directory.mkdir();save(directory/"geometry.json",case["geometry"])
    frames=case["frames"];n,k=case["columns"],case["layers"];count=n*k
    with (directory/"fields.csv").open("w") as stream:
        w=csv.writer(stream);w.writerow(("time","item","index","value"))
        for t,frame in enumerate(frames):
            clock=epoch+frame["time_s"]/86400
            for item,values in case["fields"].items():
                for i in range(count):w.writerow((clock,item,i,values[t*count+i]))
            for i,z in enumerate(frame["interfaces"]):w.writerow((clock,"interfaces",i,z))
    clock={"count":len(frames),"start":epoch,"end":epoch+frames[-1]["time_s"]/86400,"units":"julian_day"}
    def entry(item,shape,dimensions,units):return {"component":"reservoir","item":item,"artifact":"fields.csv","format":"csv","variable":item,"kind":"Float64","shape":[len(frames)]+shape,"dimensions":["time"]+dimensions,"time":clock,"units":units}
    vertical=entry("interfaces",[n,k+1],["columns","interfaces"],"m");results=[]
    for item in case["fields"]:
        shear=item.startswith("shear_");r=entry(item,[n,k],["columns","layers"],"m/s" if shear else "degC")
        r.update(location="face",geometry={"artifact":"geometry.json","entity_axis":1,"layer_axis":2,"vertical_kind":1,"vertical":vertical},
            quantity={"scale":1,"offset":0 if shear or item!="temperature" else 273.15,"value_kind":1,"powers":[1,0,-1,0,0,0,0,0,0,0] if shear else [0,0,0,0,1,0,0,0,0,0]},
            description=case["report"]["uncertainty_scope"])
        if item in ("ensemble_lower_95","ensemble_upper_95","local_lower_95","local_upper_95"):r["quantity"]["offset"]=273.15
        results.append(r)
    def manifest(path,items,caption,derived=None):
        run={"id":str(uuid.uuid4()),"caption":caption,"status":"completed"}
        save(path,{"schema_version":"1.0","run":run,"components":[{"id":"reservoir","caption":caption,"final_status":"Finished"}],"results":items,"provenance":case["report"] if derived is None else {"derived_from":derived,"scope":"Marginal temperature envelope, not an independent trajectory"}})
        return run
    path=directory/"run.json";run=manifest(path,results,case["caption"]);requests=[]
    for item,label,color in [("temperature","Baseline","#0072b2"),("ensemble_lower_95","Ensemble lower 95%","#009e73"),("ensemble_upper_95","Ensemble upper 95%","#d55e00")]:
        view=path;viewrun=run
        if item!="temperature":
            view=directory/(item+".json");r=copy.deepcopy(results[0]);r["variable"]=item
            viewrun=manifest(view,[r],case["caption"]+" — "+label,str(path))
        requests.append({"type":"recorded-plot","manifest":str(view),"component":"reservoir","item":"temperature","name":label,"units":"degC","feature":0,"verticalLayer":1 if k>1 else 0,"layerCount":k,"legendLabel":label,"plotName":label,"plotColor":color,"plotPoints":True,"plotVisible":True,
            "runId":viewrun["id"],"runTitle":viewrun["caption"],"runStatus":"completed","unitSignature":{"caption":"degC","scale":1,"offset":273.15,"powers":[0,0,0,0,1,0,0,0,0,0]}})
    layer={"type":"recorded-result","manifest":str(path),"component":"reservoir","item":"temperature","name":case["caption"],"units":"degC","opacity":1,"visible":True,"verticalLayer":1 if k>1 else 0,"scene":{"shownIn3D":True,"terrainEnabled":True,"zMode":2,"offset":0},"runId":run["id"],"runTitle":run["caption"],"runStatus":"completed"}
    return path,layer,requests


def figures(root,cases):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    for row,case in enumerate(cases):
        n,k=case["columns"],case["layers"];shape=(len(case["frames"]),n,k);f={key:np.array(values).reshape(shape) for key,values in case["fields"].items()}
        time=np.array([s["time_s"] for s in case["frames"]])/3600;ax=axes[row,0]
        ax.fill_between(time,f["ensemble_lower_95"][:,0,1],f["ensemble_upper_95"][:,0,1],color="#0072b2",alpha=.2,label="Full-run sample 2.5–97.5%")
        ax.plot(time,f["temperature"][:,0,1],color="#0072b2",label="Baseline")
        if "local_lower_95" in f:
            ax.plot(time,f["local_lower_95"][:,0,1],"--",color="#d55e00",label="Local Gaussian parameter 95%")
            ax.plot(time,f["local_upper_95"][:,0,1],"--",color="#d55e00")
        ax.set(title=case["caption"],xlabel="Elapsed hours",ylabel="Layer 1 temperature (°C)");ax.legend(fontsize=8)
        ax=axes[row,1]
        names=[key for key in f if key.endswith("_response")]
        for key,color,style in zip(names,["#0072b2","#d55e00","#009e73"],["-","--","-."]):ax.plot(time,f[key][:,0,1],color=color,linestyle=style,label=key.replace("_"," "))
        ax.set(title="Responses to stated parameter changes",xlabel="Elapsed hours",ylabel="Temperature response (°C)");ax.legend(fontsize=8)
    fig.suptitle("Combined thermal, geometry and uncertainty review\nSynthetic receiver/parameters; Bellinge water forcing recorded. Observation and structural discrepancy excluded.",fontsize=12)
    fig.savefig(root/"thermal_review.png",dpi=160);fig.savefig(root/"thermal_review.pdf");plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("source",type=Path);parser.add_argument("output",type=Path)
    parser.add_argument("--samples",type=int,default=256);parser.add_argument("--library");parser.add_argument("--figures",action="store_true")
    args=parser.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False);library=Library(args.library)
    if args.samples<32:raise ValueError("At least 32 samples required")
    column=column_review(library,root,args.samples,61015);bellinge=bellinge_review(library,args.source.resolve(),root,args.samples,61016)
    save(root/"review_data.json",{"cases":[column,bellinge]})
    epoch=json.loads((args.source/"replay_inflows.json").read_text())["epoch_julian_day"]
    views=[gui_case(root,c,epoch) for c in (bellinge,column)];path,layer,requests=views[0]
    composition={"schema_version":"1.1","metadata":{"id":"combined-thermal-review","caption":"Combined thermal and uncertainty review"},"components":[{"id":"reservoir","caption":bellinge["caption"],"execution":{"mode":"open","results":str(path)}}]}
    save(root/"Thermal_review.composition.json",composition)
    catalog=[r["manifest"] for _,_,rs in views for r in rs]
    # Existing observed-profile investigations remain separate reference runs.
    for relative in ("qualified_seasonal/baseline/gui_run.json","original_seasonal/baseline/gui_run.json"):
        reference=args.source.resolve().parents[1]/relative
        if reference.is_file():catalog.append(str(reference))
    save(root/"Thermal_review.composition.composer.json",{"sidecar_version":1,"layers":[layer],"results":{"version":1,"runs":catalog,"alignment":0,"active":requests[0],"plot":{"version":1,"pinned":True,"requests":requests}},"navigation":{"version":1,"workspace":"results"}})
    save(root/"Thermal_native.composition.json",composition);save(root/"Thermal_native.composition.composer.json",{"sidecar_version":1,"layers":[layer],"results":{"version":1,"runs":[str(path)],"alignment":0}})
    report={"column":column["report"],"bellinge":bellinge["report"],"provenance":{"library_sha256":sha(library.path),"driver_sha256":sha(Path(__file__)),"recorded_water_heat_sha256":sha(args.source/"replay_inflows.json")}}
    save(root/"qualification.json",report)
    if args.figures:figures(root,(column,bellinge))
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()
