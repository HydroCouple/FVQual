"""Independent-prior ensemble and whole-run finite-difference reference."""
import argparse
import json
import math
import random
import statistics
from pathlib import Path
from .engine import Library
from .recorded import forward, read_case, manifest

def observations(result):
    return {"mean_temperature_c":result["final"]["mean_temperature_c"],
            "volume_m3":result["final"]["volume_m3"],
            **{f"mass_kg:{k}":v for k,v in result["final"]["species_mass_kg"].items()}}

def quantile(values,p):
    values=sorted(values)
    position=p*(len(values)-1)
    low=int(position)
    return values[low]+(values[min(low+1,len(values)-1)]-values[low])*(position-low)

def sensitivity(library,spec,base,parameter,nominal,steps=(1e-2,1e-3,1e-4)):
    results=[]
    for h in steps:
        plus=observations(forward(library,spec,base,{parameter:nominal+h}))
        minus=observations(forward(library,spec,base,{parameter:nominal-h}))
        results.append({"step":h,"derivative":{k:(plus[k]-minus[k])/(2*h) for k in plus}})
    return results

def ensemble(library,path,spec,output,samples=64,seed=20261010,flow_cv=0.1,
             load_cv=0.1,temperature_sd=0.5):
    if samples<2 or min(flow_cv,load_cv,temperature_sd)<0: raise ValueError("Invalid ensemble settings")
    output=Path(output)
    output.mkdir(parents=True,exist_ok=False)
    rng=random.Random(seed)
    rows=[]; failures=[]
    def lognormal(cv):
        sigma=math.sqrt(math.log1p(cv*cv))
        return math.exp(rng.gauss(-sigma*sigma/2,sigma))
    for i in range(samples):
        parameters={"flow_scale":lognormal(flow_cv),"load_scale":lognormal(load_cv),
                    "temperature_delta_c":rng.gauss(0,temperature_sd)}
        directory=output/f"sample_{i:05d}"
        directory.mkdir()
        try:
            result=forward(library,spec,path.parent,parameters)
            (directory/"result.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
            rows.append({"sample":i,"parameters":parameters,"observations":observations(result)})
        except (RuntimeError,ValueError) as error:
            failure={"sample":i,"parameters":parameters,"error":str(error)}
            failures.append(failure)
            (directory/"failure.json").write_text(json.dumps(failure,indent=2)+"\n")
    summary={}
    if rows:
        for name in rows[0]["observations"]:
            values=[r["observations"][name] for r in rows]
            summary[name]={"mean":statistics.mean(values),"sd":statistics.stdev(values) if len(values)>1 else 0,
                           "p025":quantile(values,.025),"p50":quantile(values,.5),"p975":quantile(values,.975)}
    fd={name:sensitivity(library,spec,path.parent,name,nominal)
        for name,nominal in (("flow_scale",1),("load_scale",1),("temperature_delta_c",0))}
    report={"manifest":manifest(library,path,spec),"seed":seed,"samples":samples,
            "successful":len(rows),"failed":len(failures),"failures":failures,
            "priors":{"dependence":"independent; correlated priors are not implemented",
                      "flow_scale":{"distribution":"mean-one lognormal","cv":flow_cv},
                      "load_scale":{"distribution":"mean-one lognormal","cv":load_cv},
                      "temperature_delta_c":{"distribution":"normal","mean":0,"sd":temperature_sd}},
            "prediction_bands":summary,"finite_difference_sweep":fd,
            "interpretation":"conditional on supplied SWMM forcing; no SWMM parameter gradients or full solver VJP",
            "convergence":"sampling bands are exploratory; increase sample size and compare independent seeds"}
    (output/"ensemble.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("case",type=Path); p.add_argument("output",type=Path)
    p.add_argument("--samples",type=int,default=64); p.add_argument("--seed",type=int,default=20261010)
    p.add_argument("--library")
    args=p.parse_args(); path,spec=read_case(args.case)
    result=ensemble(Library(args.library),path,spec,args.output,args.samples,args.seed)
    print(json.dumps({k:result[k] for k in ("successful","failed","prediction_bands")},indent=2))
if __name__=="__main__": main()
