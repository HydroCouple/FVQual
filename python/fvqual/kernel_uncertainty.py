"""D1 kernel JVP/VJP gates and local variance compared with Monte Carlo."""
import argparse
import json
import math
import random
import statistics
from pathlib import Path
from .engine import Library

def experiment(library,samples=20000,seed=20261010):
    rng=random.Random(seed)
    cases=[
        ("decay",0,[2,1e-6,3600],[.02,1e-8,0]),
        ("column_heat",1,[10,200,86400,10,1000,4186],[.05,10,0,.05,0,0]),
        ("mix",2,[1000,10,.1,1.2],[.5,.02,.005,.02])]
    report={"scope":"smooth reference kernels only; independent Gaussian input uncertainty",
            "samples":samples,"seed":seed,"kernels":[]}
    for name,kind,x,sd in cases:
        gradient=library.kernel(kind,x,cotangent=1)
        tangent=[s*rng.gauss(0,1) for s in sd]
        jv=library.kernel(kind,x,direction=tangent)
        w=.73
        pullback=library.kernel(kind,x,cotangent=w)
        adjoint_error=abs(w*jv-sum(a*b for a,b in zip(pullback,tangent)))
        finite_differences=[]
        for h in (1e-2,1e-3,1e-4):
            plus=library.kernel(kind,[a+h*b for a,b in zip(x,tangent)])
            minus=library.kernel(kind,[a-h*b for a,b in zip(x,tangent)])
            finite_differences.append({"epsilon":h,"central":(plus-minus)/(2*h),"jvp":jv})
        predicted=sum((g*s)**2 for g,s in zip(gradient,sd))
        values=[library.kernel(kind,[a+s*rng.gauss(0,1) for a,s in zip(x,sd)]) for _ in range(samples)]
        actual=statistics.variance(values)
        report["kernels"].append({"name":name,"mean_inputs":x,"input_sd":sd,
            "value":library.kernel(kind,x),"vjp":gradient,"adjoint_error":adjoint_error,
            "finite_difference_sweep":finite_differences,"local_variance":predicted,
            "monte_carlo_variance":actual,"variance_ratio":actual/predicted})
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("output",type=Path)
    p.add_argument("--library"); p.add_argument("--samples",type=int,default=20000)
    a=p.parse_args()
    if a.samples<2: raise ValueError("At least two samples required")
    result=experiment(Library(a.library),a.samples)
    with a.output.open("x") as f: f.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k["name"]:k["variance_ratio"] for k in result["kernels"]},indent=2))
if __name__=="__main__": main()
