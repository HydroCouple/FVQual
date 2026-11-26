"""Reproducible runs with integrated donor exports, receipts and closed budgets."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time
from .engine import Library, Mesh, Simulation

SCHEMA="fvqual.recorded-inflows/1.0"
UNITS={"time":"elapsed_s","water":"m3","heat":"J_relative_to_0_degC","species":"kg"}

def read_case(path):
    path=Path(path).resolve()
    spec=json.loads(path.read_text())
    if spec.get("schema")!=SCHEMA or spec.get("units")!=UNITS:
        raise ValueError("Unsupported recorded forcing schema or unit convention")
    previous=0.0
    ids=set()
    for frame in spec["frames"]:
        if frame["id"] in ids or frame["start_s"]!=previous or not frame["end_s"]>previous:
            raise ValueError("Frames need unique IDs and contiguous increasing intervals from zero")
        ids.add(frame["id"])
        previous=frame["end_s"]
    if not spec["frames"]: raise ValueError("No recorded intervals")
    return path,spec

def open_mesh(library, spec, base):
    mesh=spec["mesh"]
    if "package" in mesh: return Mesh.read(library,base/mesh["package"])
    if mesh.get("units")!="m" or mesh.get("elevation_kind")!="bed-elevation" or not mesh.get("vertical_datum"):
        raise ValueError("Synthetic mesh needs metre units and explicit bed elevation datum")
    return Mesh.from_arrays(library,mesh["x"],mesh["y"],mesh["z"],mesh["offsets"],mesh["nodes"])

def transformed_frames(spec, parameters):
    scale=parameters.get("flow_scale",1.0)
    load=parameters.get("load_scale",1.0)
    delta=parameters.get("temperature_delta_c",0.0)
    if not math.isfinite(scale) or not math.isfinite(load) or scale<0 or load<0 or not math.isfinite(delta):
        raise ValueError("Invalid inflow uncertainty parameters")
    rho_cp=spec.get("options",{}).get("referenceDensity",1000)*spec.get("options",{}).get("heatCapacity",4186)
    frames=copy.deepcopy(spec["frames"])
    for frame in frames:
        for data in frame["inflows"]:
            data["volume_m3"]*=scale
            data["heat_j"]=data["heat_j"]*scale+rho_cp*data["volume_m3"]*delta
            data["species_mass_kg"]={k:v*scale*load for k,v in data.get("species_mass_kg",{}).items()}
    return frames

def forward(library,spec,base,parameters=None):
    parameters=parameters or {}
    allowed={"flow_scale","load_scale","temperature_delta_c","surface_heat_flux_w_m2"}
    if set(parameters)-allowed: raise ValueError("Unknown uncertain parameter")
    options=dict(spec.get("options",{}))
    if "surface_heat_flux_w_m2" in parameters:
        options["surfaceHeatFlux"]=parameters["surface_heat_flux_w_m2"]
    start=time.perf_counter()
    with open_mesh(library,spec,base) as mesh:
        with Simulation(mesh,layers=spec.get("layers",4),temperature=spec.get("temperature_c",10),
                        species=spec.get("species",{}),options=options) as simulation:
            initial_volume=simulation.volume
            initial_heat=simulation.mass("temperature")*simulation.options.referenceDensity*simulation.options.heatCapacity
            initial_mass={k:simulation.mass(k) for k in simulation.species}
            receipts=[]
            water=0.0; heat=0.0; masses={k:0.0 for k in simulation.species}
            history=[]
            for frame in transformed_frames(spec,parameters):
                receipt=simulation.apply(frame)
                receipts.append(receipt)
                water+=receipt["acceptedVolumeM3"]
                heat+=receipt["acceptedHeatJ"]
                for k,v in receipt["species_mass_kg"].items(): masses[k]+=v
                history.append({"elapsed_s":simulation.time,"volume_m3":simulation.volume,
                                "mean_temperature_c":simulation.mean_temperature,
                                "species_mass_kg":{k:simulation.mass(k) for k in simulation.species}})
            final=history[-1]
            final_heat=simulation.mass("temperature")*simulation.options.referenceDensity*simulation.options.heatCapacity
            heat+=simulation.options.surfaceHeatFlux*simulation.area*simulation.time
            budgets={"water_m3":final["volume_m3"]-initial_volume-water,
                     "heat_j":final_heat-initial_heat-heat,
                     "species_kg":{k:simulation.mass(k)-initial_mass[k]-masses[k] for k in masses}}
            return {"parameters":parameters,"initial_volume_m3":initial_volume,"final":final,
                    "receipts":receipts,"history":history,"budget_residuals":budgets,
                    "wall_seconds":time.perf_counter()-start,
                    "coverage":"CPU, always wet, fixed topology; one-way integrated donor exports"}

def fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def manifest(library,path,spec):
    repository=Path(__file__).resolve().parents[2]
    def git(*args):
        result=subprocess.run(["git","-C",str(repository),*args],capture_output=True,text=True,check=False)
        return result.stdout.strip() if result.returncode==0 else None
    value={"schema":"fvqual.run-manifest/1.0","case_sha256":fingerprint(path),
           "library":str(library.path),"library_sha256":fingerprint(library.path),
           "fvqual_version":library.version,"revision":git("rev-parse","HEAD"),
           "working_changes":git("status","--short"),"units":UNITS,
           "epoch_utc":spec.get("epoch_utc"),"solver_options":spec.get("options",{}),
           "derivative_coverage":"smooth reference kernels only; whole-run sensitivities use finite differences"}
    if "package" in spec["mesh"]:
        package=path.parent/spec["mesh"]["package"]
        value["mesh_manifest_sha256"]=fingerprint(package)
        package_data=json.loads(package.read_text())
        if "mesh" in package_data: value["mesh_geometry_sha256"]=fingerprint(package.parent/package_data["mesh"])
    return value

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--library")
    args=parser.parse_args()
    path,spec=read_case(args.case)
    args.output.mkdir(parents=True,exist_ok=False)
    library=Library(args.library)
    result=forward(library,spec,path.parent)
    result["manifest"]=manifest(library,path,spec)
    (args.output/"result.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"final":result["final"],"budget_residuals":result["budget_residuals"]},indent=2))

if __name__=="__main__": main()
