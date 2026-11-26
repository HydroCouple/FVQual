"""Verified live SWMM surface boundary to a mixed-mesh FVQual reservoir."""
import argparse
import copy
import json
from pathlib import Path
from .engine import Library, Simulation
from .recorded import open_mesh, read_case, forward, fingerprint, manifest
from .swmm import SurfaceBoundaryDonor

def run(library,swmm_library,case,model,output):
    path,spec=read_case(case)
    output=Path(output); output.mkdir(parents=True,exist_ok=False)
    mapping={"TSS":{"target":"tracer","kg_per_native_mass":.001}}
    with open_mesh(library,spec,path.parent) as mesh:
        with Simulation(mesh,layers=spec["layers"],temperature=spec["temperature_c"],
                        species=spec["species"],options=spec.get("options")) as reservoir:
            with SurfaceBoundaryDonor(swmm_library,model,output/"donor",0,0,mapping) as donor:
                cp=reservoir.options.referenceDensity*reservoir.options.heatCapacity
                water0=donor.volume+reservoir.volume
                mass0=donor.species_inventory("TSS")*.001+reservoir.mass("tracer")
                heat0=(donor.species_inventory("__TEMPERATURE__")+reservoir.mass("temperature"))*cp
                donor.prescribe_outflow(.01)
                frames=[]; receipts=[]
                for target in (1.375,2.875,4.125,5.5):
                    frame=donor.advance(target)
                    receipt=reservoir.apply(frame)
                    donor.accept(receipt["id"])
                    frames.append(frame); receipts.append(receipt)
                residuals={"water_m3":donor.volume+reservoir.volume-water0,
                           "species_kg":donor.species_inventory("TSS")*.001+reservoir.mass("tracer")-mass0,
                           "heat_j":(donor.species_inventory("__TEMPERATURE__")+reservoir.mass("temperature"))*cp-heat0}
                for name,total in (("water_m3",water0),("species_kg",mass0),("heat_j",heat0)):
                    if abs(residuals[name])>1e-10*max(1,abs(total)):
                        raise RuntimeError("Combined budget failed: "+name)
                live={"volume_m3":reservoir.volume,"mean_temperature_c":reservoir.mean_temperature,
                      "tracer_kg":reservoir.mass("tracer")}
    recorded=copy.deepcopy(spec); recorded["frames"]=frames
    recorded["source"]="Accepted native SWMM 2D boundary exports, already debited by donor"
    recorded["donor_model_sha256"]=fingerprint(model)
    recorded["donor_library_sha256"]=fingerprint(swmm_library)
    (output/"recorded_event.json").write_text(json.dumps(recorded,indent=2)+"\n")
    replay=forward(library,recorded,path.parent)
    differences={"water_m3":live["volume_m3"]-replay["final"]["volume_m3"],
                 "temperature_c":live["mean_temperature_c"]-replay["final"]["mean_temperature_c"],
                 "species_kg":live["tracer_kg"]-replay["final"]["species_mass_kg"]["tracer"]}
    if any(abs(v)>1e-12 for v in differences.values()):
        raise RuntimeError("Recorded replay differs from live coupling")
    report={"manifest":manifest(library,path,spec),"swmm_library_sha256":fingerprint(swmm_library),
            "swmm_model_sha256":fingerprint(model),"combined_budget_residuals":residuals,
            "live_recorded_difference":differences,"receipts":receipts,"reservoir_final":live,
            "scope":"one-way SWMM surface boundary; native 1D outfall component and stage feedback remain separate gates"}
    (output/"coupling.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("case",type=Path); p.add_argument("model",type=Path); p.add_argument("output",type=Path)
    p.add_argument("--swmm-library",required=True); p.add_argument("--library")
    a=p.parse_args(); result=run(Library(a.library),a.swmm_library,a.case,a.model,a.output)
    print(json.dumps({"combined_budget_residuals":result["combined_budget_residuals"],
                      "live_recorded_difference":result["live_recorded_difference"]},indent=2))
if __name__=="__main__": main()
