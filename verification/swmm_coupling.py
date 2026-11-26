#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Qualify the native SWMM component through Composer and its public engine API.

All artifacts are written under tests/output by default. Bellinge uses its
original runoff/network and an explicitly synthetic receiving basin. No plan
or source model is modified. No third-party Python packages are required.
"""
import argparse
import copy
import csv
import ctypes as C
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def sections(path):
    result, section = {}, ""
    for line in path.read_text().splitlines():
        line = line.split(";", 1)[0].strip()
        if line.startswith("["):
            section = line.upper()
            result.setdefault(section, [])
        elif line:
            result.setdefault(section, []).append(line.split())
    return result


def generate_bellinge(args):
    source = args.bellinge.resolve()
    blocks = sections(source)
    coordinates = {row[0]: (float(row[1]), float(row[2])) for row in blocks["[COORDINATES]"]}
    outfalls = sorted(row[0] for row in blocks["[OUTFALLS]"])
    folder = args.output / "bellinge"
    folder.mkdir(parents=True, exist_ok=True)
    # Rebase file references only. IGNORE_2D is applied through the public API
    # by the wrapper before initialization; the source network stays intact.
    text = source.read_text()
    text = re.sub(r'FILE\s+"([^"\n]+)"', lambda m: 'FILE "' + str((source.parent / m[1]).resolve()) + '"', text)
    text = re.sub(r'(?m)^(\s*(?:2D_)?MESH_FILE\s+)([^;\n]+)', lambda m: m[1] + '"' + str((source.parent / m[2].strip().strip('"')).resolve()) + '"', text)
    section, rebased = "", []
    for line in text.splitlines():
        if line.strip().startswith("["): section=line.strip().upper()
        if section=="[2D_MESH_FILE]" and line.strip().upper().startswith("FILE "):
            reference=line.strip().split(None,1)[1].strip().strip('"')
            line='FILE "'+str((source.parent/reference).resolve())+'"'
        rebased.append(line)
    (folder / "network.inp").write_text("\n".join(rebased)+"\n")
    options = {row[0]: row[1] for row in blocks["[OPTIONS]"]}
    start = dt.datetime.strptime(options["START_DATE"] + " " + options["START_TIME"], "%m/%d/%Y %H:%M:%S")
    epoch = 2440587.5 + (start - dt.datetime(1970, 1, 1)).total_seconds() / 86400
    # Location is for review, not surveyed bathymetry or an inferred destination.
    x = round(max(v[0] for v in coordinates.values()) + 100, 2)
    y = round(sum(coordinates[o][1] for o in outfalls) / len(outfalls), 2)
    xy = [(x, y), (x+50, y), (x+50, y+100), (x, y+100), (x+100, y), (x+100, y+100)]
    mesh = ";; UNITS: SI (m)\n;; SOURCE_CRS: EPSG:25832\n[2D_VERTICES]\n"
    mesh += "".join(f"{a:.2f} {b:.2f} -5\n" for a, b in xy)
    mesh += "[2D_TRIANGLES]\n1 4 2 .03 5\n4 5 2 .03 5\n[2D_QUADS]\n0 1 2 3 .03 5\n"
    (folder / "receiver.2dm").write_text(mesh)
    package = folder / "receiver.mesh.json"
    if not package.exists():
        subprocess.run([str(args.fvqual), "import-mesh", str(folder / "receiver.2dm"), str(package),
                        "--crs", "EPSG:25832", "--datum", "synthetic-review-bed-minus-5m",
                        "--elevation-kind", "bed-elevation"], check=True, capture_output=True, text=True)
    geometry = dict(x=[a for a, b in xy], y=[b for a, b in xy], z=[-5.]*6,
                    offsets=[0, 3, 6, 10], nodes=[1, 4, 2, 4, 5, 2, 0, 1, 2, 3])
    faces = json.loads(package.read_text())["face_ids"]
    mapping = []
    for i, name in enumerate(outfalls):
        c = i % 3
        mapping.append(dict(id=name, temperature={"mode": "fixed", "value_degC": 18,
            "provenance": "Synthetic 18 Celsius inflow; Bellinge has no inlet temperature observations in this model"},
            placement="surface", targets=[dict(column=c, face_id=faces[c], weight=.75),
                                          dict(column=(c+1)%3, face_id=faces[(c+1)%3], weight=.25)]))
    config = dict(schema="openswmm.reservoir-coupling/1.0", instance_id="openswmm", model="network.inp",
                  role="network_only", ownership="swmm-network/fvqual-receiver", epoch_julian_day=epoch,
                  exchange_times_s=[7.,15.,31.,60.,120.,300.,600.,1800.,3600.,7200.,10800.,14400.,15000.,15600.,16200.,16800.,17400.,18000.,19800.,21600.], receiver_columns=3,
                  receiver_geometry=geometry, species={}, outfalls=mapping, output_directory="donor_output")
    write(folder / "configuration.json", config)
    composition = json.loads((ROOT / "examples/swmm_reservoir/live/composition.json").read_text())
    composition["metadata"] = {"id": "bellinge-live-coupling", "caption": "Bellinge SWMM network → synthetic mixed-mesh receiving basin"}
    receiver = composition["components"][1]["arguments"]
    receiver["mesh"] = {"@uri": "receiver.mesh.json"}
    receiver["geometry_values"]["values"] = [100, 100, 5, 10, 7, 21600]
    receiver["epoch_julian_day"] = {"values": [epoch]}
    receiver["conservative_species"] = {}
    receiver["maximum_time_step_s"] = {"values": [2.]}
    write(folder / "composition.json", composition)
    write(folder / "composition.composer.json", {"sidecar_version": 1,
          "components": {"openswmm": {"x": -300, "y": -70}, "reservoir": {"x": 180, "y": -70}}})
    with (folder / "mapping.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["outfall", "source_x_m", "source_y_m", "receiver_column", "face_id", "weight", "temperature_degC", "basis"])
        for row in mapping:
            for target in row["targets"]:
                writer.writerow([row["id"], *coordinates[row["id"]], target["column"], target["face_id"], target["weight"], 18, "synthetic review mapping and temperature"])
    # Keep these layers in the model's projected CRS. Mapping lines express
    # ownership/weights for software review, not physical pipe alignments.
    def geojson(name, features):
        write(folder / name, {"type": "FeatureCollection", "name": name.removesuffix(".geojson"),
            "crs": {"type": "name", "properties": {"name": "EPSG:25832"}}, "features": features})

    def feature(kind, points, properties):
        return {"type": "Feature", "properties": properties,
                "geometry": {"type": kind, "coordinates": points}}

    vertices = {}
    for row in blocks.get("[VERTICES]", []):
        vertices.setdefault(row[0], []).append([float(row[1]), float(row[2])])
    network = []
    for kind in ("CONDUITS", "PUMPS", "WEIRS", "ORIFICES", "OUTLETS"):
        for row in blocks.get("[" + kind + "]", []):
            if row[1] in coordinates and row[2] in coordinates:
                network.append(feature("LineString", [coordinates[row[1]], *vertices.get(row[0], []), coordinates[row[2]]],
                    {"id": row[0], "kind": kind, "from_node": row[1], "to_node": row[2], "basis": "source Bellinge network"}))
    geojson("bellinge_network.geojson", network)
    geojson("bellinge_outfalls.geojson", [feature("Point", coordinates[name],
        {"outfall": name, "basis": "source Bellinge outfall"}) for name in outfalls])
    centroids = []
    for a, b in zip(geometry["offsets"], geometry["offsets"][1:]):
        nodes = geometry["nodes"][a:b]
        centroids.append([sum(xy[n][0] for n in nodes)/len(nodes), sum(xy[n][1] for n in nodes)/len(nodes)])
    geojson("synthetic_transfer_mapping.geojson", [feature("LineString", [coordinates[row["id"]], centroids[target["column"]]],
        {"outfall": row["id"], "receiver_column": target["column"], "face_id": target["face_id"],
         "weight": target["weight"], "temperature_degC": 18, "basis": "synthetic review mapping; not a surveyed destination"})
        for row in mapping for target in row["targets"]])
    write(folder / "provenance.json", {"source_model": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "role": "Native Bellinge runoff and 1D network; native 2D disabled before initialization",
        "receiver": "Synthetic 100 by 100 m basin, flat -5 m bed, two triangles and one quad; not Bellinge bathymetry",
        "mapping": "All nine source outfalls assigned explicitly for software qualification; weights and destination require field review",
        "temperature": "Prescribed synthetic 18 Celsius; native Bellinge heat transport is absent",
        "species": "Bellinge source contains no pollutants; kg species are qualified separately in the synthetic native-heat case",
        "clock": "Original civil model start 2012-06-29 00:01 retained as a shared Julian-day epoch; source timezone unspecified",
        "window": "Six-hour coupling review within the original 48-hour input; native report may label the early stop interrupted",
        "map_layers": ["bellinge_network.geojson", "bellinge_outfalls.geojson", "synthetic_transfer_mapping.geojson"],
        "exchange_times_s": config["exchange_times_s"]})
    return folder / "composition.json"


class Receipt(C.Structure):
    _fields_ = [("exported_m3", C.c_double), ("reversed_m3", C.c_double), ("exported_heat_j", C.c_double), ("flags", C.c_uint)]


class Engine:
    def __init__(self, library):
        self.dll = C.CDLL(str(library))
        p, d, i, s = C.c_void_p, C.c_double, C.c_int, C.c_char_p
        signatures = {"swmm_engine_create": (p, []), "swmm_engine_destroy": (None, [p]),
          "swmm_engine_open": (i, [p, s, s, s, p]), "swmm_options_set": (i, [p, s, s]),
          "swmm_engine_initialize": (i, [p]), "swmm_engine_start": (i, [p, i]),
          "swmm_engine_advance_to": (i, [p, d, C.POINTER(d)]), "swmm_engine_end": (i, [p]),
          "swmm_engine_close": (i, [p]), "swmm_get_last_error_msg": (s, [p]),
          "swmm_pollutant_count": (i, [p]), "swmm_node_index": (i, [p, s]),
          "swmm_node_get_export_receipt": (i, [p, i, C.POINTER(Receipt), C.POINTER(d), i]),
          "swmm_get_routing_total": (i, [p, i, C.POINTER(d)])}
        for name, (result, arguments) in signatures.items():
            fn = getattr(self.dll, name); fn.restype = result; fn.argtypes = arguments

    def check(self, handle, result):
        if result:
            raise RuntimeError(self.dll.swmm_get_last_error_msg(handle).decode() + f" ({result})")

    def run(self, config, folder):
        folder.mkdir(parents=True, exist_ok=True)
        handle = self.dll.swmm_engine_create()
        history, previous = [], {}
        try:
            self.check(handle, self.dll.swmm_engine_open(handle, str(config["model_path"]).encode(), str(folder/"donor.rpt").encode(), str(folder/"donor.out").encode(), None))
            self.check(handle, self.dll.swmm_options_set(handle, b"IGNORE_2D", b"YES"))
            self.check(handle, self.dll.swmm_engine_initialize(handle))
            self.check(handle, self.dll.swmm_engine_start(handle, 1))
            np = self.dll.swmm_pollutant_count(handle)
            for end in config["exchange_times_s"]:
                actual = C.c_double()
                self.check(handle, self.dll.swmm_engine_advance_to(handle, end, C.byref(actual)))
                assert actual.value == end
                exports = {}
                for row in config["outfalls"]:
                    receipt, amounts = Receipt(), (C.c_double * np)()
                    index = self.dll.swmm_node_index(handle, row["id"].encode())
                    self.check(handle, self.dll.swmm_node_get_export_receipt(handle, index, C.byref(receipt), amounts, np))
                    before = previous.get(row["id"], (0, 0, [0]*np))
                    volume = receipt.exported_m3-before[0]
                    heat = receipt.exported_heat_j-before[1]
                    if row["temperature"]["mode"] == "fixed":
                        heat = volume * 1000 * 4186 * row["temperature"]["value_degC"]
                    assert receipt.reversed_m3 == 0
                    exports[row["id"]] = dict(water=volume, heat=heat, species=[amounts[p]-before[2][p] for p in range(np)])
                    previous[row["id"]] = (receipt.exported_m3, receipt.exported_heat_j, list(amounts))
                native = C.c_double()
                self.check(handle, self.dll.swmm_get_routing_total(handle, 6, C.byref(native)))
                # This API returns internal ft3. All source outfalls are mapped
                # in both review fixtures, so their summed receipt equals it.
                close(sum(value[0] for value in previous.values()), native.value * .028316846592, 1e-9)
                history.append(exports)
            self.check(handle, self.dll.swmm_engine_end(handle))
            return history
        finally:
            self.dll.swmm_engine_close(handle)
            self.dll.swmm_engine_destroy(handle)


def close(a, b, tolerance=1e-8):
    assert math.isfinite(a) and math.isfinite(b) and abs(a-b) <= tolerance + 1e-12*abs(b), (a, b)


def run_host(args, path, success=True):
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen", OPENSWMM_2D_BACKEND="cpu")
    command = [str(args.composer), "--components", str(args.fvqual_components), "--components", str(args.swmm_components), "--run", str(path)]
    result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=180)
    path.with_suffix(".log").write_text(result.stdout + result.stderr)
    assert (result.returncode == 0) == success, result.stdout + result.stderr


def rows(path):
    result = {}
    for row in csv.DictReader(path.open()):
        key = (float(row["time"]), row["item"], int(row["index"]))
        assert key not in result, "duplicate snapshot item identity"
        result[key] = float(row["value"])
    return result


def qualify(args, composition_path):
    folder = composition_path.parent
    composition = json.loads(composition_path.read_text())
    config_path = folder / composition["components"][0]["arguments"]["configuration"]["@uri"]
    config = json.loads(config_path.read_text())
    config["model_path"] = (folder / config["model"]).resolve()
    run_host(args, composition_path)
    manifest = json.loads((folder / "live_output/run.json").read_text())
    assert all(c["final_status"] == "Finished" for c in manifest["components"])
    donor_folder = max((folder / config["output_directory"]).iterdir(), key=lambda p: p.stat().st_mtime_ns)
    packets = [json.loads((donor_folder / f"export_{i}.json").read_text()) for i in range(len(config["exchange_times_s"]))]
    for i, packet in enumerate(packets):
        ack = json.loads((donor_folder / f"accepted_{i}.json").read_text())
        assert ack["transfer_id"] == packet["id"] and ack["end_s"] == packet["end_s"]
    native = Engine(args.engine_library).run(config, folder / "direct_native")
    max_mapping_error = 0
    for exports, packet in zip(native, packets):
        for mapping in config["outfalls"]:
            items = [r for r in packet["inflows"] if r["source_object"] == mapping["id"]]
            assert len(items) == len(mapping["targets"])
            expected = exports[mapping["id"]]
            for item, target in zip(items, mapping["targets"]):
                assert item["column"] == target["column"]
                close(item["volume_m3"], expected["water"]*target["weight"], 1e-10)
                close(item["heat_j"], expected["heat"]*target["weight"], 1e-6)
            max_mapping_error = max(max_mapping_error, abs(sum(r["volume_m3"] for r in items)-expected["water"]))
            masses = [sum(r["species_mass_kg"].get(name, 0) for r in items) for name in config["species"].values()]
            for a, b in zip(masses, expected["species"]): close(a, b, 1e-12)
    series = json.loads((donor_folder / "inflow_series.json").read_text())
    # Retain the receiver's original initial conditions. Donor exports supply
    # amounts, not the reservoir's pre-existing constituent concentrations.
    series["species"] = {}
    write(folder / "replay_inflows.json", series)
    replay = copy.deepcopy(composition)
    replay["components"] = [replay["components"][1]]
    replay["components"][0]["arguments"]["recorded_inflows"] = {"@uri": "replay_inflows.json"}
    replay["connections"] = []
    replay["writers"][0]["path"] = "replay_output/fields.csv"
    replay["run"]["manifest"] = "replay_output/run.json"
    write(folder / "replay.composition.json", replay)
    run_host(args, folder / "replay.composition.json")
    live_rows, replay_rows = rows(folder / "live_output/fields.csv"), rows(folder / "replay_output/fields.csv")
    differences = [abs(live_rows[key]-value) for key, value in replay_rows.items()]
    assert max(differences) == 0, max(differences)
    water_residual = heat_residual = 0
    # sigma defaults to uniform thirds; known initial field [22,16,12].
    initial_water = 200*5 if not args.bellinge or folder.name != "bellinge" else 10000*5
    initial_heat = initial_water*(22+16+12)/3*1000*4186
    cumulative_water = cumulative_heat = 0
    for packet in packets:
        day = config["epoch_julian_day"] + packet["end_s"]/86400
        cumulative_water += sum(r["volume_m3"] for r in packet["inflows"])
        cumulative_heat += sum(r["heat_j"] for r in packet["inflows"])
        close(live_rows[day, "accepted_water", 0], cumulative_water)
        close(live_rows[day, "accepted_heat", 0], cumulative_heat, 1e-5)
        residual = live_rows[day, "water_volume", 0] - initial_water - cumulative_water
        water_residual = max(water_residual, abs(residual))
        close(residual, 0, 1e-7)
        residual = live_rows[day, "heat_inventory", 0] - initial_heat - cumulative_heat
        heat_residual = max(heat_residual, abs(residual))
        close(residual, 0, max(.01, abs(initial_heat)*1e-12))
        for name in config["species"].values():
            close(live_rows[day, "accepted_"+name, 0], sum(sum(r["species_mass_kg"][name] for r in p["inflows"]) for p in packets if p["end_s"]<=packet["end_s"]), 1e-12)
            close(live_rows[day, "inventory_"+name, 0] - initial_water*.002, live_rows[day, "accepted_"+name, 0], 1e-9)
    assert cumulative_water>0, "Qualification must include a nonzero donor export"
    review = dict(status="passed", steps=len(packets), outfalls=len(config["outfalls"]),
        exchange_times_s=config["exchange_times_s"], mapped_water_m3=cumulative_water, mapped_heat_j=cumulative_heat,
        native_mapping_max_water_error_m3=max_mapping_error, maximum_water_inventory_residual_m3=water_residual,
        maximum_heat_inventory_residual_j=heat_residual, maximum_heat_inventory_relative_residual=heat_residual/(abs(initial_heat)+abs(cumulative_heat)), replay_compared_values=len(differences), maximum_live_replay_difference=max(differences),
        all_exports_acknowledged=True, manifest=str(folder/"live_output/run.json"), configuration=str(config_path), donor_artifacts=str(donor_folder))
    write(folder / "qualification.json", review)
    return review


def failure_gates(args):
    source=ROOT/"examples/swmm_reservoir/live"
    original=json.loads((source/"composition.json").read_text())
    configuration=json.loads((source/"configuration.json").read_text())
    run_folder=args.output/"failure_gates"/dt.datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    reviews={}
    cases={"epoch":"same physical epoch", "weights":"weights must sum", "mesh_order":"mesh geometry/order differs",
           "outfall":"Unknown or duplicate", "dual_forcing":"not both", "count_units":"count concentrations",
           "receiver_species":"invalid species identity", "reversal":"Outfall reversal", "multiple_receivers":"one direct acknowledged receiver"}
    for fault,reason in cases.items():
        folder=run_folder/fault
        config=copy.deepcopy(configuration)
        config["model"]=str((source/"network.inp").resolve())
        config["output_directory"]="donor_output"
        composition=copy.deepcopy(original)
        receiver=composition["components"][1]["arguments"]
        receiver["mesh"]={"@uri":str((source.parent/"reservoir.mesh.json").resolve())}
        if fault=="epoch": receiver["epoch_julian_day"]={"values":[config["epoch_julian_day"]+1]}
        if fault=="weights": config["outfalls"][0]["targets"][0]["weight"]=.1
        if fault=="mesh_order": config["receiver_geometry"]["x"][0]+=1
        if fault=="outfall": config["outfalls"][0]["id"]="nonexistent_outfall"
        if fault=="dual_forcing": receiver["recorded_inflows"]={"@uri":str((source/"replay_inflows.json").resolve())}
        if fault=="count_units":
            folder.mkdir(parents=True,exist_ok=True)
            model=folder/"counts.inp"
            model.write_text((source/"network.inp").read_text().replace("A MG/L", "A #/L"))
            config["model"]=str(model.resolve())
        if fault=="receiver_species": receiver["conservative_species"]={"tracer":.002}
        if fault=="reversal":
            folder.mkdir(parents=True,exist_ok=True)
            model=folder/"backflow.inp"
            model.write_text((source/"network.inp").read_text().replace("FLOW_ROUTING KINWAVE","FLOW_ROUTING DYNWAVE").replace("O 9 FREE NO","O 9 FIXED 100 NO"))
            config["model"]=str(model.resolve())
        if fault=="multiple_receivers":
            extra=copy.deepcopy(composition["components"][1]);extra["id"]="extra_receiver"
            composition["components"].append(extra)
            link=copy.deepcopy(composition["connections"][0]);link["to"]["component"]="extra_receiver"
            composition["connections"].append(link)
        write(folder/"configuration.json",config)
        write(folder/"composition.json",composition)
        run_host(args,folder/"composition.json",success=False)
        text=(folder/"composition.log").read_text()
        assert reason.lower() in text.lower(), text
        exports=list(folder.glob("donor_output/*/export_*.json"))
        acknowledgments=list(folder.glob("donor_output/*/accepted_*.json"))
        assert len(exports)==(1 if fault=="receiver_species" else 0), fault
        assert not acknowledgments, fault
        native_boundaries=list(folder.glob("donor_output/*/native_boundary_*.json"))
        assert len(native_boundaries)==(1 if fault in ("receiver_species","reversal") else 0), fault
        if fault=="reversal":
            audit=json.loads(native_boundaries[0].read_text())
            assert sum(row["cumulative_reversed_m3"] for row in audit["outfalls"])>0
        reviews[fault]={"status":"passed", "exports_before_failure":len(exports),
            "acknowledgments":len(acknowledgments), "log":str(folder/"composition.log")}
    write(run_folder/"qualification.json",reviews)
    return reviews


def factory_isolation(args):
    source=ROOT/"examples/swmm_reservoir/live"
    folder=args.output/"factory_isolation"/dt.datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    components=[]
    for name in ("first","second"):
        config=json.loads((source/"configuration.json").read_text())
        config.update(instance_id=name,model=str((source/"network.inp").resolve()),output_directory=name+"_output")
        write(folder/(name+".json"),config)
        components.append({"id":name,"info":{"component_info_id":"openswmm.model"},"arguments":{"configuration":{"@uri":name+".json"}}})
    write(folder/"composition.json",{"schema_version":"1.1","components":components,"workflow":{"strategy":"time_stepped"},
        "writers":[],"run":{"manifest":"run.json"}})
    run_host(args,folder/"composition.json")
    first=next((folder/"first_output").iterdir());second=next((folder/"second_output").iterdir())
    a=json.loads((first/"inflow_series.json").read_text());b=json.loads((second/"inflow_series.json").read_text())
    assert len(a["frames"])==len(b["frames"])==5
    for x,y in zip(a["frames"],b["frames"]):
        assert x["id"]!=y["id"] and x["inflows"]==y["inflows"]
    result={"status":"passed","instances":2,"shared_factory":True,"independent_transfer_ids":True,
        "matching_native_exports":True,"manifest":str(folder/"run.json")}
    write(folder/"qualification.json",result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--composer", type=Path, required=True)
    parser.add_argument("--fvqual", type=Path, default=ROOT/"build/apps/fvqual")
    parser.add_argument("--fvqual-components", type=Path, default=ROOT/"build/components")
    parser.add_argument("--swmm-components", type=Path, required=True)
    parser.add_argument("--engine-library", type=Path, required=True)
    parser.add_argument("--bellinge", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT/"tests/output/swmm_coupling_phase3")
    args = parser.parse_args()
    for key in ("composer", "fvqual", "fvqual_components", "swmm_components", "engine_library", "output"):
        setattr(args, key, getattr(args, key).resolve())
    reviews = {"synthetic_native_heat_species": qualify(args, ROOT/"examples/swmm_reservoir/live/composition.json")}
    reviews["failure_gates"]=failure_gates(args)
    reviews["factory_isolation"]=factory_isolation(args)
    if args.bellinge:
        reviews["bellinge"] = qualify(args, generate_bellinge(args))
    write(args.output / "qualification.json", reviews)
    print(json.dumps(reviews, indent=2))


if __name__ == "__main__":
    main()
