# SPDX-License-Identifier: MIT
"""Recorded scientific conversion preserves independent field values/axes."""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location("circulation_review",Path(__file__).parent.parent/"circulation_review.py")
review=importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


def example():
    geometry={"node_x":[0,0,1,1,2,2],"node_y":[0,1,0,1,0,1],"node_z":[-2]*6,
              "face_nodes":[0,2,3,1,2,4,5,3],"face_offsets":[0,4,8]}
    frames=[]
    for time,values in ((0,[10,11,12,13]),(60,[20,21,22,23])):
        frames.append({"time_s":time,"temperature":values,"tracer":[.1,.2,.3,.4],
            "velocity_x":[1,2,3,4],"velocity_y":[0,0,0,0],"pressure_x":[0]*4,
            "interfaces":[0,-1,-2,0,-1,-2],"viscosity":[0,1,0,0,2,0],
            "diffusivity":[0,3,0,0,4,0],"scalars":{"wind_work_J":12,"mode_impulse_x_kg_m_s":3}})
    return {"id":"fixture","columns":2,"layers":2,"geometry":geometry,"frames":frames}


class CirculationConversion(unittest.TestCase):
    def test_values_and_interface_coefficients_keep_column_layer_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=example();path=review.make_run(case,root,"fixture-sha")
            manifest=json.loads(path.read_text())
            with (path.parent/"fields.csv").open() as stream: rows=list(csv.DictReader(stream))
            temperatures=[float(r["value"]) for r in rows if r["item"]=="temperature"]
            self.assertEqual(temperatures,[10,11,12,13,20,21,22,23])
            self.assertEqual([float(r["value"]) for r in rows if r["item"]=="viscosity"],[0,1,0,2]*2)
            temp=next(r for r in manifest["results"] if r["item"]=="temperature")
            self.assertEqual(temp["shape"],[2,2,2]);self.assertEqual(temp["geometry"]["vertical"]["shape"],[2,2,3])
            self.assertEqual(temp["quantity"]["offset"],273.15)
            impulse=next(r for r in manifest["results"] if r["item"]=="mode_impulse_x_kg_m_s")
            self.assertEqual(impulse["quantity"]["powers"],[1,1,-1,0,0,0,0,0,0,0])
            self.assertEqual(len(set(float(r["time"]) for r in rows)),2)
            review.make_composition(root,{"baseline":path})
            sidecar=json.loads((root/"basin_review.composition.composer.json").read_text())
            self.assertEqual(sidecar["layers"][0]["runId"],manifest["run"]["id"])

    def test_corrupt_interface_and_time_records_are_refused(self):
        for mutation in (lambda c:c["frames"][1]["interfaces"].pop(),
                         lambda c:c["frames"][1].update(time_s=0),
                         lambda c:c["frames"][1]["temperature"].__setitem__(0,float("nan"))):
            case=copy.deepcopy(example());mutation(case)
            with tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError): review.make_run(case,Path(directory),"fixture-sha")


if __name__=="__main__": unittest.main()
