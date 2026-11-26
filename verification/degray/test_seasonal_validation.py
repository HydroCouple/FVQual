# SPDX-License-Identifier: MIT
"""Observation-operator checks against small independent linear fields."""
import json
import csv
from pathlib import Path
import tempfile
import unittest
from seasonal_validation import sample, compare, read_casts, metrics
from review_artifacts import make_run

class ObservationOperator(unittest.TestCase):
    def test_split_recorded_fields_preserve_all_values_and_vertical_axes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            def table(name,records):
                with (root/name).open("w") as stream:
                    writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
            (root/"driver.json").write_text(json.dumps({"columns":2,"layers":1,"geometry_scope":"synthetic test"}))
            table("degray_mesh.csv",[{"column":c,"node_index":i,"x":x,"y":y}
                for c in (0,1) for i,(x,y) in enumerate(((c,0),(c+1,0),(c+1,1),(c,1)))])
            table("degray_field.csv",[{"jday":day,"column":c,"layer":0,"top_m":0,"bottom_m":-1,"temperature_C":value+c}
                for day,value in ((64.5,10),(65.5,20)) for c in (0,1)])
            table("degray_daily.csv",[dict(jday=day,surface_temperature_C=10,bottom_temperature_C=10,
                thermocline_depth_m=0,surface_elevation_m=0,volume_m3=2,heat_content_J=100,
                water_independent_residual_m3=0) for day in (64.5,65.5)])
            table("degray_budgets.csv",[dict(jday=day,budget_water_residual=0,budget_heat_residual=0) for day in (64.5,65.5)])
            table("degray_momentum.csv",[dict(jday=day,wind_work_J=day,mode_impulse_x_kg_m_s=-day) for day in (64.5,65.5)])
            make_run(root,"Shared field fixture")
            with (root/"gui_fields.csv").open() as stream: shared=list(csv.DictReader(stream))
            make_run(root,"Separate field fixture",split_fields=True)
            manifest=json.loads((root/"gui_run.json").read_text())
            for result in manifest["results"]+[manifest["results"][0]["geometry"]["vertical"]]:
                name=result["item"]
                with (root/result["artifact"]).open() as stream: actual=list(csv.DictReader(stream))
                self.assertEqual(actual,[row for row in shared if row["item"]==name])
            temperature=manifest["results"][0]
            self.assertEqual(temperature["artifact"],"gui_temperature.csv")
            self.assertEqual(temperature["geometry"]["vertical"]["artifact"],"gui_degray_interfaces.csv")
            momentum=next(r for r in manifest["results"] if r["item"]=="mode_impulse_x_kg_m_s")
            self.assertEqual(momentum["quantity"]["powers"],[1,1,-1,0,0,0,0,0,0,0])

    def test_exact_time_and_bracketing_preserve_linear_depth_time_field(self):
        profiles = {100.: [(1.,10.),(5.,18.)], 101.: [(1.,14.),(5.,22.)]}
        self.assertEqual(sample(profiles,100.,3.),14.)
        self.assertEqual(sample(profiles,100.25,3.),15.)
        self.assertIsNone(sample(profiles,99.9,3.))
        self.assertIsNone(sample(profiles,100.,0.))
        self.assertIsNone(sample({100.:profiles[100.],102.:profiles[101.]},101.,3.))

    def test_compare_separates_reserved_dates_and_explains_exclusions(self):
        profiles = {77.5: [(1.,11.),(3.,15.)], 273.5: [(1.,12.),(3.,16.)]}
        casts = [{"jday":63.5,"depth_m":[0,2],"temperature_c":[10,10]},
                 {"jday":77.5,"depth_m":[0,1,2],"temperature_c":[10,10,12]},
                 {"jday":273.5,"depth_m":[1,2],"temperature_c":[10,14]}]
        summary, points = compare(profiles,casts)
        self.assertEqual(summary["development"]["n"],2)
        self.assertEqual(summary["reserved_review"]["n"],2)
        self.assertEqual(summary["development"]["bias_c"],1.)
        self.assertEqual(summary["reserved_review"]["bias_c"],1.)
        self.assertEqual(len(summary["excluded"]),2)
        self.assertEqual(metrics([3.,4.])["rmse_c"],(12.5)**.5)

    def test_strict_observed_parser_rejects_duplicate_depth_and_truncation(self):
        header="CNO DAYO NRO\n22 77.5 2\nZO C2O\n"
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/"casts.npt"
            path.write_text(header+"0 10 2 14\n")
            self.assertEqual(read_casts(path)[0]["temperature_c"],[10,14])
            for data in ("0 10 0 14\n","0 10\nCNO DAYO NRO\n22 80 1\nZO C2O\n0 11\n"):
                path.write_text(header+data)
                with self.assertRaises(ValueError): read_casts(path)

if __name__=="__main__": unittest.main()
