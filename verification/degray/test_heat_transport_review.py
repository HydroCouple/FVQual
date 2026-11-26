# SPDX-License-Identifier: MIT
import csv
import json
from pathlib import Path
import tempfile
import unittest

from heat_transport_review import CAPACITY, TERMS, evaluate


class SavedHeatAudit(unittest.TestCase):
    def fixture(self, directory, corrupt=None):
        (directory/"driver.json").write_text(json.dumps({"columns":1,"layers":1,"dam_column":0}))
        records=[]
        for day,energy in ((64.5,0),(65.5,100)):
            r={"jday":day,"column":0,"layer":0,"top_m":0,"bottom_m":-1,"volume_m3":100,
               "heat_J":10*100*CAPACITY+energy,"initial_heat_J":10*100*CAPACITY,
               "u_m_s":0,"v_m_s":0,"baroclinic_x_m_s2":0,"baroclinic_y_m_s2":0,
               "kz_bottom_m2_s":0}
            r.update({t:energy if t=="surface_J" else 0 for t in TERMS})
            r.update({f"upward_{process}_{side}_J":0 for process in ("advection","diffusion") for side in ("top","bottom")})
            records.append(r)
        if corrupt:corrupt(records[-1])
        with (directory/"degray_heat_transport.csv").open("w") as stream:
            writer=csv.DictWriter(stream,fieldnames=records[0]);writer.writeheader();writer.writerows(records)

    def test_temperature_attribution_closes_known_surface_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp);self.fixture(directory)
            series,_,review=evaluate(directory)
            self.assertTrue(review["audit_passed"])
            self.assertAlmostEqual(series[-1]["bottom_temperature_C"],10+100/(100*CAPACITY))
            self.assertAlmostEqual(series[-1]["surface_C"],100/(100*CAPACITY))

    def test_corrupt_heat_or_interface_flux_cannot_qualify(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            self.fixture(directory,lambda r:r.update(heat_J=r["heat_J"]+1e9))
            with self.assertRaisesRegex(ValueError,"attribution does not close"):evaluate(directory)
            self.fixture(directory,lambda r:r.update(upward_advection_top_J=1e9))
            _,_,review=evaluate(directory)
            self.assertFalse(review["audit_passed"])


if __name__=="__main__":unittest.main()
