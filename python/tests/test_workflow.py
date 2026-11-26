"""C ABI integration gates; artifacts remain in tests/output for review."""
import copy
import json
import math
import os
from pathlib import Path
import unittest
import uuid

from fvqual.engine import Library, Mesh, Simulation
from fvqual.recorded import forward, read_case
from fvqual.uncertainty import ensemble, sensitivity
from fvqual.coupled import run as coupled_run

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples/swmm_reservoir"


class WorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = Library()
        cls.path, cls.spec = read_case(EXAMPLE / "recorded_event.json")
        cls.output = ROOT / "tests/output" / ("python_" + uuid.uuid4().hex)
        cls.output.mkdir(parents=True)

    def test_owned_simulation_outlives_mesh_and_replays_receipts(self):
        mesh = Mesh.channel(self.library, columns=2)
        simulation = Simulation(mesh, species={"tracer": {"units": "kg/m3"}})
        mesh.close()
        mesh.close()
        with simulation:
            initial = simulation.volume
            token = simulation.save()
            frame = {"id": "one", "start_s": 0, "end_s": .1,
                     "inflows": [{"column": 0, "volume_m3": .001,
                                  "heat_j": .001 * 1000 * 4186 * 12,
                                  "species_mass_kg": {"tracer": .0002}}]}
            receipt = simulation.apply(frame)
            final = (simulation.volume, simulation.mean_temperature, simulation.mass("tracer"))
            self.assertFalse(receipt["replayed"])
            self.assertTrue(simulation.apply(frame)["replayed"])
            self.assertEqual(final, (simulation.volume, simulation.mean_temperature, simulation.mass("tracer")))
            changed = copy.deepcopy(frame)
            changed["inflows"][0]["volume_m3"] *= 2
            with self.assertRaises(RuntimeError):
                simulation.apply(changed)
            self.assertEqual(simulation.time, .1)
            simulation.restore(token)
            self.assertEqual(simulation.volume, initial)
            self.assertEqual(simulation.time, 0)
            self.assertFalse(simulation.apply(frame)["replayed"])
            self.assertEqual(final, (simulation.volume, simulation.mean_temperature, simulation.mass("tracer")))
            simulation.release(token)
            with self.assertRaises(RuntimeError):
                simulation.restore(token)
        simulation.close()
        with self.assertRaises(RuntimeError):
            _ = simulation.time

    def test_recorded_mixed_mesh_budgets_and_parameter_reference(self):
        result = forward(self.library, self.spec, self.path.parent)
        (self.output / "recorded.json").write_text(json.dumps(result, indent=2))
        self.assertAlmostEqual(result["final"]["volume_m3"], 1000.07, places=10)
        self.assertAlmostEqual(result["final"]["species_mass_kg"]["tracer"], .014, places=12)
        self.assertLess(abs(result["budget_residuals"]["water_m3"]), 1e-9)
        self.assertLess(abs(result["budget_residuals"]["heat_j"]), 1e-3)
        self.assertLess(abs(result["budget_residuals"]["species_kg"]["tracer"]), 1e-12)
        for row in sensitivity(self.library, self.spec, self.path.parent, "load_scale", 1):
            self.assertAlmostEqual(row["derivative"]["mass_kg:tracer"], .014, places=10)
            self.assertEqual(row["derivative"]["volume_m3"], 0)

    def test_kernel_pullbacks_match_directional_differences(self):
        for kind, values in ((0, [2, 1e-6, 3600]),
                             (1, [10, 200, 86400, 10, 1000, 4186]),
                             (2, [1000, 10, .1, 1.2])):
            direction = [max(abs(v), .01) * .001 for v in values]
            tangent = self.library.kernel(kind, values, direction=direction)
            pullback = self.library.kernel(kind, values, cotangent=.73)
            self.assertAlmostEqual(.73 * tangent, sum(a*b for a,b in zip(pullback, direction)), places=12)
            for epsilon in (.01, .001, .0001):
                plus = self.library.kernel(kind, [a + epsilon*b for a,b in zip(values, direction)])
                minus = self.library.kernel(kind, [a - epsilon*b for a,b in zip(values, direction)])
                self.assertLess(abs((plus-minus)/(2*epsilon)-tangent), 1e-8)
        with self.assertRaises(RuntimeError):
            self.library.kernel(0, [1, -1, 1], cotangent=1)

    def test_ensemble_seed_reproduces_bands_and_counts_every_run(self):
        a = ensemble(self.library, self.path, self.spec, self.output / "ensemble_a", samples=8, seed=42)
        b = ensemble(self.library, self.path, self.spec, self.output / "ensemble_b", samples=8, seed=42)
        self.assertEqual(a["prediction_bands"], b["prediction_bands"])
        self.assertEqual(a["successful"], 8)
        self.assertEqual(a["failed"], 0)
        self.assertEqual(len(list((self.output / "ensemble_a").glob("sample_*/result.json"))), 8)
        self.assertEqual(a["successful"] + a["failed"], a["samples"])

    @unittest.skipUnless(os.environ.get("SWMM_ENGINE_LIBRARY"), "live SWMM library is optional")
    def test_live_swmm_closes_combined_budgets_and_matches_recording(self):
        result = coupled_run(self.library, os.environ["SWMM_ENGINE_LIBRARY"], self.path,
                             EXAMPLE / "surface_donor.inp", self.output / "live")
        self.assertLess(abs(result["combined_budget_residuals"]["water_m3"]), 1e-9)
        self.assertLess(abs(result["combined_budget_residuals"]["species_kg"]), 1e-12)
        self.assertLess(abs(result["combined_budget_residuals"]["heat_j"]), 1e-3)
        self.assertTrue(all(value == 0 for value in result["live_recorded_difference"].values()))


if __name__ == "__main__":
    unittest.main()
