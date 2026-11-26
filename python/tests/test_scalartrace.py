"""Actual-solver C ABI, observation transpose and conditional covariance gates."""
import ctypes as ct
import math
import unittest
from fvqual import Library,Mesh,Simulation
from fvqual.engine import D
from fvqual.scalar_uncertainty import FixedObservation,local_covariance


class ScalarTraceTests(unittest.TestCase):
    def simulation(self,library,mesh):
        return Simulation(mesh,layers=3,species={"tracer":{"units":"kg/m3","initial":.002}},
                          options={"advectionScheme":0,"maxTimeStep":5,"backgroundDiffusivity":.001})
    def frame(self,start=0,end=10):
        return {"id":str(start),"start_s":start,"end_s":end,"inflows":[
            {"column":0,"volume_m3":.1,"heat_j":.1*1000*4186*10,"placement":"surface",
             "species_mass_kg":{"tracer":.0003}}]}
    def test_actual_coupled_derivative_snapshot_replay_and_lifetime(self):
        library=Library()
        with Mesh.channel(library,columns=2) as mesh:
            with self.simulation(library,mesh) as sim:
                sim.begin_scalar_trace("tracer");token=sim.save()
                sim.apply(self.frame());snap=sim.snapshot("tracer")
                trace=sim.capture_scalar_trace();self.assertEqual(trace.states[-1],snap["values"])
                before=(sim.time,sim.mass("tracer"),sim.mass("temperature"))
                tangent=trace.jvp(sources=trace.sources)
                self.assertAlmostEqual(sum(a*b for a,b in zip(tangent[-1],snap["volumes"])),.0003,places=13)
                for a,b in zip(trace.evaluate()[-1],snap["values"]): self.assertAlmostEqual(a,b,places=15)
                sim.restore(token);sim.apply(self.frame());sim.release(token)
                with sim.capture_scalar_trace() as replay: self.assertEqual(replay.states,trace.states)
                self.assertEqual(before,(sim.time,sim.mass("tracer"),sim.mass("temperature")))
        with trace:
            self.assertTrue(trace.evaluate());self.assertGreater(trace.storage_bytes,0)
        with self.assertRaises(RuntimeError):trace.jvp()
    def test_mapping_dot_product_covariance_and_atomic_abi_errors(self):
        library=Library()
        with Mesh.channel(library,columns=2) as mesh, self.simulation(library,mesh) as sim:
            sim.begin_scalar_trace("tracer");sim.apply(self.frame())
            with sim.capture_scalar_trace() as trace:
                last=trace.steps
                observation=FixedObservation((last+1,trace.cells),[
                    [(last,0,.25),(last,1,.75)],[(last,0,.5),(last-1,0,.5)],
                    [(last,0,1),(last,0,1)]],scale=1000,offset=7)
                initial=[.0001]*trace.cells;sources=[.0002]*trace.source_count
                tangent=trace.jvp(initial,sources,.2);w=[.3,-.7,.2]
                bar=trace.vjp(observation.vjp(w))
                left=sum(a*b for a,b in zip(w,observation.jvp(tangent)))
                right=sum(a*b for a,b in zip(initial,bar["initial"]))+sum(a*b for a,b in zip(sources,bar["sources"]))+.2*bar["diffusion_multiplier"]
                self.assertAlmostEqual(left,right,places=13)
                result=local_covariance(trace,[{"initial":initial,"sources":sources}],observation)
                column=observation.jvp(trace.jvp(initial,sources))
                for i in range(3): self.assertAlmostEqual(result["covariance"][i][i],column[i]**2,places=14)
                # C ABI must validate ALL arrays before writing any output.
                count=trace.cells*(trace.steps+1);out=(D*count)(*([42]*count))
                a=(D*trace.cells)(*([0]*trace.cells));b=(D*trace.source_count)(*([0]*trace.source_count))
                code=library.dll.fvq_scalar_trace_jvp(trace.handle,a,len(a),b,len(b),math.nan,out,count)
                self.assertNotEqual(code,0);self.assertEqual(list(out),[42]*count)
                with self.assertRaises(RuntimeError):trace.evaluate(diffusion_multiplier=0)
                with self.assertRaises(ValueError):observation.vjp([1])
    def test_unsupported_modes_and_interval_rollback(self):
        library=Library()
        with Mesh.channel(library,columns=2) as mesh:
            with Simulation(mesh,species={"tracer":{"units":"kg/m3"}}) as sim:
                with self.assertRaises(RuntimeError):sim.begin_scalar_trace("tracer")
                with self.assertRaises(RuntimeError):sim.begin_scalar_trace("temperature")
            with self.simulation(library,mesh) as sim:
                sim.begin_scalar_trace("tracer",max_steps=1)
                initial=sim.snapshot("tracer")
                with self.assertRaisesRegex(RuntimeError,"rolled back"):sim.apply(self.frame())
                self.assertEqual(sim.snapshot("tracer"),initial)
                with sim.capture_scalar_trace() as trace: self.assertEqual(trace.steps,0)
                sim.end_scalar_trace();sim.apply(self.frame())


if __name__=="__main__":unittest.main()
