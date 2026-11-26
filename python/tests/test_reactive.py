"""Native multispecies C ABI, owned history and independent forward oracles."""
import ctypes as ct
import math
import unittest
from fvqual import Library,Mesh,Simulation
from fvqual.engine import D,I64


class ReactiveTests(unittest.TestCase):
    def simulation(self,library,mesh,k=.025,initial=1.4):
        deck=f"[SPECIES]\nA kg/m3 {initial}\nB kg/m3 .2\n[CONSTANTS]\nk {k}\n[RATES]\nA WATER -k*A*(1+B)\nB WATER k*A*(1+B)\n"
        return Simulation(mesh,layers=2,options={"advectionScheme":0,"maxTimeStep":2},operations={"reaction_deck":deck})
    def test_native_sweep_adjoint_ownership_and_checkpoint(self):
        library=Library()
        with Mesh.channel(library,columns=2) as mesh:
            with self.simulation(library,mesh) as sim:
                sim.begin_reactive_trace();sim.advance(10);token=sim.save();sim.advance(10)
                trace=sim.capture_reactive_trace();self.assertEqual(trace.field(trace.states,"A")[-1],sim.snapshot("A")["values"])
                replay=trace.evaluate()
                for a,b in zip(replay[-1],trace.states[-1]):self.assertAlmostEqual(a,b,places=13)
                sim.restore(token);sim.advance(10);sim.release(token)
                with sim.capture_reactive_trace() as again:self.assertEqual(again.states,trace.states)
            with trace:
                tangent=trace.jvp(initial={"A":[.2]*trace.cells},constants={"k":.005})
                w=[[math.sin(i+.4*t) for i in range(trace.state_size)] for t in range(trace.steps+1)]
                bar=trace.vjp(w);left=sum(a*b for x,y in zip(w,tangent) for a,b in zip(x,y))
                right=.2*sum(bar["initial"]["A"])+.005*bar["constants"]["k"]
                self.assertAlmostEqual(left,right,places=11)
                for eps in (.01,.001,.0001):
                    runs=[]
                    for sign in (1,-1):
                        with self.simulation(library,mesh,.025+sign*eps*.005,1.4+sign*eps*.2) as fresh:
                            fresh.advance(20);runs.append(fresh.snapshot("A")["values"])
                    for a,b,d in zip(*runs,trace.field(tangent,"A")[-1]):self.assertLess(abs((a-b)/(2*eps)-d),eps*eps*.005+2e-8)
                with self.assertRaises(ValueError):trace.jvp(constants={"typo":1})
        with self.assertRaises(RuntimeError):trace.evaluate()
    def test_abi_errors_leave_every_output_untouched(self):
        library=Library()
        with Mesh.channel(library,columns=1) as mesh,self.simulation(library,mesh) as sim:
            sim.begin_reactive_trace();sim.advance(4)
            with sim.capture_reactive_trace() as trace:
                initial=(D*trace.state_size)(*([0]*trace.state_size));parameters=(D*1)(math.nan);states=(D*trace._count)(*([17]*trace._count))
                self.assertNotEqual(library.dll.fvq_reactive_trace_jvp(trace.handle,initial,len(initial),parameters,1,states,len(states)),0)
                self.assertEqual(list(states),[17]*len(states))
                w=(D*trace._count)(*([1]*trace._count));a=(D*trace.state_size)(*([23]*trace.state_size));p=(D*1)(29);peak=I64(31)
                self.assertNotEqual(library.dll.fvq_reactive_trace_vjp(trace.handle,w,len(w),a,len(a),p,0,ct.byref(peak)),0)
                self.assertEqual(list(a),[23]*len(a));self.assertEqual(p[0],29);self.assertEqual(peak.value,31)
                sim.end_reactive_trace();sim.set_reaction_constant("k",.03);sim.begin_reactive_trace();sim.advance(2)
                with sim.capture_reactive_trace() as changed:self.assertEqual(changed.parameters,{"k":.03})
    def test_unsupported_reactions_and_atomic_interval_rollback(self):
        library=Library()
        with Mesh.channel(library,columns=1) as mesh:
            with self.simulation(library,mesh) as sim:
                sim.begin_reactive_trace(max_steps=1);initial=sim.snapshot("A")
                frame={"id":"rollback","start_s":0,"end_s":6,"inflows":[{"column":0,"volume_m3":.01,"heat_j":418600.,"species_mass_kg":{"A":.01,"B":0.}}]}
                with self.assertRaisesRegex(RuntimeError,"rolled back"):sim.apply(frame)
                self.assertEqual(initial,sim.snapshot("A"))
                with sim.capture_reactive_trace() as trace:self.assertEqual(trace.steps,0)
            with Simulation(mesh,operations={"reaction_deck":"[RATES]\n__TEMPERATURE__ WATER 1\n"},options={"advectionScheme":0}) as bad:
                with self.assertRaises(RuntimeError):bad.begin_reactive_trace()


if __name__=="__main__":unittest.main()
