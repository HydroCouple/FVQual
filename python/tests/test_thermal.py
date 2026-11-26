"""Native C ABI thermal, moving-depth observation and replay qualification."""
import ctypes as ct
import math
import unittest
from fvqual import Library,ThermalColumn,SigmaDepthObservation
from fvqual.engine import D


class ThermalTests(unittest.TestCase):
    def test_geometry_observation_and_true_transpose_match_native_differences(self):
        with ThermalColumn(Library(),{"layers":6,"duration_s":300,"dt_s":10}) as model:
            self.assertTrue(model.derivative_supported)
            observation=SigmaDepthObservation((model.time_count,model.state_size),6,10,[(30,2.6),(30,7.4),(15,2.6),(30,2.6)])
            dx=[.1]*6+[.002*math.sin(i) for i in range(12)]
            dp={"wind_speed_m_s":.3,"depth_scale":.025,"shortwave_w_m2":10}
            tangent=model.jvp(dx,dp);predicted=observation.jvp(model.states,tangent,.025)
            for epsilon in (.01,.001,.0001,.00001):
                samples=[]
                for sign in (1,-1):
                    initial=[x+sign*epsilon*v for x,v in zip(model.initial,dx)]
                    parameters={key:model.parameters[key]+sign*epsilon*value for key,value in dp.items()}
                    actual=model.evaluate(initial,parameters)
                    samples.append(observation.apply(actual,1+sign*epsilon*.025))
                for a,b,d in zip(*samples,predicted):self.assertLess(abs((a-b)/(2*epsilon)-d),3e-7)
            weights=[.3,-.7,.2,.4];mapped=observation.vjp(model.states,weights);bar=model.vjp(mapped["states"])
            left=sum(a*b for a,b in zip(weights,predicted))
            right=sum(a*b for a,b in zip(dx,bar["initial"]))+sum(v*bar["parameters"][k] for k,v in dp.items())+.025*mapped["depth_scale"]
            self.assertAlmostEqual(left,right,places=12)
            # Omitting the geometry part is a deliberate derivative mutation.
            self.assertGreater(abs(left-(right-.025*mapped["depth_scale"])),.001)
            self.assertEqual(bar["stats"]["recomputedSteps"],30)
    def test_branch_and_abi_failure_leave_inputs_and_outputs_intact(self):
        library=Library()
        with ThermalColumn(library,{"layers":4,"initial":[20]*4+[0]*8,"duration_s":10}) as model:
            self.assertFalse(model.derivative_supported)
            state=model.states
            with self.assertRaisesRegex(RuntimeError,"branch"):model.jvp(parameters={"depth_scale":.01})
            self.assertEqual(model.evaluate(),state)
            n=model.state_size*model.time_count
            x=(D*model.state_size)(*([0]*model.state_size));p=(D*model.parameter_count)(*([0]*model.parameter_count));out=(D*n)(*([42]*n))
            self.assertNotEqual(library.dll.fvq_column_jvp(model.handle,x,len(x),p,len(p),out,n),0)
            self.assertEqual(list(out),[42]*n)
        with self.assertRaises(RuntimeError):model.evaluate()
        for config in ({"layers":2.5},{"manning":.03},{"weather_constants":{"unsupported":1}}):
            with self.assertRaises(RuntimeError):ThermalColumn(library,config)
    def test_observation_refuses_extrapolation_and_stencil_changes(self):
        o=SigmaDepthObservation((2,12),4,10,[(1,3.75)])
        values=[[20,18,16,14]+[0]*8]*2
        self.assertEqual(o.apply(values),[18])
        with self.assertRaisesRegex(ValueError,"stencil"):o.jvp(values,values,.01)
        with self.assertRaises(ValueError):SigmaDepthObservation((2,12),4,10,[(1,.1)]).apply(values)
        with self.assertRaises(ValueError):o.apply(values,0)


if __name__=="__main__":unittest.main()
