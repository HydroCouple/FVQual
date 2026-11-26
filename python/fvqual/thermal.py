"""Owned controlled thermal-column experiment using the native FVQual solver.

Analytic discrete derivatives include thermal/density/Richardson/shear feedback,
depth/area scaling and blocked adjoint replay. This is a restricted closed
column, not a full horizontal reservoir or SWMM adjoint.
"""
import ctypes as ct
import json
import math
from .engine import D,I,I64,P


class ColumnSize(ct.Structure):
    _fields_=[(name,I64) for name in ("stateSize","parameterCount","timeCount","checkpointCount")]+[
        ("forwardDifference",D),("minimumBuoyancyMargin",D),("derivativeSupported",I)]


class AdjointStats(ct.Structure):
    _fields_=[(name,I64) for name in ("recomputedSteps","peakTapeNodes","peakBlockStates")]


class ThermalColumn:
    def __init__(self,library,configuration=None):
        self.library,self.handle=library,P()
        f=library._function
        f("fvq_column_create",[ct.c_char_p,ct.POINTER(P)])
        f("fvq_column_destroy",[P],None);f("fvq_column_size",[P,ct.POINTER(ColumnSize)])
        f("fvq_column_metadata",[P],ct.c_char_p)
        f("fvq_column_baseline",[P,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(D),I64])
        for method in ("evaluate","jvp"):
            f("fvq_column_"+method,[P,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(D),I64])
        f("fvq_column_vjp",[P,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(AdjointStats)])
        library.check(library.dll.fvq_column_create(json.dumps(configuration or {},allow_nan=False).encode(),ct.byref(self.handle)))
        try:
            size=ColumnSize();library.check(library.dll.fvq_column_size(self.handle,ct.byref(size)))
            self.state_size,self.parameter_count,self.time_count=size.stateSize,size.parameterCount,size.timeCount
            self.layers=self.state_size//3;self.checkpoint_count=size.checkpointCount
            self.forward_difference=size.forwardDifference;self.minimum_buoyancy_margin=size.minimumBuoyancyMargin
            self.derivative_supported=bool(size.derivativeSupported)
            self.metadata=json.loads(library.dll.fvq_column_metadata(self.handle));self.parameter_names=self.metadata["parameters"]
            x=(D*self.state_size)();p=(D*self.parameter_count)();times=(D*self.time_count)()
            library.check(library.dll.fvq_column_baseline(self.handle,x,len(x),p,len(p),times,len(times)))
            self.initial=list(x);self.parameters=dict(zip(self.parameter_names,p));self.times=list(times)
            self.states=self.evaluate()
        except BaseException:
            self.close();raise
    def _ready(self):
        if not self.handle:raise RuntimeError("Thermal column is closed")
    def _parameter_array(self,parameters,baseline):
        if set(parameters)-set(self.parameter_names):raise ValueError("Unknown thermal parameter")
        values=dict(self.parameters) if baseline else dict.fromkeys(self.parameter_names,0.)
        values.update(parameters)
        return (D*self.parameter_count)(*(values[name] for name in self.parameter_names))
    def _action(self,method,initial,parameters):
        self._ready()
        if len(initial)!=self.state_size:raise ValueError("Thermal initial state dimensions differ")
        x=(D*self.state_size)(*initial);p=self._parameter_array(parameters,method=="evaluate")
        out=(D*(self.state_size*self.time_count))()
        self.library.check(getattr(self.library.dll,"fvq_column_"+method)(self.handle,x,len(x),p,len(p),out,len(out)))
        return [list(out[i:i+self.state_size]) for i in range(0,len(out),self.state_size)]
    def evaluate(self,initial=None,parameters=None):
        """Full native solver rerun, including any changed mixing regime."""
        return self._action("evaluate",self.initial if initial is None else initial,parameters or {})
    def jvp(self,initial=None,parameters=None):
        return self._action("jvp",[0.]*self.state_size if initial is None else initial,parameters or {})
    def vjp(self,state_cotangents):
        self._ready()
        if len(state_cotangents)!=self.time_count or any(len(row)!=self.state_size for row in state_cotangents):
            raise ValueError("Thermal cotangent dimensions differ")
        w=(D*(self.state_size*self.time_count))(*(x for row in state_cotangents for x in row))
        x=(D*self.state_size)();p=(D*self.parameter_count)();stats=AdjointStats()
        self.library.check(self.library.dll.fvq_column_vjp(self.handle,w,len(w),x,len(x),p,len(p),ct.byref(stats)))
        return {"initial":list(x),"parameters":dict(zip(self.parameter_names,p)),
                "stats":{name:getattr(stats,name) for name,_ in stats._fields_}}
    def close(self):
        if self.handle:self.library.dll.fvq_column_destroy(self.handle);self.handle=P()
    def __enter__(self):return self
    def __exit__(self,*_):self.close()
    def __del__(self):self.close()


class SigmaDepthObservation:
    """Fixed instrument depths sampled on a vertically scaled sigma column.

    Requests are (time_index, depth_below_surface_m). Outputs are temperatures.
    Interpolation weights and their geometry derivatives are included. No
    extrapolation; geometry derivatives refuse interpolation-stencil boundaries.
    """
    def __init__(self,shape,layers,base_depth,requests):
        self.shape=tuple(shape);self.layers=layers;self.depth=base_depth;self.requests=list(requests)
        if len(shape)!=2 or shape[0]<1 or shape[1]!=layers*3 or layers<1 or not math.isfinite(base_depth) or base_depth<=0:
            raise ValueError("Invalid sigma observation dimensions")
        for t,d in self.requests:
            if not isinstance(t,int) or not 0<=t<shape[0] or not math.isfinite(d) or d<0:raise ValueError("Invalid instrument time/depth")
        if not self.requests:raise ValueError("No instrument observations")
    def _rows(self,scale,derivative=False):
        if not math.isfinite(scale) or scale<=0:raise ValueError("Invalid geometry scale")
        result=[]
        for t,depth in self.requests:
            if self.layers==1:
                if depth>self.depth*scale:raise ValueError("Instrument outside the water column")
                result.append((t,0,0.,0.));continue
            q=depth*self.layers/(self.depth*scale)-.5
            if not 0<=q<self.layers-1:raise ValueError("Instrument depth outside centre support; extrapolation refused")
            k=math.floor(q);fraction=q-k
            if derivative and min(fraction,1-fraction)<1e-8:raise ValueError("Instrument on a changing interpolation stencil boundary")
            result.append((t,k,fraction,-depth*self.layers/(self.depth*scale*scale)))
        return result
    def _validate(self,states):
        if len(states)!=self.shape[0] or any(len(row)!=self.shape[1] for row in states):raise ValueError("Observation state shape differs")
        if any(not math.isfinite(x) for row in states for x in row):raise ValueError("Nonfinite observation state")
    def apply(self,states,depth_scale=1):
        self._validate(states)
        return [states[t][k] if self.layers==1 else (1-w)*states[t][k]+w*states[t][k+1]
                for t,k,w,_ in self._rows(depth_scale)]
    def jvp(self,states,tangents,depth_direction=0,depth_scale=1):
        self._validate(states);self._validate(tangents)
        if not math.isfinite(depth_direction):raise ValueError("Nonfinite geometry direction")
        return [tangents[t][k] if self.layers==1 else (1-w)*tangents[t][k]+w*tangents[t][k+1]+dw*(states[t][k+1]-states[t][k])*depth_direction
                for t,k,w,dw in self._rows(depth_scale,depth_direction!=0)]
    def vjp(self,states,cotangents,depth_scale=1):
        self._validate(states)
        if len(cotangents)!=len(self.requests) or any(not math.isfinite(x) for x in cotangents):raise ValueError("Invalid observation cotangents")
        out=[[0.]*self.shape[1] for _ in range(self.shape[0])];geometry=0.
        for (t,k,w,dw),bar in zip(self._rows(depth_scale,True),cotangents):
            out[t][k]+=(1-w)*bar
            if self.layers>1:out[t][k+1]+=w*bar;geometry+=dw*(states[t][k+1]-states[t][k])*bar
        return {"states":out,"depth_scale":geometry}
