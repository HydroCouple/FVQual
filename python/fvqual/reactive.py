"""Native reactive transport derivatives conditional on prescribed hydraulics.

Rows contain all selected species blocks, each in face-major layer order.
Constants keep the reaction deck's units. No SWMM, thermal or geometry gradient
is implied. A trace owns immutable solver history and outlives the simulation.
"""
import ctypes as ct
import json
from .engine import D,I,I64,P


class Size(ct.Structure):
    _fields_=[(name,I64) for name in ("cells","species","steps","parameters")]


def bind(library):
    f=library._function
    f("fvq_sim_begin_reactive_trace",[P,I64]);f("fvq_sim_end_reactive_trace",[P])
    f("fvq_sim_capture_reactive_trace",[P,ct.POINTER(P)])
    f("fvq_reactive_trace_destroy",[P],None)
    f("fvq_reactive_trace_size",[P,ct.POINTER(Size)])
    f("fvq_reactive_trace_metadata",[P],ct.c_char_p)
    f("fvq_reactive_trace_baseline",[P,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(D),I64])
    for action in ("evaluate","jvp"):
        f("fvq_reactive_trace_"+action,[P,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(D),I64])
    f("fvq_reactive_trace_vjp",[P,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(I64)])
    f("fvq_sim_set_reaction_constant",[P,ct.c_char_p,D])


class ReactiveTrace:
    def __init__(self,library,handle):
        self.library,self.handle=library,handle
        try:
            bind(library);size=Size();library.check(library.dll.fvq_reactive_trace_size(handle,ct.byref(size)))
            self.cells,self.steps=size.cells,size.steps
            self.metadata=json.loads(library.dll.fvq_reactive_trace_metadata(handle))
            self.species=self.metadata["species"];self.parameter_names=self.metadata["parameter_names"]
            self.state_size=self.cells*len(self.species);self._count=self.state_size*(self.steps+1)
            x=(D*self._count)();p=(D*size.parameters)();t=(D*(self.steps+1))()
            library.check(library.dll.fvq_reactive_trace_baseline(handle,x,len(x),p,len(p),t,len(t)))
            self.states=self._reshape(x);self.initial=self.states[0][:]
            self.parameters=dict(zip(self.parameter_names,p));self.times=list(t)
        except BaseException:
            self.close();raise
    def _ready(self):
        if not self.handle:raise RuntimeError("Reactive trace is closed")
    def _reshape(self,x):return [list(x[i:i+self.state_size]) for i in range(0,len(x),self.state_size)]
    def field(self,states,species):
        if species not in self.species:raise ValueError("Unknown trace species")
        if len(states)!=self.steps+1 or any(len(row)!=self.state_size for row in states):raise ValueError("Reactive state dimensions differ")
        start=self.species.index(species)*self.cells
        return [row[start:start+self.cells] for row in states]
    def _initial(self,value,default):
        if value is None:return list(default)
        if isinstance(value,dict):
            if set(value)-set(self.species):raise ValueError("Unknown initial species")
            result=list(default)
            for name,values in value.items():
                if len(values)!=self.cells:raise ValueError("Initial species dimensions differ")
                start=self.species.index(name)*self.cells;result[start:start+self.cells]=values
            return result
        if len(value)!=self.state_size:raise ValueError("Initial state dimensions differ")
        return list(value)
    def _action(self,action,initial,constants):
        self._ready();default=self.initial if action=="evaluate" else [0.]*self.state_size
        a=self._initial(initial,default);p=dict(self.parameters) if action=="evaluate" else dict.fromkeys(self.parameter_names,0.)
        if constants:
            if set(constants)-set(p):raise ValueError("Unknown kinetic constant")
            p.update(constants)
        av=(D*len(a))(*a);pv=(D*len(p))(*(p[name] for name in self.parameter_names));out=(D*self._count)()
        self.library.check(getattr(self.library.dll,"fvq_reactive_trace_"+action)(self.handle,av,len(av),pv,len(pv),out,len(out)))
        return self._reshape(out)
    def evaluate(self,initial=None,constants=None):return self._action("evaluate",initial,constants)
    def jvp(self,initial=None,constants=None):return self._action("jvp",initial,constants)
    def vjp(self,cotangents):
        self._ready()
        if len(cotangents)!=self.steps+1 or any(len(row)!=self.state_size for row in cotangents):raise ValueError("Cotangents must have shape [steps+1,species*cells]")
        w=(D*self._count)(*(v for row in cotangents for v in row));x=(D*self.state_size)();p=(D*len(self.parameter_names))();peak=I64()
        self.library.check(self.library.dll.fvq_reactive_trace_vjp(self.handle,w,len(w),x,len(x),p,len(p),ct.byref(peak)))
        return {"initial":{name:list(x[i*self.cells:(i+1)*self.cells]) for i,name in enumerate(self.species)},
                "constants":dict(zip(self.parameter_names,p)),"peak_cell_tape_nodes":peak.value}
    def close(self):
        if self.handle:self.library.dll.fvq_reactive_trace_destroy(self.handle);self.handle=P()
    def __enter__(self):return self
    def __exit__(self,*_):self.close()
    def __del__(self):self.close()
