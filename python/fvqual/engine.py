"""Owned handles and conservative interval transfers. SI units throughout."""
import ctypes as ct
import math
import json
import os
from pathlib import Path

D = ct.c_double
I = ct.c_int
I64 = ct.c_int64
P = ct.c_void_p

class Options(ct.Structure):
    _fields_ = [(name, D) for name in (
        "gravity referenceDensity coriolis neutralViscosity backgroundViscosity "
        "backgroundDiffusivity convectiveMixing horizontalViscosity horizontalDiffusivity"
    ).split()] + [("advectionScheme", I)] + [(name, D) for name in (
        "windStressX windStressY surfaceHeatFlux heatCapacity manningN cflNumber minTimeStep maxTimeStep"
    ).split()]

class Sigma(ct.Structure):
    _fields_ = [("distribution",I),("layerCount",I),("theta",D),("blend",D),
                ("minLayerThickness",D),("interfaces",ct.POINTER(D)),("interfaceCount",I)]

class InflowSpecies(ct.Structure):
    _fields_ = [("id",ct.c_char_p),("massKg",D)]

class Inflow(ct.Structure):
    _fields_ = [("column",I64),("volumeM3",D),("heatJ",D),("placement",I),
                ("elevationM",D),("species",ct.POINTER(InflowSpecies)),("speciesCount",I64)]

class Receipt(ct.Structure):
    _fields_ = [(name,D) for name in (
        "requestedVolumeM3 acceptedVolumeM3 requestedHeatJ acceptedHeatJ").split()] + [("replayed",I)]

class TraceSize(ct.Structure):
    _fields_=[(name,I64) for name in ("cells","steps","sources","storageBytes")]

class Library:
    def __init__(self, path=None):
        if path is None:
            path = os.environ.get("FVQUAL_LIBRARY")
        if path is None:
            root = Path(__file__).resolve().parents[2]
            candidates = [root/"build/src/libfvqual-core.dylib", root/"build/src/libfvqual-core.so",
                          root/"build/bin/fvqual-core.dll"]
            path = next((p for p in candidates if p.exists()), None)
        if path is None:
            raise FileNotFoundError("Set FVQUAL_LIBRARY to the shared FVQual core library")
        self.path = Path(path).resolve()
        self.dll = ct.CDLL(str(self.path))
        self._bind()
    def _function(self, name, args, result=I):
        f = getattr(self.dll, name)
        f.argtypes, f.restype = args, result
        return f
    def _bind(self):
        f=self._function
        f("fvq_last_error",[],ct.c_char_p)
        f("fvq_version",[],ct.c_char_p)
        f("fvq_mesh_create_channel",[I,D,D,D,ct.POINTER(P)])
        f("fvq_mesh_read_package",[ct.c_char_p,ct.POINTER(P)])
        f("fvq_mesh_read_ugrid",[ct.c_char_p,ct.c_char_p,ct.POINTER(P)])
        f("fvq_mesh_create",[ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),I64,
                             ct.POINTER(I64),ct.POINTER(I64),I64,ct.POINTER(P)])
        f("fvq_mesh_destroy",[P],None)
        f("fvq_mesh_total_area",[P],D)
        f("fvq_options_default",[ct.POINTER(Options)],None)
        f("fvq_sigma_default",[ct.POINTER(Sigma)],None)
        f("fvq_sim_create",[P,ct.POINTER(Sigma),ct.POINTER(Options),ct.POINTER(P)])
        f("fvq_sim_destroy",[P],None)
        f("fvq_sim_add_temperature",[P,D,ct.POINTER(I)])
        f("fvq_sim_add_species",[P,ct.c_char_p,ct.c_char_p,D,ct.POINTER(I)])
        f("fvq_sim_initialize",[P])
        f("fvq_sim_advance",[P,D])
        f("fvq_sim_time",[P],D)
        f("fvq_sim_column_count",[P],I64)
        f("fvq_sim_layer_count",[P],I)
        f("fvq_sim_cell_count",[P],I64)
        f("fvq_sim_scalar_snapshot",[P,I,ct.POINTER(D),ct.POINTER(D),I64,ct.POINTER(D),I64])
        f("fvq_sim_total_volume",[P],D)
        f("fvq_sim_species_mass",[P,I,ct.POINTER(D)])
        f("fvq_sim_apply_inflow_frame",[P,ct.c_char_p,D,D,ct.POINTER(Inflow),I64,ct.POINTER(Receipt)])
        f("fvq_sim_transfer_species_receipt",[P,ct.c_char_p,ct.c_char_p,ct.POINTER(D)])
        f("fvq_sim_save_state",[P,ct.POINTER(ct.c_uint64)])
        f("fvq_sim_restore_state",[P,ct.c_uint64])
        f("fvq_sim_release_state",[P,ct.c_uint64])
        f("fvq_kernel_value",[I,ct.POINTER(D),I,ct.POINTER(D)])
        f("fvq_kernel_jvp",[I,ct.POINTER(D),ct.POINTER(D),I,ct.POINTER(D)])
        f("fvq_kernel_vjp",[I,ct.POINTER(D),I,D,ct.POINTER(D)])
        f("fvq_sim_begin_scalar_trace",[P,I,I64])
        f("fvq_sim_end_scalar_trace",[P])
        f("fvq_sim_capture_scalar_trace",[P,ct.POINTER(P)])
        f("fvq_scalar_trace_destroy",[P],None)
        f("fvq_scalar_trace_size",[P,ct.POINTER(TraceSize)])
        f("fvq_scalar_trace_baseline",[P,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(D),I64])
        f("fvq_scalar_trace_source_id",[P,I64],ct.c_char_p)
        for name in ("evaluate","jvp"):
            f("fvq_scalar_trace_"+name,[P,ct.POINTER(D),I64,ct.POINTER(D),I64,D,ct.POINTER(D),I64])
        f("fvq_scalar_trace_vjp",[P,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(D),I64,ct.POINTER(D)])
    def check(self, code):
        if code:
            raise RuntimeError(self.dll.fvq_last_error().decode("utf8",errors="replace"))
    @property
    def version(self):
        return self.dll.fvq_version().decode()
    def kernel(self, kind, x, direction=None, cotangent=None):
        """D1 smooth kernels: decay=0, fully mixed column heat=1, conservative mix=2.
        No full FVQual step, adaptive schedule or SWMM derivative is implied.
        """
        values=(D*len(x))(*x)
        if direction is not None and cotangent is not None:
            raise ValueError("Choose a JVP direction or a VJP cotangent")
        if cotangent is not None:
            result=(D*len(x))()
            self.check(self.dll.fvq_kernel_vjp(kind,values,len(x),cotangent,result))
            return list(result)
        result=D()
        if direction is None:
            self.check(self.dll.fvq_kernel_value(kind,values,len(x),ct.byref(result)))
        else:
            if len(direction)!=len(x): raise ValueError("Tangent dimension mismatch")
            v=(D*len(x))(*direction)
            self.check(self.dll.fvq_kernel_jvp(kind,values,v,len(x),ct.byref(result)))
        return result.value

class Mesh:
    def __init__(self, library, handle):
        self.library,self.handle=library,handle
    @classmethod
    def channel(cls, library, columns=4, length=100, width=20, bed=-5):
        h=P()
        library.check(library.dll.fvq_mesh_create_channel(columns,length,width,bed,ct.byref(h)))
        return cls(library,h)
    @classmethod
    def read(cls, library, path):
        path=Path(path).resolve()
        h=P()
        if path.suffix==".json":
            code=library.dll.fvq_mesh_read_package(os.fsencode(path),ct.byref(h))
        else:
            code=library.dll.fvq_mesh_read_ugrid(os.fsencode(path),None,ct.byref(h))
        library.check(code)
        return cls(library,h)
    @classmethod
    def from_arrays(cls, library, x, y, z, offsets, nodes):
        if len(x)!=len(y) or len(x)!=len(z): raise ValueError("Node array lengths differ")
        if len(offsets)<2 or offsets[0]!=0 or offsets[-1]!=len(nodes):
            raise ValueError("Invalid CSR offsets")
        if any(a>b for a,b in zip(offsets,offsets[1:])):
            raise ValueError("CSR offsets are not monotone")
        arrays=[(D*len(a))(*a) for a in (x,y,z)]
        off=(I64*len(offsets))(*offsets); nd=(I64*len(nodes))(*nodes); h=P()
        library.check(library.dll.fvq_mesh_create(*arrays,len(x),off,nd,len(offsets)-1,ct.byref(h)))
        return cls(library,h)
    @property
    def area(self):
        if not self.handle: raise RuntimeError("Mesh is closed")
        return self.library.dll.fvq_mesh_total_area(self.handle)
    def close(self):
        if self.handle:
            self.library.dll.fvq_mesh_destroy(self.handle); self.handle=P()
    def __enter__(self): return self
    def __exit__(self,*_): self.close()
    def __del__(self): self.close()

class Simulation:
    """One owned simulation; methods on the same instance must be serialized."""
    placements={"uniform":0,"surface":1,"bed":2,"elevation":3,"density-matched":4}
    def __init__(self, mesh, layers=4, temperature=10, species=None, options=None, operations=None):
        self.library=mesh.library
        self.handle=P()
        self.species={}
        if not mesh.handle: raise RuntimeError("Mesh is closed")
        sigma=Sigma(); opts=Options()
        self.library.dll.fvq_sigma_default(ct.byref(sigma))
        self.library.dll.fvq_options_default(ct.byref(opts))
        sigma.layerCount=layers
        sigma.minLayerThickness=0
        allowed={name for name,_ in Options._fields_}
        for name,value in (options or {}).items():
            if name not in allowed: raise ValueError("Unknown solver option: "+name)
            if not math.isfinite(value): raise ValueError("Nonfinite option: "+name)
            setattr(opts,name,value)
        self.options=opts
        self.area=mesh.area
        self.library.check(self.library.dll.fvq_sim_create(mesh.handle,ct.byref(sigma),ct.byref(opts),ct.byref(self.handle)))
        try:
            index=I()
            self.library.check(self.library.dll.fvq_sim_add_temperature(self.handle,temperature,ct.byref(index)))
            self.temperature_index=index.value
            for name,definition in (species or {}).items():
                self.library.check(self.library.dll.fvq_sim_add_species(self.handle,name.encode(),
                    definition["units"].encode(),definition.get("initial",0),ct.byref(index)))
                self.species[name]=index.value
            if operations is not None:
                self.library._function("fvq_sim_configure_json",[P,ct.c_char_p])
                self.library._function("fvq_sim_species_count",[P],I)
                self.library._function("fvq_sim_species_id",[P,I],ct.c_char_p)
                self.library.check(self.library.dll.fvq_sim_configure_json(self.handle,json.dumps(operations,allow_nan=False).encode()))
                for i in range(self.library.dll.fvq_sim_species_count(self.handle)):
                    if i!=self.temperature_index:self.species[self.library.dll.fvq_sim_species_id(self.handle,i).decode()]=i
            self.library.check(self.library.dll.fvq_sim_initialize(self.handle))
        except BaseException:
            self.close()
            raise
    def _ready(self):
        if not self.handle: raise RuntimeError("Simulation is closed")
    @property
    def time(self):
        self._ready(); return self.library.dll.fvq_sim_time(self.handle)
    @property
    def volume(self):
        self._ready(); return self.library.dll.fvq_sim_total_volume(self.handle)
    def mass(self, species):
        self._ready(); value=D()
        index=self.temperature_index if species=="temperature" else self.species[species]
        self.library.check(self.library.dll.fvq_sim_species_mass(self.handle,index,ct.byref(value)))
        return value.value
    @property
    def mean_temperature(self): return self.mass("temperature")/self.volume
    def advance(self,duration):
        self._ready(); self.library.check(self.library.dll.fvq_sim_advance(self.handle,duration))
    def budget(self,name):
        self._ready();self.library._function("fvq_sim_budget_value",[P,ct.c_char_p,ct.POINTER(D)])
        value=D();self.library.check(self.library.dll.fvq_sim_budget_value(self.handle,name.encode(),ct.byref(value)));return value.value
    def set_reaction_constant(self,name,value):
        self._ready();from .reactive import bind
        bind(self.library);self.library.check(self.library.dll.fvq_sim_set_reaction_constant(self.handle,name.encode(),value))
    def begin_reactive_trace(self,max_steps=10000):
        self._ready();from .reactive import bind
        bind(self.library);self.library.check(self.library.dll.fvq_sim_begin_reactive_trace(self.handle,max_steps))
    def end_reactive_trace(self):
        self._ready();from .reactive import bind
        bind(self.library);self.library.check(self.library.dll.fvq_sim_end_reactive_trace(self.handle))
    def capture_reactive_trace(self):
        self._ready();from .reactive import bind,ReactiveTrace
        bind(self.library);handle=P();self.library.check(self.library.dll.fvq_sim_capture_reactive_trace(self.handle,ct.byref(handle)));return ReactiveTrace(self.library,handle)
    def apply(self,frame):
        self._ready()
        records=[]
        keep=[]
        for data in frame["inflows"]:
            species=(InflowSpecies*len(data.get("species_mass_kg",{})))(
                *(InflowSpecies(k.encode(),v) for k,v in sorted(data.get("species_mass_kg",{}).items())))
            keep.append(species)
            records.append(Inflow(data["column"],data["volume_m3"],data["heat_j"],
                self.placements[data.get("placement","uniform")],data.get("elevation_m",0),species,len(species)))
        transfers=(Inflow*len(records))(*records); result=Receipt()
        self.library.check(self.library.dll.fvq_sim_apply_inflow_frame(self.handle,frame["id"].encode(),
            frame["start_s"],frame["end_s"],transfers,len(records),ct.byref(result)))
        output={name:getattr(result,name) for name,_ in Receipt._fields_}
        output["id"]=frame["id"]
        output["species_mass_kg"]={}
        for name in self.species:
            value=D()
            # Empty frames have no species receipts.
            if any(name in d.get("species_mass_kg",{}) for d in frame["inflows"]):
                self.library.check(self.library.dll.fvq_sim_transfer_species_receipt(self.handle,
                    frame["id"].encode(),name.encode(),ct.byref(value)))
            output["species_mass_kg"][name]=value.value
        return output
    def save(self):
        self._ready(); token=ct.c_uint64()
        self.library.check(self.library.dll.fvq_sim_save_state(self.handle,ct.byref(token)))
        return token.value
    def begin_scalar_trace(self, species, max_steps=10000):
        """Actual passive upwind transport; unsupported physics is refused."""
        self._ready()
        index=self.temperature_index if species=="temperature" else self.species[species]
        self.library.check(self.library.dll.fvq_sim_begin_scalar_trace(self.handle,index,max_steps))
    def snapshot(self,species):
        self._ready()
        n=self.library.dll.fvq_sim_cell_count(self.handle)
        c=self.library.dll.fvq_sim_column_count(self.handle);k=self.library.dll.fvq_sim_layer_count(self.handle)
        values=(D*n)();volumes=(D*n)();interfaces=(D*(c*(k+1)))()
        index=self.temperature_index if species=="temperature" else self.species[species]
        self.library.check(self.library.dll.fvq_sim_scalar_snapshot(self.handle,index,values,volumes,n,interfaces,len(interfaces)))
        return {"time_s":self.time,"columns":c,"layers":k,"values":list(values),"volumes":list(volumes),"interfaces":list(interfaces)}
    def end_scalar_trace(self):
        self._ready();self.library.check(self.library.dll.fvq_sim_end_scalar_trace(self.handle))
    def capture_scalar_trace(self):
        self._ready();handle=P()
        self.library.check(self.library.dll.fvq_sim_capture_scalar_trace(self.handle,ct.byref(handle)))
        return ScalarTrace(self.library,handle)
    def restore(self,token):
        self._ready(); self.library.check(self.library.dll.fvq_sim_restore_state(self.handle,token))
    def release(self,token):
        self._ready(); self.library.check(self.library.dll.fvq_sim_release_state(self.handle,token))
    def close(self):
        if self.handle:
            self.library.dll.fvq_sim_destroy(self.handle); self.handle=P()
    def __enter__(self): return self
    def __exit__(self,*_): self.close()
    def __del__(self): self.close()

class ScalarTrace:
    """Owned immutable actual-solver trace. Geometry, hydraulics, placement and
    schedule frozen. Derivatives at baseline for initial concentrations,
    per-step inflow concentrations and a recorded vertical diffusion multiplier.
    No thermal, kinetic, geometry or SWMM parameter gradients are implied.
    States are nested [initial + completed steps][cells]. Sources are flattened
    in step order; source_ids identify repeated interval/source contributions.
    """
    def __init__(self,library,handle):
        self.library,self.handle=library,handle
        try:
            info=TraceSize();library.check(library.dll.fvq_scalar_trace_size(handle,ct.byref(info)))
            self.cells,self.steps,self.source_count,self.storage_bytes=(info.cells,info.steps,info.sources,info.storageBytes)
            self._state_count=self.cells*(self.steps+1)
            states=(D*self._state_count)();sources=(D*self.source_count)();times=(D*(self.steps+1))()
            library.check(library.dll.fvq_scalar_trace_baseline(handle,states,len(states),sources,len(sources),times,len(times)))
            self.states=self._reshape(states);self.sources=list(sources);self.times=list(times)
            self.source_ids=[library.dll.fvq_scalar_trace_source_id(handle,i).decode() for i in range(self.source_count)]
        except BaseException:
            self.close();raise
    def _ready(self):
        if not self.handle: raise RuntimeError("Scalar trace is closed")
    def _reshape(self,values):
        return [list(values[i:i+self.cells]) for i in range(0,len(values),self.cells)]
    def _action(self,method,initial,sources,diffusion):
        self._ready()
        if len(initial)!=self.cells or len(sources)!=self.source_count: raise ValueError("Trace input dimensions differ")
        a=(D*self.cells)(*initial);b=(D*self.source_count)(*sources);out=(D*self._state_count)()
        self.library.check(getattr(self.library.dll,"fvq_scalar_trace_"+method)(self.handle,a,len(a),b,len(b),diffusion,out,len(out)))
        return self._reshape(out)
    def evaluate(self,initial=None,sources=None,diffusion_multiplier=1):
        return self._action("evaluate",self.states[0] if initial is None else initial,
                            self.sources if sources is None else sources,diffusion_multiplier)
    def jvp(self,initial=None,sources=None,diffusion=0):
        return self._action("jvp",[0]*self.cells if initial is None else initial,
                            [0]*self.source_count if sources is None else sources,diffusion)
    def vjp(self,state_cotangents):
        self._ready()
        if len(state_cotangents)!=self.steps+1 or any(len(row)!=self.cells for row in state_cotangents):
            raise ValueError("Cotangents must have shape [steps+1,cells]")
        w=(D*self._state_count)(*(x for row in state_cotangents for x in row))
        a=(D*self.cells)();b=(D*self.source_count)();d=D()
        self.library.check(self.library.dll.fvq_scalar_trace_vjp(self.handle,w,len(w),a,len(a),b,len(b),ct.byref(d)))
        return {"initial":list(a),"sources":list(b),"diffusion_multiplier":d.value}
    def close(self):
        if self.handle: self.library.dll.fvq_scalar_trace_destroy(self.handle);self.handle=P()
    def __enter__(self): return self
    def __exit__(self,*_): self.close()
    def __del__(self): self.close()
