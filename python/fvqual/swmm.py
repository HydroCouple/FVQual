"""One-way SWMM surface-boundary exports through accepted cumulative receipts.
This is a C ABI bridge, not a native HydroCouple SWMM component or a 1D outfall
adapter. Donor exports have already been debited; no SWMM extraction is added.
"""
import ctypes as ct
import json
import math
import os
from pathlib import Path
from .engine import D, I, P

class NativeReceipt(ct.Structure):
    _fields_=[(name,D) for name in (
        "requested_m3 applied_m3 rejected_m3 requested_heat_j applied_heat_j rejected_heat_j "
        "last_applied_m3_s last_interval_seconds").split()]

class SurfaceBoundaryDonor:
    def __init__(self,library,model,output,cell,edge,species_mapping,column=0):
        self.handle=P()
        self.path=Path(library).resolve()
        self.dll=ct.CDLL(str(self.path))
        self.cell,self.edge,self.column=cell,edge,column
        self.pending=None
        self.time=0.0
        self.sequence=0
        self.output=Path(output)
        self.output.mkdir(parents=True,exist_ok=False)
        self._bind()
        self.handle=self.dll.swmm_engine_create()
        if not self.handle: raise RuntimeError("SWMM engine creation failed")
        try:
            self.check(self.dll.swmm_engine_open(self.handle,os.fsencode(Path(model).resolve()),
                os.fsencode(self.output/"donor.rpt"),os.fsencode(self.output/"donor.out"),None))
            self.check(self.dll.swmm_engine_initialize(self.handle))
            self.check(self.dll.swmm_engine_start(self.handle,0))
            flags=ct.c_uint()
            self.check(self.dll.swmm_coupling_capabilities(self.handle,ct.byref(flags)))
            if flags.value & (4|512) != (4|512):
                raise ValueError("Live bridge requires surface heat and accepted boundary receipt capabilities")
            count=I()
            self.check(self.dll.swmm_2d_species_count(self.handle,ct.byref(count)))
            self.names=[]
            for i in range(count.value):
                name=ct.create_string_buffer(512)
                self.check(self.dll.swmm_2d_species_name(self.handle,i,name,len(name)))
                self.names.append(name.value.decode())
            required=set(self.names)-{"__TEMPERATURE__","__WATER_AGE__"}
            if set(species_mapping)!=required:
                raise ValueError("Map every SWMM transported pollutant/MSX species explicitly")
            self.mapping=species_mapping
            targets=[m["target"] for m in self.mapping.values()]
            if len(targets)!=len(set(targets)):
                raise ValueError("Use distinct target IDs; aggregation is not implemented")
            for m in self.mapping.values():
                if not math.isfinite(m["kg_per_native_mass"]) or m["kg_per_native_mass"]<=0:
                    raise ValueError("Species conversion must be explicit and positive")
            self.previous=self.receipt()
        except BaseException:
            self.close()
            raise
    def _function(self,name,args,result=I):
        f=getattr(self.dll,name); f.argtypes=args; f.restype=result
    def _bind(self):
        f=self._function
        f("swmm_engine_create",[],P)
        f("swmm_engine_open",[P,ct.c_char_p,ct.c_char_p,ct.c_char_p,ct.c_char_p])
        f("swmm_engine_initialize",[P]); f("swmm_engine_start",[P,I])
        f("swmm_engine_close",[P]); f("swmm_engine_destroy",[P],None)
        f("swmm_coupling_capabilities",[P,ct.POINTER(ct.c_uint)])
        f("swmm_coupling_advance_to",[P,D,ct.POINTER(D)])
        f("swmm_2d_species_count",[P,ct.POINTER(I)])
        f("swmm_2d_species_name",[P,I,ct.c_char_p,I])
        f("swmm_2d_get_total_volume",[P,ct.POINTER(D)])
        f("swmm_2d_get_species_ledger",[P,I,I,ct.POINTER(D)])
        f("swmm_coupling_get_boundary_receipt",[P,I,I,I,ct.POINTER(NativeReceipt),
                                              ct.POINTER(D),ct.POINTER(D),I])
        f("swmm_2d_set_runtime_boundary",[P,I,I,I,D,ct.POINTER(D),I])
    def check(self,code):
        if code: raise RuntimeError(f"SWMM API error {code}")
    @property
    def volume(self):
        value=D(); self.check(self.dll.swmm_2d_get_total_volume(self.handle,ct.byref(value)))
        return value.value
    def species_inventory(self,name):
        value=D(); self.check(self.dll.swmm_2d_get_species_ledger(self.handle,self.names.index(name),0,ct.byref(value)))
        return value.value
    def prescribe_outflow(self,flow_m3_s):
        """Synthetic/reference prescription. Production uses authored boundaries."""
        if self.pending is not None: raise RuntimeError("Unaccepted export is pending")
        if not math.isfinite(flow_m3_s) or flow_m3_s<0: raise ValueError("Outward flow must be nonnegative")
        concentrations=(D*len(self.names))()
        self.check(self.dll.swmm_2d_set_runtime_boundary(self.handle,self.cell,self.edge,2,
            flow_m3_s,concentrations,len(self.names)))
    def receipt(self):
        receipt=NativeReceipt()
        requested=(D*len(self.names))(); applied=(D*len(self.names))()
        self.check(self.dll.swmm_coupling_get_boundary_receipt(self.handle,0,self.cell,self.edge,
                   ct.byref(receipt),requested,applied,len(self.names)))
        return receipt,list(applied)
    def advance(self,end_s):
        if self.pending is not None: raise RuntimeError("Accept or retry the pending export before advancing SWMM")
        if not math.isfinite(end_s) or end_s<=self.time: raise ValueError("Invalid donor exchange time")
        actual=D()
        self.check(self.dll.swmm_coupling_advance_to(self.handle,end_s,ct.byref(actual)))
        if abs(actual.value-end_s)>1e-8: raise RuntimeError("SWMM did not land on exchange time")
        current=self.receipt()
        water=-(current[0].applied_m3-self.previous[0].applied_m3)
        heat=-(current[0].applied_heat_j-self.previous[0].applied_heat_j)
        if water < -1e-12:
            raise ValueError("Boundary reversal requires a feedback/extraction protocol")
        if water < 0: water=0.0  # only floating-point cancellation at zero
        masses={}
        for native,m in self.mapping.items():
            index=self.names.index(native)
            mass=-(current[1][index]-self.previous[1][index])*m["kg_per_native_mass"]
            if mass < -1e-12: raise ValueError("Negative species export")
            masses[m["target"]]=max(0,mass)
        frame={"id":f"swmm-surface/{self.cell}/{self.edge}/{self.sequence}",
               "start_s":self.time,"end_s":end_s,
               "inflows":[{"column":self.column,"volume_m3":water,"heat_j":heat,
                            "placement":"uniform","species_mass_kg":masses}]}
        # Write the accepted donor debit before the receiver can fail. The
        # donor has no rollback contract; a failed receiver must stop the run.
        (self.output/f"export_{self.sequence:05d}.json").write_text(json.dumps(frame,indent=2,allow_nan=False)+"\n")
        self.pending=(frame,current)
        return frame
    def accept(self,transfer_id):
        if self.pending is None or transfer_id!=self.pending[0]["id"]:
            raise ValueError("No matching pending export")
        frame,current=self.pending
        self.previous=current
        self.time=frame["end_s"]
        self.sequence+=1
        self.pending=None
    def close(self):
        if self.handle:
            self.dll.swmm_engine_close(self.handle); self.dll.swmm_engine_destroy(self.handle); self.handle=P()
    def __enter__(self): return self
    def __exit__(self,*_): self.close()
    def __del__(self): self.close()
