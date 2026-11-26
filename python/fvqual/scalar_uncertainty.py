"""Local uncertainty on qualified passive-scalar traces; no calibrated posterior.

Modes are columns of an input covariance factor L, including any prescribed
correlations. Hydraulics/geometry and their observation weights are frozen.
Observation noise, geometry uncertainty and structural discrepancy are separate
inputs to a scientific assessment, not silently included in these bands.
"""
import math


class FixedObservation:
    """Sparse linear observations of [steps+1,cells] with scale/offset.

    Each row contains (step, cell, weight). Repeated terms/fan-out accumulate.
    Use fixed depth/time interpolation, fixed conservative mapping or inventory
    weights. Moving observation geometry requires additional derivatives.
    """
    def __init__(self,shape,rows,scale=1,offset=0):
        self.shape=tuple(shape)
        if len(self.shape)!=2 or any(not isinstance(n,int) or n<1 for n in self.shape):
            raise ValueError("Observation shape must be two positive dimensions")
        if not math.isfinite(scale) or not math.isfinite(offset): raise ValueError("Nonfinite unit conversion")
        self.scale,self.offset=scale,offset
        self.rows=[]
        for row in rows:
            terms=[]
            for step,cell,weight in row:
                if not isinstance(step,int) or not isinstance(cell,int) or not 0<=step<shape[0] or not 0<=cell<shape[1] or not math.isfinite(weight):
                    raise ValueError("Observation index or weight invalid")
                terms.append((step,cell,weight))
            if not terms: raise ValueError("Empty observation row")
            self.rows.append(terms)
        if not self.rows: raise ValueError("No observations")
    def _sample(self,states,offset):
        if len(states)!=self.shape[0] or any(len(row)!=self.shape[1] for row in states):
            raise ValueError("Observation state dimensions differ")
        if any(not math.isfinite(x) for row in states for x in row): raise ValueError("Nonfinite observation state")
        result=[offset+self.scale*math.fsum(states[t][c]*w for t,c,w in row) for row in self.rows]
        if any(not math.isfinite(x) for x in result): raise ValueError("Observation overflows")
        return result
    def apply(self,states): return self._sample(states,self.offset)
    def jvp(self,tangents): return self._sample(tangents,0)
    def vjp(self,cotangents):
        if len(cotangents)!=len(self.rows) or any(not math.isfinite(x) for x in cotangents):
            raise ValueError("Invalid observation cotangents")
        result=[[0.]*self.shape[1] for _ in range(self.shape[0])]
        for row,cotangent in zip(self.rows,cotangents):
            for t,c,w in row: result[t][c]+=self.scale*w*cotangent
        if any(not math.isfinite(x) for row in result for x in row): raise ValueError("Observation pullback overflows")
        return result


def local_covariance(trace,modes,observations):
    """Cy=J L L^T J^T. Modes contain initial/sources/diffusion directions."""
    if not modes: raise ValueError("At least one uncertainty mode is required")
    columns=[observations.jvp(trace.jvp(**mode)) for mode in modes]
    n=len(observations.rows)
    covariance=[[math.fsum(col[i]*col[j] for col in columns) for j in range(n)] for i in range(n)]
    return {"mean":observations.apply(trace.states),"covariance":covariance,
            "standard_deviation":[math.sqrt(max(0,covariance[i][i])) for i in range(n)],
            "mode_sensitivities":columns,
            "scope":"local passive-scalar parameter covariance; fixed hydraulic trajectory and observation weights"}
