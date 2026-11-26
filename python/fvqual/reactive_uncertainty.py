"""Local conditional uncertainty and observation-design diagnostics.

Modes are columns of a covariance factor over initial species and kinetic
constants. Hydraulic/thermal forcing and structural discrepancy are prescribed.
NumPy is optional and is required only for the Gaussian/SVD diagnostic.
"""
import math


def local_covariance(trace,modes,observations):
    if not modes:raise ValueError("At least one covariance-factor mode is required")
    columns=[observations.jvp(trace.jvp(**mode)) for mode in modes]
    n=len(observations.rows)
    covariance=[[math.fsum(col[i]*col[j] for col in columns) for j in range(n)] for i in range(n)]
    return {"mean":observations.apply(trace.states),"mode_sensitivities":columns,"covariance":covariance,
            "standard_deviation":[math.sqrt(max(0,covariance[i][i])) for i in range(n)],
            "scope":"local reactive-species parameter covariance conditional on prescribed hydraulics, heat, geometry and inputs"}


def gaussian_diagnostics(local,observation_sd,observed=None):
    """Whitened sensitivity SVD and local Gaussian update in prior-mode coordinates.

    Independent Gaussian observation errors are an explicit assumption. The
    posterior is conditional on this local linearization and the supplied data;
    it is not a calibrated predictive interval or a field-validation result.
    Without observations, covariance describes an illustrative survey design.
    """
    import numpy as np
    j=np.asarray(local["mode_sensitivities"],dtype=float).T
    mean=np.asarray(local["mean"],dtype=float);sd=np.asarray(observation_sd,dtype=float)
    if j.ndim!=2 or j.shape[0]!=mean.size or sd.shape!=mean.shape or np.any(sd<=0) or not np.all(np.isfinite(j)) or not np.all(np.isfinite(mean)) or not np.all(np.isfinite(sd)):
        raise ValueError("Finite sensitivity columns and positive observation SDs are required")
    a=j/sd[:,None];singular=np.linalg.svd(a,compute_uv=False)
    tolerance=max(a.shape)*np.finfo(float).eps*(singular[0] if len(singular) else 0)
    rank=int(np.sum(singular>tolerance));cov=np.linalg.solve(np.eye(j.shape[1])+a.T@a,np.eye(j.shape[1]))
    norms=np.linalg.norm(a,axis=0);denominator=norms[:,None]*norms[None,:]
    cosine=np.divide(a.T@a,denominator,out=np.zeros((len(norms),len(norms))),where=denominator>0)
    latent=np.zeros(j.shape[1])
    if observed is not None:
        data=np.asarray(observed,dtype=float)
        if data.shape!=mean.shape or not np.all(np.isfinite(data)):raise ValueError("Observed values must match the observation vector")
        latent=cov@a.T@((data-mean)/sd)
    parameter_cov=j@cov@j.T
    return {"singular_values":singular.tolist(),"numerical_rank":rank,"mode_count":j.shape[1],
            "condition_number":float(singular[0]/singular[-1]) if rank==j.shape[1] else None,
            "sensitivity_cosines":cosine.tolist(),"posterior_mode_covariance":cov.tolist(),"posterior_mode_mean":latent.tolist(),
            "posterior_observation_mean":(mean+j@latent).tolist(),"posterior_parameter_sd":np.sqrt(np.maximum(0,np.diag(parameter_cov))).tolist(),
            "observation_sd":sd.tolist(),"data_used":observed is not None,
            "scope":"local conditional Gaussian calculation in prior covariance-factor coordinates; independent stated observation noise; structural discrepancy excluded"}
