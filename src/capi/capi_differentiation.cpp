// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "capi_internal.hpp"
#include "fvqual/capi/fvqual_differentiation.h"
#include <array>
#include <cmath>
#include <stdexcept>
struct FVQ_ScalarTrace_ {
    FVQual::Transport::ScalarTrajectory trajectory;
    std::vector<std::string> ids;
};
namespace {
void array(const void *p,int64_t count,size_t expected) {
    if (count<0 || static_cast<size_t>(count)!=expected || (!p && count))
        throw std::invalid_argument("scalar trace array is null or has the wrong length");
}
const FVQual::Transport::ScalarTrajectory &trajectory(const FVQ_ScalarTrace *t) {
    if (!t) throw std::invalid_argument("scalar trace is null");return t->trajectory;
}
void dimensions(const FVQ_ScalarTrace *t,const void *initial,int64_t cells,
    const void *sources,int64_t sourceCount,const void *states,int64_t stateCount) {
    const auto &r=trajectory(t);array(initial,cells,r.initial.size());
    array(sources,sourceCount,t->ids.size());array(states,stateCount,(r.count+1)*r.initial.size());
}
struct Evaluation { double value; std::array<double,6> gradient{}; };
Evaluation evaluate(FVQ_Kernel kind, const double *x, int n) {
    const int expected = kind == FVQ_KERNEL_DECAY ? 3 :
        kind == FVQ_KERNEL_COLUMN_HEAT ? 6 : kind == FVQ_KERNEL_MIX ? 4 : -1;
    if (!x || n != expected) throw std::invalid_argument("invalid kernel or input count");
    for (int i = 0; i < n; ++i)
        if (!std::isfinite(x[i])) throw std::invalid_argument("kernel inputs must be finite");
    Evaluation result;
    if (kind == FVQ_KERNEL_DECAY) {
        if (x[1] < 0 || x[2] < 0) throw std::invalid_argument("negative decay rate or duration");
        const double e = std::exp(-x[1]*x[2]);
        result.value = x[0]*e;
        result.gradient[0] = e;
        result.gradient[1] = -x[2]*result.value;
        result.gradient[2] = -x[1]*result.value;
    } else if (kind == FVQ_KERNEL_COLUMN_HEAT) {
        if (x[2] < 0 || x[3] <= 0 || x[4] <= 0 || x[5] <= 0)
            throw std::invalid_argument("heat kernel requires positive depth, density and heat capacity");
        const double denominator = x[3]*x[4]*x[5];
        const double delta = x[1]*x[2]/denominator;
        result.value = x[0]+delta;
        result.gradient[0] = 1;
        result.gradient[1] = x[2]/denominator;
        result.gradient[2] = x[1]/denominator;
        result.gradient[3] = -delta/x[3];
        result.gradient[4] = -delta/x[4];
        result.gradient[5] = -delta/x[5];
    } else {
        const double volume = x[0]+x[2];
        if (x[0] < 0 || x[2] < 0 || !(volume > 0))
            throw std::invalid_argument("mix kernel requires nonnegative volumes and positive total");
        result.value = (x[0]*x[1]+x[3])/volume;
        result.gradient[0] = (x[1]-result.value)/volume;
        result.gradient[1] = x[0]/volume;
        result.gradient[2] = -result.value/volume;
        result.gradient[3] = 1/volume;
    }
    if (!std::isfinite(result.value)) throw std::invalid_argument("kernel value overflows");
    for (int i = 0; i < n; ++i)
        if (!std::isfinite(result.gradient[i])) throw std::invalid_argument("kernel derivative overflows");
    return result;
}
template<class F> FVQ_ErrorCode checked(F f) {
    return FVQual::CApi::guard([&]() -> FVQ_ErrorCode {
        try { f(); return FVQual::CApi::succeed(); }
        catch (const std::invalid_argument &e) { return FVQual::CApi::fail(FVQ_ERR_INVALID_ARGUMENT, e.what()); }
    });
}
}
extern "C" {
FVQ_ErrorCode fvq_sim_begin_scalar_trace(FVQ_Simulation *s,int species,int64_t limit) {
    return checked([&] {
        if (!s || !s->initialized || limit<1 || limit>10000)
            throw std::invalid_argument("scalar recording requires initialized simulation and step limit 1..10000");
        std::string message;
        if (!s->solver.beginScalarTrace(species,static_cast<size_t>(limit),message)) throw std::invalid_argument(message);
    });
}
FVQ_ErrorCode fvq_sim_end_scalar_trace(FVQ_Simulation *s) {
    return checked([&] {
        if (!s || !s->initialized) throw std::invalid_argument("initialized simulation required");
        s->solver.endScalarTrace();
    });
}
FVQ_ErrorCode fvq_sim_capture_scalar_trace(const FVQ_Simulation *s,FVQ_ScalarTrace **out) {
    return checked([&] {
        if (!s || !s->initialized || !out || s->solver.scalarTrace().initial.empty())
            throw std::invalid_argument("begin scalar recording before capturing an owned trace");
        auto t=std::make_unique<FVQ_ScalarTrace>();t->trajectory=s->solver.scalarTrace();
        for (const auto *step:t->trajectory.steps())
            t->ids.insert(t->ids.end(),step->sourceIds.begin(),step->sourceIds.end());
        *out=t.release();
    });
}
void fvq_scalar_trace_destroy(FVQ_ScalarTrace *t) { delete t; }
FVQ_ErrorCode fvq_scalar_trace_size(const FVQ_ScalarTrace *t,FVQ_TraceSize *out) {
    return checked([&] {
        const auto &r=trajectory(t);if (!out) throw std::invalid_argument("trace size output is null");
        *out={static_cast<int64_t>(r.initial.size()),static_cast<int64_t>(r.count),
            static_cast<int64_t>(t->ids.size()),static_cast<int64_t>(r.bytes())};
    });
}
FVQ_ErrorCode fvq_scalar_trace_baseline(const FVQ_ScalarTrace *t,double *states,int64_t stateCount,
    double *sources,int64_t sourceCount,double *times,int64_t timeCount) {
    return checked([&] {
        const auto &r=trajectory(t);
        array(states,stateCount,(r.count+1)*r.initial.size());array(sources,sourceCount,t->ids.size());
        array(times,timeCount,r.count+1);
        const auto sourceValues=r.sourceValues();const auto history=r.steps();
        std::copy(r.initial.begin(),r.initial.end(),states);times[0]=r.startTime;
        size_t offset=r.initial.size();
        for (size_t k=0;k<history.size();++k) {
            std::copy(history[k]->after.begin(),history[k]->after.end(),states+offset);
            offset+=r.initial.size();times[k+1]=history[k]->time+history[k]->dt;
        }
        std::copy(sourceValues.begin(),sourceValues.end(),sources);
    });
}
const char *fvq_scalar_trace_source_id(const FVQ_ScalarTrace *t,int64_t source) {
    if (!t || source<0 || static_cast<size_t>(source)>=t->ids.size()) {
        FVQual::CApi::fail(FVQ_ERR_INVALID_ARGUMENT,"invalid trace source index");return nullptr;
    }
    FVQual::CApi::succeed();return t->ids[source].c_str();
}
FVQ_ErrorCode fvq_scalar_trace_evaluate(const FVQ_ScalarTrace *t,const double *initial,int64_t cells,
    const double *sources,int64_t sourceCount,double scale,double *states,int64_t stateCount) {
    return checked([&] {
        dimensions(t,initial,cells,sources,sourceCount,states,stateCount);
        const auto result=t->trajectory.evaluate({initial,static_cast<size_t>(cells)},
            {sources,static_cast<size_t>(sourceCount)},scale);
        std::copy(result.begin(),result.end(),states);
    });
}
FVQ_ErrorCode fvq_scalar_trace_jvp(const FVQ_ScalarTrace *t,const double *initial,int64_t cells,
    const double *sources,int64_t sourceCount,double direction,double *states,int64_t stateCount) {
    return checked([&] {
        dimensions(t,initial,cells,sources,sourceCount,states,stateCount);
        const auto result=t->trajectory.jvp({initial,static_cast<size_t>(cells)},
            {sources,static_cast<size_t>(sourceCount)},direction);
        std::copy(result.begin(),result.end(),states);
    });
}
FVQ_ErrorCode fvq_scalar_trace_vjp(const FVQ_ScalarTrace *t,const double *states,int64_t stateCount,
    double *initial,int64_t cells,double *sources,int64_t sourceCount,double *scale) {
    return checked([&] {
        dimensions(t,initial,cells,sources,sourceCount,states,stateCount);
        if (!scale) throw std::invalid_argument("diffusion cotangent output is null");
        const auto result=t->trajectory.vjp({states,static_cast<size_t>(stateCount)});
        std::copy(result.initial.begin(),result.initial.end(),initial);
        std::copy(result.sources.begin(),result.sources.end(),sources);*scale=result.diffusionMultiplier;
    });
}
FVQ_ErrorCode fvq_kernel_value(FVQ_Kernel k, const double *x, int n, double *value) {
    if (!x || !value) return FVQual::CApi::fail(FVQ_ERR_NULL_ARGUMENT, "kernel input/output is null");
    return checked([&] { const auto e = evaluate(k,x,n); *value = e.value; });
}
FVQ_ErrorCode fvq_kernel_jvp(FVQ_Kernel k, const double *x, const double *v, int n, double *product) {
    if (!x || !v || !product) return FVQual::CApi::fail(FVQ_ERR_NULL_ARGUMENT, "kernel direction/output is null");
    return checked([&] {
        const auto e = evaluate(k,x,n);
        double result = 0;
        for (int i = 0; i < n; ++i) {
            if (!std::isfinite(v[i])) throw std::invalid_argument("nonfinite tangent");
            result += e.gradient[i]*v[i];
        }
        if (!std::isfinite(result)) throw std::invalid_argument("tangent product overflows");
        *product = result;
    });
}
FVQ_ErrorCode fvq_kernel_vjp(FVQ_Kernel k, const double *x, int n, double w, double *product) {
    if (!x || !product) return FVQual::CApi::fail(FVQ_ERR_NULL_ARGUMENT, "kernel cotangent/output is null");
    return checked([&] {
        if (!std::isfinite(w)) throw std::invalid_argument("nonfinite cotangent");
        auto e = evaluate(k,x,n);
        for (int i = 0; i < n; ++i) {
            e.gradient[i] *= w;
            if (!std::isfinite(e.gradient[i])) throw std::invalid_argument("cotangent product overflows");
        }
        for (int i = 0; i < n; ++i) product[i] = e.gradient[i];
    });
}
}
