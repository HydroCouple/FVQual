// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/hydro/internalmode.h"
#include "fvqual/parallel.h"
#include <algorithm>

namespace FVQual::Hydro {
void InternalModeSolver::applyHorizontalDiffusion(double dt, std::span<double> values)
{
    const auto &mesh = *m_mesh;
    const auto &h = mesh.horizontal();
    const auto &offsets = h.faceEdgeOffsets(), &edges = h.faceEdges();
    const int layers = mesh.layerCount();
    const int64_t columns = mesh.columnCount();
    auto neighbour = [&](int64_t c, int64_t e) {
        return h.edgeLeftFace()[e] == c ? h.edgeRightFace()[e] : h.edgeLeftFace()[e];
    };
    auto conductance = [&](int64_t e, size_t a, size_t b) {
        return m_options.horizontalDiffusivity * h.edgeLength()[e] *
            .5 * (mesh.thickness()[a] + mesh.thickness()[b]) / h.edgeNormalDistance()[e];
    };
    // The low-order update is a convex average at this limit. Subcycling
    // protects direct step() callers as well as the adaptive advance() path.
    double stable = dt;
    for (int64_t c=0;c<columns;++c) for (int k=0;k<layers;++k) {
        const size_t a=mesh.cell(c,k); double sum=0;
        for (auto i=offsets[c];i<offsets[c+1];++i) {
            const auto e=edges[i], other=neighbour(c,e);
            if (other!=Mesh::kNoFace) sum+=conductance(e,a,mesh.cell(other,k));
        }
        if (sum>0) stable=std::min(stable,.45*mesh.cellVolume()[a]/sum);
    }
    double remaining=dt;
    while (remaining>0) {
        const double step=std::min(remaining,stable);
        std::copy(values.begin(),values.end(),m_scalarWork.begin());
        // Weighted least squares is exact for linear fields on a conforming
        // stencil. Closed boundaries supply the homogeneous normal flux.
        FVQual_PARALLEL_FOR
        for (int64_t c=0;c<columns;++c) for (int k=0;k<layers;++k) {
            const size_t a=mesh.cell(c,k);
            double xx=0,xy=0,yy=0,bx=0,by=0;
            for (auto i=offsets[c];i<offsets[c+1];++i) {
                const auto e=edges[i], other=neighbour(c,e);
                double dx=h.edgeNormalX()[e],dy=h.edgeNormalY()[e],delta=0,weight=1;
                if (other!=Mesh::kNoFace) {
                    dx=h.faceCentroidX()[other]-h.faceCentroidX()[c];
                    dy=h.faceCentroidY()[other]-h.faceCentroidY()[c];
                    weight=1/(dx*dx+dy*dy);
                    delta=m_scalarWork[mesh.cell(other,k)]-m_scalarWork[a];
                }
                xx+=weight*dx*dx; xy+=weight*dx*dy; yy+=weight*dy*dy;
                bx+=weight*dx*delta; by+=weight*dy*delta;
            }
            const double determinant=xx*yy-xy*xy;
            if (determinant>1e-12*(xx+yy)*(xx+yy)) {
                m_diffusionGradX[a]=(yy*bx-xy*by)/determinant;
                m_diffusionGradY[a]=(xx*by-xy*bx)/determinant;
            } else {
                m_diffusionGradX[a]=xx+yy>0?bx/(xx+yy):0;
                m_diffusionGradY[a]=xx+yy>0?by/(xx+yy):0;
            }
        }
        auto correction = [&](int64_t c,int64_t other,int64_t e,size_t a,size_t b) {
            const double direction=h.edgeLeftFace()[e]==c?1:-1;
            const double dx=h.faceCentroidX()[other]-h.faceCentroidX()[c] -
                direction*h.edgeNormalX()[e]*h.edgeNormalDistance()[e];
            const double dy=h.faceCentroidY()[other]-h.faceCentroidY()[c] -
                direction*h.edgeNormalY()[e]*h.edgeNormalDistance()[e];
            return -step*conductance(e,a,b)*.5*((m_diffusionGradX[a]+m_diffusionGradX[b])*dx +
                                               (m_diffusionGradY[a]+m_diffusionGradY[b])*dy);
        };
        // Limit the correction as an equal-and-opposite mass transfer.
        // Clipping concentrations afterwards would destroy conservation.
        FVQual_PARALLEL_FOR
        for (int64_t c=0;c<columns;++c) for (int k=0;k<layers;++k) {
            const size_t a=mesh.cell(c,k); const double old=m_scalarWork[a];
            double inflow=0,positive=0,negative=0,lo=old,hi=old;
            for (auto i=offsets[c];i<offsets[c+1];++i) {
                const auto e=edges[i], other=neighbour(c,e); if (other==Mesh::kNoFace) continue;
                const size_t b=mesh.cell(other,k);
                lo=std::min(lo,m_scalarWork[b]); hi=std::max(hi,m_scalarWork[b]);
                inflow+=conductance(e,a,b)*(m_scalarWork[b]-old);
                const double anti=correction(c,other,e,a,b);
                positive+=std::max(anti,0.0); negative+=std::max(-anti,0.0);
            }
            const double volume=mesh.cellVolume()[a];
            const double low=old+step*inflow/volume;
            m_diffusionLow[a]=low;
            m_diffusionPlus[a]=positive>0?std::clamp((hi-low)*volume/positive,0.0,1.0):1;
            m_diffusionMinus[a]=negative>0?std::clamp((low-lo)*volume/negative,0.0,1.0):1;
        }
        FVQual_PARALLEL_FOR
        for (int64_t c=0;c<columns;++c) for (int k=0;k<layers;++k) {
            const size_t a=mesh.cell(c,k); double inflow=0;
            for (auto i=offsets[c];i<offsets[c+1];++i) {
                const auto e=edges[i], other=neighbour(c,e); if (other==Mesh::kNoFace) continue;
                const size_t b=mesh.cell(other,k);
                const double anti=correction(c,other,e,a,b);
                const double fraction=anti>0?std::min(m_diffusionPlus[a],m_diffusionMinus[b]):
                    std::min(m_diffusionMinus[a],m_diffusionPlus[b]);
                inflow+=fraction*anti;
            }
            values[a]=m_diffusionLow[a]+inflow/mesh.cellVolume()[a];
        }
        remaining=step==remaining?0:remaining-step;
    }
}
} // namespace FVQual::Hydro
