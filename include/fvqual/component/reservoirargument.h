// SPDX-License-Identifier: MIT
#ifndef FVQUAL_RESERVOIRARGUMENT_H
#define FVQUAL_RESERVOIRARGUMENT_H
#include "hydrocouplesdk/data/argument1d.h"
#include "hydrocouplesdk/data/exchangeitems1d.h"
#include "hydrocouplespatial.h"
namespace FVQual::Component {
// An optional named override distinguishes an omitted value from its default,
// preserving compositions authored with the legacy geometry arrays.
template<class T> class ReservoirArgument : public HydroCouple::SDK::Argument1D<T> {
public:
    using Base = HydroCouple::SDK::Argument1D<T>;
    using Base::Base;
    bool supplied = false;
    bool initialize(const std::string &v, HydroCouple::IArgument::ArgumentInputType t,
                    std::string &m) override {
        if (!Base::initialize(v,t,m)) return false;
        supplied = true; return true;
    }
    bool initialize(const HydroCouple::IComponentDataItem &v, std::string &m) override {
        if (!Base::initialize(v,m)) return false;
        supplied = true; return true;
    }
};
class ColumnOutput : public HydroCouple::SDK::Output1DDouble,
                     public virtual HydroCouple::Spatial::IPolyhedralSurfaceComponentDataItem {
public:
    ColumnOutput(std::string_view id, HydroCouple::SDK::Dimension *dim, int n,
                 HydroCouple::SDK::ValueDefinition *q, HydroCouple::IModelComponent *mc,
                 HydroCouple::Spatial::IPolyhedralSurface *surface)
        : HydroCouple::SDK::Output1DDouble(id,dim,n,q,mc), m_surface(surface), m_dimension(dim) {}
    HydroCouple::Spatial::IPolyhedralSurface *polyhedralSurface() const override {return m_surface;}
    HydroCouple::Spatial::MeshLocation location() const override {return HydroCouple::Spatial::MeshLocation::Face;}
    HydroCouple::Spatial::SpatialDataType meshDataType() const override {return HydroCouple::Spatial::SpatialDataType::Scalar;}
    HydroCouple::Spatial::VectorBasis vectorBasis() const override {return HydroCouple::Spatial::VectorBasis::Unknown;}
    HydroCouple::IDimension *entityDimension() const override {return m_dimension;}
private:
    HydroCouple::Spatial::IPolyhedralSurface *m_surface;
    HydroCouple::SDK::Dimension *m_dimension;
};
}
#endif
