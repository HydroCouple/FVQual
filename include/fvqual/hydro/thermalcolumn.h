// SPDX-License-Identifier: MIT
#ifndef FVQUAL_HYDRO_THERMALCOLUMN_H
#define FVQUAL_HYDRO_THERMALCOLUMN_H
#include "fvqual/hydro/internalmode.h"
#include <map>
namespace FVQual::Hydro {
//! Controlled native InternalModeSolver experiment. One closed, always-wet
//! column; zero bed friction/Coriolis, no sources/reactions/settling. Outputs
//! are temperature and velocity departures from the thickness-weighted mean.
//! Density, Richardson mixing, wind shear and meteorological heat feed back.
//! This is not a horizontal-reservoir or coupled SWMM adjoint.
struct FVQual_EXPORT ThermalColumnConfiguration {
    int layers=8;
    double depth=10,area=10000,bedAreaFraction=.35,dt=30,duration=3600;
    int checkpointStride=0; //!< 0 chooses ceil(sqrt(steps)); 1 stores every step.
    InternalModeOptions options;
    Meteorology weather;
    ThermalColumnConfiguration();
};
struct ColumnPullback {
    std::vector<double> initial,parameters;
    size_t recomputedSteps=0,peakTapeNodes=0,peakBlockStates=0;
};
class FVQual_EXPORT ThermalColumnTrace {
public:
    static const std::vector<std::string> &parameterNames();
    static std::vector<double> defaultParameters(const ThermalColumnConfiguration &);
    static std::vector<double> defaultInitial(const ThermalColumnConfiguration &);
    ThermalColumnTrace(ThermalColumnConfiguration,std::vector<double> initial={},std::vector<double> parameters={});
    [[nodiscard]] std::vector<double> evaluate(std::span<const double> initial,std::span<const double> parameters) const;
    [[nodiscard]] std::vector<double> jvp(std::span<const double> initialDirection,std::span<const double> parameterDirection) const;
    [[nodiscard]] ColumnPullback vjp(std::span<const double> stateCotangents) const;
    [[nodiscard]] const std::vector<double> &initial() const {return m_initial;}
    [[nodiscard]] const std::vector<double> &parameters() const {return m_parameters;}
    [[nodiscard]] const std::vector<double> &times() const {return m_times;}
    [[nodiscard]] size_t checkpointCount() const {return m_checkpoints.size();}
    [[nodiscard]] double maximumForwardDifference() const {return m_forwardDifference;}
    [[nodiscard]] double minimumBuoyancyMargin() const {return m_buoyancyMargin;}
    [[nodiscard]] bool derivativeSupported() const {return m_supported;}
    [[nodiscard]] const ThermalColumnConfiguration &configuration() const {return m_config;}
private:
    ThermalColumnConfiguration m_config;
    std::vector<double> m_initial,m_parameters,m_times;
    std::vector<std::pair<double,double>> m_checks;
    std::map<size_t,std::vector<double>> m_checkpoints;
    double m_forwardDifference=0,m_buoyancyMargin=1;
    bool m_supported=true;
    void validateReplay(size_t,std::span<const double>) const;
};
}
#endif
