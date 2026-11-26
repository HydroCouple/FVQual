// SPDX-License-Identifier: MIT
#ifndef FVQUAL_HYDRO_PROCESSBUDGET_H
#define FVQUAL_HYDRO_PROCESSBUDGET_H

#include "fvqual/forcing/sources.h"
#include <array>
#include <vector>

namespace FVQual::Hydro {
//! Cumulative amounts since initialization. Scalars retain concentration * m3;
//! thermal surface terms are J. In/out water amounts are positive magnitudes.
struct ProcessBudget {
    struct Flow {
        double waterIn = 0, waterOut = 0;
        std::vector<double> scalarIn, scalarOut;
    };
    double initialWater = 0;
    std::vector<double> initialScalar, settling, reactions;
    std::array<Flow, static_cast<size_t>(Forcing::SourceOrigin::Count)> flows;
    // Net shortwave, incoming longwave, back radiation, latent, sensible,
    // prescribed surface flux. Positive heats water; latent is energy only.
    std::array<double, 6> surfaceHeat{};
};

//! Optional cumulative heat audit on fixed cell identities since initialize.
//! All amounts are J using rho_ref * Cp * temperature_C * volume. Interfaces
//! are positive upward. Cells follow the moving sigma grid, not fixed depths.
//! Velocities/diffusivity must be read separately; this is an energy ledger.
struct HeatTransportAudit {
    std::vector<double> initial, horizontalAdvection, verticalAdvection, sources;
    std::vector<double> horizontalDiffusion, verticalDiffusion, surface, other;
    std::vector<double> upwardAdvection, upwardDiffusion;
};

//! Cumulative per-column audit of the implicit vertical momentum operator.
//! Impulses are kg*m/s, work/dissipation are J; positive drag/viscosity terms
//! remove energy. Backward Euler's numerical dissipation is accounted for.
//! Mode reconciliation is measured separately. This is not a full mechanical
//! energy ledger for the pressure, horizontal transport or external subcycle.
//! Relative residuals include old/new absolute inventory scales, so round-off
//! in nearly unchanged velocities is not divided by a vanishing step impulse.
struct MomentumColumnAudit {
    double windImpulseX=0, windImpulseY=0, bedImpulseX=0, bedImpulseY=0;
    double verticalResidualX=0, verticalResidualY=0;
    double windWork=0, bedDissipation=0, viscosityDissipation=0;
    double implicitDissipation=0, verticalEnergyResidual=0;
    double maxRelativeImpulseResidual=0, maxRelativeEnergyResidual=0;
    double maxAbsoluteImpulseResidual=0, maxAbsoluteEnergyResidual=0;
    double modeImpulseX=0, modeImpulseY=0, modeEnergy=0;
};
}
#endif
