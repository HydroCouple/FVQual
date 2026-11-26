// SPDX-License-Identifier: MIT
#include "fvqual/hydro/thermalcolumn.h"
#include <gtest/gtest.h>
#include <cmath>
#include <numeric>
using namespace FVQual::Hydro;
namespace {
double dot(const std::vector<double> &a,const std::vector<double> &b) {return std::inner_product(a.begin(),a.end(),b.begin(),0.0);}
}
TEST(ThermalColumn, NativeHeatShearAndGeometryDirectionalSweep) {
    ThermalColumnConfiguration c;c.duration=300;c.dt=10;c.layers=6;
    ThermalColumnTrace trace(c);ASSERT_TRUE(trace.derivativeSupported());EXPECT_LT(trace.maximumForwardDifference(),1e-11);
    std::vector<double> dx(trace.initial().size()),dp(trace.parameters().size(),0);
    for(size_t i=0;i<dx.size();++i)dx[i]=i<6?.2*std::cos(i):.002*std::sin(i);
    dp[0]=.2;dp[2]=.3;dp[4]=10;dp[6]=.001;dp[8]=2e-6;dp[10]=.1;dp[11]=.03;dp[12]=.02;dp[13]=20;dp[14]=.001;dp[16]=1e-4;dp[17]=.04;
    const auto tangent=trace.jvp(dx,dp);
    for(double epsilon:{.01,.001,.0001,.00001}) {
        std::vector<std::vector<double>> runs;
        for(double sign:{1.,-1.}) {
            auto x=trace.initial(),p=trace.parameters();
            for(size_t i=0;i<x.size();++i)x[i]+=sign*epsilon*dx[i];
            for(size_t i=0;i<p.size();++i)p[i]+=sign*epsilon*dp[i];
            runs.push_back(trace.evaluate(x,p));
        }
        double difference=0;for(size_t i=0;i<tangent.size();++i)difference=std::max(difference,std::abs((runs[0][i]-runs[1][i])/(2*epsilon)-tangent[i]));
        EXPECT_LT(difference,2e-7+epsilon*epsilon*.01) << epsilon;
    }
    std::vector<double> w(tangent.size());for(size_t i=0;i<w.size();++i)w[i]=std::sin(i*.39);
    const auto bar=trace.vjp(w);
    EXPECT_NEAR(dot(w,tangent),dot(dx,bar.initial)+dot(dp,bar.parameters),2e-12);
    EXPECT_EQ(bar.recomputedSteps,30u);EXPECT_GT(bar.peakTapeNodes,100u);
}
TEST(ThermalColumn, FullyMixedHeatHasIndependentDepthCapacityAndEnergyDerivatives) {
    ThermalColumnConfiguration c;c.layers=1;c.duration=100;c.dt=10;c.depth=4;c.area=100;c.bedAreaFraction=.2;
    c.weather.albedo=0;c.weather.waterEmissivity=0;c.weather.longwaveReflection=1;c.weather.psychrometric=0;
    auto p=ThermalColumnTrace::defaultParameters(c);p[0]=0;p[4]=200;p[14]=p[15]=0;
    ThermalColumnTrace trace(c,{10,0,0},p);const auto values=trace.evaluate(trace.initial(),p);
    const double volume=240,rho=1000,cp=p[13],heat=200*100*100;
    EXPECT_NEAR(values[values.size()-3],10+heat/(rho*cp*volume),2e-13);
    std::vector<double> dp(p.size(),0);dp[11]=1;auto gradient=trace.jvp(std::vector<double>(3,0),dp);
    EXPECT_NEAR(gradient[gradient.size()-3],-heat/(rho*cp*volume),2e-13);
    // E=rho*cp*V*T; depth dependence includes changing storage, not just dT.
    const double dEnergy=rho*cp*volume*(gradient[gradient.size()-3]+values[values.size()-3]);
    EXPECT_NEAR(dEnergy,rho*cp*volume*10,1e-5);
    dp.assign(p.size(),0);dp[17]=1;gradient=trace.jvp(std::vector<double>(3,0),dp);
    for(double value:gradient)EXPECT_NEAR(value,0,1e-13); // area and volume scale together
}
TEST(ThermalColumn, CheckpointPoliciesAgreeAndRecordingDoesNotChangePhysics) {
    ThermalColumnConfiguration c;c.duration=240;c.dt=10;c.layers=4;c.checkpointStride=1;ThermalColumnTrace full(c);
    auto baseline=full.evaluate(full.initial(),full.parameters());std::vector<double> w(baseline.size(),0);
    for(size_t i=0;i<w.size();++i)w[i]=std::cos(i*.1);
    const auto expected=full.vjp(w);
    c.checkpointStride=7;ThermalColumnTrace blocked(c);const auto actual=blocked.vjp(w);
    EXPECT_LT(blocked.checkpointCount(),full.checkpointCount());EXPECT_EQ(blocked.evaluate(blocked.initial(),blocked.parameters()),baseline);
    for(size_t i=0;i<actual.initial.size();++i)EXPECT_NEAR(actual.initial[i],expected.initial[i],2e-12);
    for(size_t i=0;i<actual.parameters.size();++i)EXPECT_NEAR(actual.parameters[i],expected.parameters[i],2e-10);
    EXPECT_LE(actual.peakBlockStates,8u);EXPECT_EQ(actual.recomputedSteps,24u);
    EXPECT_THROW((void)blocked.vjp(std::vector<double>{1}),std::invalid_argument);
    EXPECT_EQ(blocked.evaluate(blocked.initial(),blocked.parameters()),baseline);
}
TEST(ThermalColumn, UniformDensityBranchRefusesGradientsButPermitsForwardEnsembles) {
    ThermalColumnConfiguration c;c.layers=4;c.duration=10;c.dt=10;
    std::vector<double> x(12,0);std::fill(x.begin(),x.begin()+4,20);
    ThermalColumnTrace trace(c,x);EXPECT_FALSE(trace.derivativeSupported());
    EXPECT_THROW((void)trace.jvp(std::vector<double>(12,0),std::vector<double>(18,0)),std::invalid_argument);
    EXPECT_FALSE(trace.evaluate(x,trace.parameters()).empty());
    c.options.manningN=.03;EXPECT_THROW((void)ThermalColumnTrace(c),std::invalid_argument);
}
