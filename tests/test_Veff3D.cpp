#include <gtest/gtest.h>

#include <complex>

#include "testparams.h"
#include "effpot.h"


/* Test params in 3D theory (without A0). 
These should correspond to inputs 
MH2 = 600, sinTheta = 0.12, a2 = 4.0, b3 = -75, b4 = 0.5
using 1-loop MS-bar matching and NNLO DR matching at T = 118. 
DR matching scale is the usual ~7.055T */
ParameterMap SetTestParameters3D() {

    ParameterMap params3D;
    params3D["RGScale"] = 118.0;
    params3D["g1sq"] = 15.2742;
    params3D["g2sq"] = 46.6236;
    params3D["msqPhi"] = -834.085;
    params3D["lambda"] = 29.9012;
    params3D["b1"] = -1.05681e+6;
    params3D["b2"] = 233676.;
    params3D["b3"] = -414.537;
    params3D["b4"] = 162.743;
    params3D["a1"] = 2862.32;
    params3D["a2"] = 686.951;

    params3D["c05"] = 0.0165248;
    params3D["c23"] = 0.119285;
    params3D["c41"] = 0.100293;
    params3D["c60"] = 0.00308239;
    params3D["c06"] = 0.00292394;
    params3D["c42"] = 0.0343714;
    params3D["c24"] = 0.0316515;

    return params3D;
}


TEST(Veff3DTest, TreeLevel) {

    ParameterMap params3D = SetTestParameters3D();

    // background fields
    double v = 0.00451;
    double x = -3.42;

    EffPot<std::complex<double>> effPot(params3D);

    // Will not include dim-6 

    std::complex<double> V0 = effPot.EvaluatePotentialAsDouble(v, x, ELoopOrderVeff::tree, false);

    EXPECT_DOUBLE_EQ(V0.real(), 4991967.6045201151);
    EXPECT_DOUBLE_EQ(V0.imag(), 0.0);
}


TEST(Veff3DTest, OneLoop) {

    ParameterMap params3D = SetTestParameters3D();

    // background fields
    double v = 0.00451;
    double x = -3.42;

    EffPot<std::complex<double>> effPot(params3D);

    // Will not include dim-6 

    std::complex<double> V0 = effPot.EvaluatePotentialAsDouble(v, x, ELoopOrderVeff::loop1, false);

    EXPECT_DOUBLE_EQ(V0.real(), 1829768.8946448057);
    EXPECT_DOUBLE_EQ(V0.imag(), 7510.8544202073263);
}


TEST(Veff3DTest, TwoLoop) {

    ParameterMap params3D = SetTestParameters3D();

    // background fields
    double v = 0.00451;
    double x = -3.42;

    EffPot<std::complex<double>> effPot(params3D);

    // Will not include dim-6 

    std::complex<double> V0 = effPot.EvaluatePotentialAsDouble(v, x, ELoopOrderVeff::loop2, false);

    EXPECT_DOUBLE_EQ(V0.real(), 2036804.0742510404);
    EXPECT_DOUBLE_EQ(V0.imag(), 99606.765149759405);
}
