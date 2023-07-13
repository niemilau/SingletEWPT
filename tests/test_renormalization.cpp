#include <gtest/gtest.h>

#include "testparams.h"
#include "renormalization.h"

/* This file checks that the T=0 matching between inputs and renormalized parameters is OK */


// Check that experimental constants are 'correct', ie. match the values used in original benchmarks
TEST(ExperimentalInputTest, PhysicalConstants) {

    double relativeError = 1e-15;

    // Expected values
    const double MH = 125.10;
    const double MW = 80.379;
    const double MZ = 91.1876;
    const double Mt = 172.76;

    // Strong fine-structure constant in MS-bar at Z-pole
    const double alphaQCD = 0.1181;
    // QCD coupling constant squared
    const double gs2 = 4.0*PI*alphaQCD;
    // Fermi constant (GeV^-2)
    const double Gf = 1.1663787e-5;

    // 'Tree-level' SU(2) gauge coupling from Fermi constant and W mass
    const double g0sq = 4.0 * sqrt(2.0) * Gf * MW * MW;
    
    EXPECT_NEAR(ExperimentalInput::MH, MH, MH * relativeError);
    EXPECT_NEAR(ExperimentalInput::Mt, Mt, Mt * relativeError);
    EXPECT_NEAR(ExperimentalInput::MW, MW, MW * relativeError);
    EXPECT_NEAR(ExperimentalInput::MZ, MZ, MZ * relativeError);

    EXPECT_NEAR(ExperimentalInput::alphaQCD, alphaQCD, alphaQCD * relativeError);
    EXPECT_NEAR(ExperimentalInput::gs2, gs2, gs2 * relativeError);
    EXPECT_NEAR(ExperimentalInput::Gf, Gf, Gf * relativeError);
    EXPECT_NEAR(ExperimentalInput::g0sq, g0sq, g0sq * relativeError);
}



TEST(RenormalizationTest, TreeLevelMatching) {

    // Calculate MS-bar parameters at Z pole, tree level matching

    double inputScale = ExperimentalInput::MZ;

    // This has our test input parameters built in
    TestParams testParams;
    
    Renormalization renorm(testParams.inputParams, inputScale);

    ParameterMap MSParams = renorm.CalcMS(ELoopOrderMS::tree);

    // Expected params. These are from the old Mathematica implementation 
    // apart from a1, lambda which are from C++ (tested version)
    double ytsq = 0.984625240053557;
    double g1sq = 0.122353601948581;
    double g2sq = 0.426284721044574;
    double g3sq = 1.484088369555819;
    double msqPhi = -10304.324928;
    double lambda = 0.16997072342780253;
    double b1 = -5.05038380332448e+6;
    double b2 = 233793.127269626;
    double b3 = -75.0;
    double b4 = 0.5;
    double a1 = 333.22605590844131;
    double a2 = 4.0;

    EXPECT_DOUBLE_EQ(MSParams["ytsq"], ytsq);
    EXPECT_DOUBLE_EQ(MSParams["g1sq"], g1sq);
    EXPECT_DOUBLE_EQ(MSParams["g2sq"], g2sq);
    EXPECT_DOUBLE_EQ(MSParams["g3sq"], g3sq);
    EXPECT_DOUBLE_EQ(MSParams["msqPhi"], msqPhi);
    EXPECT_DOUBLE_EQ(MSParams["lambda"], lambda);
    EXPECT_DOUBLE_EQ(MSParams["b1"], b1);
    EXPECT_DOUBLE_EQ(MSParams["b2"], b2);
    EXPECT_DOUBLE_EQ(MSParams["b3"], b3);
    EXPECT_DOUBLE_EQ(MSParams["b4"], b4);
    EXPECT_DOUBLE_EQ(MSParams["a1"], a1);
    EXPECT_DOUBLE_EQ(MSParams["a2"], a2);
}



TEST(RenormalizationTest, OneLoopMatching) {

    // Calculate MS-bar parameters at Z pole, 1-loop matching

    double inputScale = ExperimentalInput::MZ;

    // This has our test input parameters built in
    TestParams testParams;
    
    Renormalization renorm(testParams.inputParams, inputScale);

    ParameterMap MSParams = renorm.CalcMS(ELoopOrderMS::loop1);

    // Expected params. These are from tested version of this code
    double ytsq = 0.95666898910014875;
    double g1sq = 0.128238967691677;
    double g2sq = 0.4239175257087634;
    double g3sq = 1.484088369555819;
    double msqPhi = -23138.27795660917;
    double lambda = 0.11036933569847558;
    double b1 = -1.12169010951581e+7;
    double b2 = 207370.91674140835;
    double b3 = -75.0;
    double b4 = 0.5;
    double a1 = 208.21665733826168;
    double a2 = 4.0;

    EXPECT_DOUBLE_EQ(MSParams["ytsq"], ytsq);
    EXPECT_DOUBLE_EQ(MSParams["g1sq"], g1sq);
    EXPECT_DOUBLE_EQ(MSParams["g2sq"], g2sq);
    EXPECT_DOUBLE_EQ(MSParams["g3sq"], g3sq);
    EXPECT_DOUBLE_EQ(MSParams["msqPhi"], msqPhi);
    EXPECT_DOUBLE_EQ(MSParams["lambda"], lambda);
    EXPECT_DOUBLE_EQ(MSParams["b1"], b1);
    EXPECT_DOUBLE_EQ(MSParams["b2"], b2);
    EXPECT_DOUBLE_EQ(MSParams["b3"], b3);
    EXPECT_DOUBLE_EQ(MSParams["b4"], b4);
    EXPECT_DOUBLE_EQ(MSParams["a1"], a1);
    EXPECT_DOUBLE_EQ(MSParams["a2"], a2);
}

