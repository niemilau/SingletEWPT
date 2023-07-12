#include <gtest/gtest.h>

#include "testparams.h"
#include "renormalization.h"

/* This file checks that the T=0 matching between inputs and renormalized parameters is OK */

TEST(RenormalizationTest, TreeLevelMatching) {

    // Calculate MS-bar parameters at Z pole
    double inputScale = ExperimentalInput::MZ;
    //Renormalization renorm(scanner.currentInput, inputScale);

    //ParameterMap MSParams = renorm.CalcMS(scanner.loopOrderMS);

}

