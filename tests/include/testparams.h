#ifndef TESTPARAMS_H
#define TESTPARAMS_H

#include "common.h"

struct TestParams {

    double relativeError = 1e-15;

    // Input parameters for SM + singlet
    ParameterMap inputParams;

public:
    TestParams() {
        inputParams["mh2"] = 600.0;
    }

    double mh2 = 600.0;
    double a2 = 4.0;
    double sinTheta = 0.12;
    double b3 = -75.0;
    double b4 = 0.5;

    double T = 143.2;
};



#endif