#ifndef TESTPARAMS_H
#define TESTPARAMS_H

#include "common.h"

struct TestParams {

    // Input parameters for SM + singlet
    ParameterMap inputParams;

    TestParams() {
        inputParams["Mh1"] = 125.10;
        inputParams["Mh2"] = 600.0;
        inputParams["a2"] = 4.0;
        inputParams["sinTheta"] = 0.12;
        inputParams["b3"] = -75.0;
        inputParams["b4"] = 0.5;
        
    }

    double T = 143.2;
};



#endif