#ifndef DIMRED_H
#define DIMRED_H

#include "common.h"


namespace DimRed {

    /* Constant cc = 1/2 * (log(8pi/9) + Zeta'(2)/Zeta(2) - 2EulerGamma) */
    constexpr double CC = -0.34872273635229090400412764915100;

    /* Perform dimensional reduction: 4D SM + singlet --> 3D SM + singlet + A0. 
    Parameter conventions are as in 2103.07467.
    loopOrder specifies accuracy: LO = 1-loop thermal masses only
                                 NLO = full O(g^4) matching including 2-loop thermal masses
                          NLONo2Loop = Same as NLO but do not include 2-loop effects.

    If bDoDim6 = true, will compute also dimension 5 and 6 coefficients.
    3D running to scale finalScale is performed here as well. */
    ParameterMap IntegrateHardModes(const double T, const ParameterMap &MSParams, 
        const ELoopOrderDR loopOrderDR, double const finalScale, const bool bDoDim6, const bool bNLOCubics);

    /* Integrate out A0, B0, C0, ie reduce 3D SM + singlet + A0 -> 3D SM + singlet. 
    No RG running is done here, so end scale is the same as the initial scale. 
    Dim-5 and 6 operators are TODO */
    ParameterMap IntegrateSoftModes(const ParameterMap &MSParams, 
		const ELoopOrderDR loopOrderDR, const bool bDoDim6, const bool bNLOCubics); 

    // Integrate out hard and soft modes, ie. do 4D SM + singlet -> 3D SM + singlet.
    ParameterMap DoFullDimRed(const double T, const ParameterMap &MSParams, 
        const ELoopOrderDR loopOrderDR, double const finalScale, const bool bDoDim6, const bool bNLOCubics);
};

#endif