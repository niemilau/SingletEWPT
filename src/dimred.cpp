#include "dimred.h"


ParameterMap DimRed::IntegrateHardModes(const double T, const ParameterMap &MSParams, 
		const ELoopOrderDR loopOrderDR, double const finalScale, const bool bDoDim6, const bool bNLOCubics) 		
{

	// TODO (if ever needed?)
	if (bNLOCubics) {
		Die("!!! NNLO cubic terms in DR not implemented yet...", 666);
	}

	// Logarithms etc
	double scale = GetFromMap(MSParams, "RGScale");
	double Lb = 2.0*log(scale/T) - 2.0*(log(4.0*PI) - EULERGAMMA);
	double Lf = Lb + 4.0*log(2.0);
	
	double cc = DimRed::CC;
	// 2loop logarithm
	double ccPlusLog = cc + log(3*T/finalScale);

	// 1/(4pi)^2
	double pi4sqInv = 1.0 / (16.0 * PI*PI);

	// How many Higgs doublets
	int Nd = 1;
	// How many fermion _generations_
	int Nf = 3;

	// 4D parameters
	double gs2 = GetFromMap(MSParams, "g3sq");
	double g2sq = GetFromMap(MSParams, "g2sq");
	double g1sq = GetFromMap(MSParams, "g1sq");
	double ytsq = GetFromMap(MSParams, "ytsq");
	double msqPhi = GetFromMap(MSParams, "msqPhi");
	double lambda = GetFromMap(MSParams, "lambda");
	double b1 = GetFromMap(MSParams, "b1");
	double b2 = GetFromMap(MSParams, "b2");
	double b3 = GetFromMap(MSParams, "b3");
	double b4 = GetFromMap(MSParams, "b4");
	double a1 = GetFromMap(MSParams, "a1");
	double a2 = GetFromMap(MSParams, "a2");

	// Begin DR matching

	ParameterMap params3D;
	params3D["RGScale"] = finalScale;

	double T2 = T*T;
	bool bIsNLO = (loopOrderDR != ELoopOrderDR::LO);

	// U(1) Debye mass squared
	double mD1sq = g1sq * T2 * (Nd / 6.0 + 5.0*Nf / 9.0);
	// SU(2) Debye mass squared
	double mD2sq = g2sq * T2 * ((4.0+Nd) / 6.0 + Nf / 3.0);
	// SU(3) Debye mass squared
	double mD3sq = gs2 * T2 * (1.0 + Nf / 3.0);

	if (bIsNLO && loopOrderDR != ELoopOrderDR::NLONo2Loop) {
		// TODO: Debye masses at 2-loops
	}
	params3D["mD1sq"] = mD1sq;
	params3D["mD2sq"] = mD2sq;
	params3D["mD3sq"] = mD3sq;
	 
	// U(1) gauge coupling squared
	
		double g1sq3D = g1sq * T;
		if (bIsNLO) {
			g1sq3D *= ( 1 + pi4sqInv * g1sq * (-Lb * Nd / 6.0 - Lf * 20.0 * Nf / 9.0) );
		}
		params3D["g1sq"] = g1sq3D;
	

	// SU(2) gauge coupling squared
	
		double g2sq3D = g2sq * T;
		if (bIsNLO) {
			g2sq3D *= ( 1 + pi4sqInv * g2sq * (Lb * (44.0 - Nd) / 6.0  + 2.0/3.0 - Lf * 4.0*Nf / 3.0));
		}
		params3D["g2sq"] = g2sq3D;
	

	// phi^2 A0^2 coupling
	
		double h3D = g2sq * T / 4.0;
		if (bIsNLO) {
			h3D *= ( 1 + pi4sqInv * ( ((44.0-Nd)/6.0 * Lb + 53/6.0 - Nd/3.0 - 4*Nf/3.0 *(Lf-1.0) ) * g2sq 
					+ g1sq / 2.0 - 6*ytsq + 12*lambda ) 
				); 
		}   
		params3D["h3"] = h3D;
	

	// phi^2 B0^2 coupling
	double hp3D =  g1sq * T / 4.0;
	if (bIsNLO) {
		hp3D *= ( 1 + pi4sqInv * ( 3*g2sq / 2.0 + (0.5 - Nd/6.0 * (2+Lb) - 20*Nf / 9.0 * (Lf-1.0)) * g1sq 
				- 34.0/3.0*ytsq + 12*lambda )
			); 
	}   
	params3D["hp3"] = hp3D;

	// phi^2 A0 B0 coupling
	
	double hpp3D = sqrt(g2sq * g1sq) * T / 2.0;
	if (bIsNLO) {
		hpp3D *= ( 1 + pi4sqInv * ( -(5+Nd)/6.0 * g2sq + (3-Nd)/6.0 * g1sq 
				+ Lb * ( (44-Nd)/12.0 *g2sq - Nd/12.0 * g1sq ) 
				- Nf * (Lf-1.0) * (2.0/3.0 * g2sq + 10.0/9.0 *g1sq) + 2*ytsq + 4*lambda)
			); 
	}   
	params3D["hpp3"] = hpp3D;
	

	// phi^2 C0 coupling
	double omega3D = 0;
	if (bIsNLO) {
		omega3D += -T * pi4sqInv * 2.0 * gs2 * ytsq;
	}   
	params3D["omega3"] = omega3D;
	// A0, B0, C0 self interactions not needed, skip...

	// Singlet - A0, B0 couplings
	double x3D = 0, xp3D = 0, y3D = 0, yp3D = 0;
	if (bIsNLO) {
		x3D = sqrt(T) * pi4sqInv * g2sq * a1;
		xp3D = sqrt(T) * pi4sqInv * g1sq * a1;
		y3D = T * pi4sqInv * g2sq * a2 * 0.5;
		yp3D = T * pi4sqInv * g1sq * a2 * 0.5;
	}
	params3D["x3"] = x3D;
	params3D["xp3"] = xp3D;
	params3D["y3"] = y3D;
	params3D["yp3"] = yp3D;

	// Higgs quartic self interaction
	double lambda3D = lambda*T;
	if (bIsNLO) {
		lambda3D += T * pi4sqInv * ((2-3*Lb) / 16.0 * (3*g2sq*g2sq + 2*g2sq*g1sq + g1sq*g1sq) 
					+ 3*ytsq*Lf * (ytsq - 2*lambda) 
					+ Lb*(3.0/2.0 * (3*g2sq + g1sq)*lambda - 12*lambda*lambda - 1.0/4.0 * a2*a2)
		);
	}
	params3D["lambda"] = lambda3D;

	// Singlet cubic self interaction
	double b33D = b3 * sqrt(T);
	if (bIsNLO) {
		b33D += sqrt(T) * pi4sqInv * (-3 * Lb * (3 * b4 * b3 + 1.0/2.0 * a1 * a2)); 
	}
	params3D["b3"] = b33D;

	// Singlet quartic self interaction
	double b43D = b4 * T;
	if (bIsNLO) {
		b43D += -T * pi4sqInv * Lb * (a2*a2 + 9*b4*b4); 
	}
	params3D["b4"] = b43D;

	// Singlet-Higgs cubic coupling
	double a13D = a1 * sqrt(T);
	if (bIsNLO) {
		a13D += sqrt(T) * pi4sqInv * (-3*Lf*ytsq*a1 + 
				Lb * (-2*a2*b3 + a1 * (3.0/4.0 * (3*g2sq + g1sq) - 6*lambda - 2*a2)) 
			); 
	}
	params3D["a1"] = a13D;

	// Singlet-Higgs quartic coupling
	
		double a23D = a2 * T;
		if (bIsNLO) {
			a23D += T * pi4sqInv * (-3 * Lf * ytsq * a2 
					+ Lb * a2 * (3.0/4.0 *(3*g2sq + g1sq) - 6*lambda - 2*a2 - 3*b4)
				);
		}
		params3D["a2"] = a23D;
	

	/* 2 loop matchings. */

	// Higgs mass parameter squared
	{
		double msqPhi3D = msqPhi + T2 / 16.0 * (3*g2sq + g1sq + 4*ytsq + 8*lambda) + T2 * a2/24.0;
		if (bIsNLO && loopOrderDR != ELoopOrderDR::NLONo2Loop) {

			// 2-loop part taken from Mathematica using CForm 
			msqPhi3D += pi4sqInv * (msqPhi*(((3*(g1sq + 3*g2sq))/4. - 6*lambda)*Lb - 3*Lf*ytsq) + 
						T*T*((g1sq*g1sq)/288. - (3*g1sq*g2sq)/16. + (167*(g2sq*g2sq))/96. + ((g1sq + 3*g2sq)*lambda)/4. + 
						((-5*(g1sq*g1sq))/48. - (3*g1sq*g2sq)/16. + (17*(g2sq*g2sq))/16. + (3*(g1sq + 3*g2sq)*lambda)/4. - 
						6*(lambda*lambda))*Lb + ((5*(g1sq*g1sq))/108. + (g2sq*g2sq)/12.)*Nf - 
						((11*g1sq)/48. + (3*g2sq)/16. + 2*gs2)*ytsq + 
						Lf*(-(((5*(g1sq*g1sq))/36. + (g2sq*g2sq)/4.)*Nf) + ((17*g1sq)/48. + (9*g2sq)/16. + 2*gs2 - 3*lambda)*ytsq + 
						(3*(ytsq*ytsq))/8.) + (((5*(g1sq*g1sq))/6. + (3*(g2sq*g2sq))/2.)*Nf + 
						((-47*g1sq)/72. - (21*g2sq)/8. + (8*gs2)/3. + 9*lambda)*ytsq - (3*(ytsq*ytsq))/2.)*log(2.0)) + 
						((-5*(g1sq3D*g1sq3D))/16. - (9*g1sq3D*g2sq3D)/8. + (39*(g2sq3D*g2sq3D))/16. + 12*g2sq3D*h3D - 6*(h3D*h3D) - 
						2*(hp3D*hp3D) - 3*(hpp3D*hpp3D) + 3*(g1sq3D + 3*g2sq3D)*lambda3D - 12*(lambda3D*lambda3D))*ccPlusLog 
	 				); 
			// Add singlet 2-loop contributions
			msqPhi3D += pi4sqInv * ( -0.25*(a1*a1*Lb) - (a2*Lb*b2)/2. - (a2*((10*a2)/3. + 2*b4 + 4*lambda)*Lb*(T*T))/16. + 
						(a2*(T*T)*((3*(g1sq + 3*g2sq)*Lb)/4. - 3*Lf*ytsq))/24. - (a23D*a23D* ccPlusLog )/2. 
					);	
		}
		params3D["msqPhi"] = msqPhi3D;
	}

	// Singlet mass parameter squared (i.e. b2)
	{
		double b23D = b2 + T2 * (a2/6.0 +  b4/4.0);
		if (bIsNLO && loopOrderDR != ELoopOrderDR::NLONo2Loop) {

			b23D += -pi4sqInv*Lb * (2*b3*b3 + a1*a1/2.0 + 2*a2*msqPhi + 3*b4*b2)
					+ pi4sqInv * T2 * ( (2+3*Lb)/24.0 * (3*g2sq + g1sq)*a2 
					- Lb * ( a2 * (lambda + 7.0/12.0*a2 + b4/2.0) + 9.0/4.0 * b4*b4)
					- 1.0/4.0 * a2*ytsq * (3*Lb - Lf))
					+ pi4sqInv*ccPlusLog * ( (3*g2sq3D + g1sq3D)*a23D - 2*a23D*a23D - 6*b43D*b43D);
		}
		params3D["b2"] = b23D;
	}

	// Singlet linear coupling
	{
		double b13D = (1.0 / sqrt(T)) * (b1 + T2/12.0 * (b3 + a1));
		if (bIsNLO && loopOrderDR != ELoopOrderDR::NLONo2Loop) {

			b13D += -pi4sqInv * Lb / sqrt(T) * (a1 * msqPhi + b3 * b2)
					+ pi4sqInv * T*sqrt(T) * ( (2+3*Lb) / 48.0 * (3*g2sq + g1sq)*a1 - Lb/2.0 * ( (lambda + 7.0/12.0*a2)*a1 
					+ (1.0/3.0 *a2 + 3.0/2.0*b4)*b3 ) - 1.0/8.0*a1*ytsq * (3*Lb-Lf) )
					- pi4sqInv * ccPlusLog*( 2*b33D*b43D - 0.5*a13D*(3*g2sq3D + g1sq3D -2*a23D) );
		}
		params3D["b1"] = b13D;
	}

	/* Dimension 5 and 6 operators. Terms with gauge couplings are NOT included: 
	1) they are small, 2) they came with some gauge dependence that we did not understand */
	double c05 = 0.0, c23 = 0.0, c41 = 0.0, c60 = 0.0, c06 = 0.0, c42 = 0.0, c24 = 0.0;
	if (bDoDim6) {
		double zeta3 = 1.2020569031595942854;
		double zeta4pi4 = zeta3 / (4*4*4*4 * PI*PI*PI*PI);
		double invSqrtT = 1.0 / sqrt(T);

		c05 = zeta4pi4*invSqrtT * (18*b3*b4*b4 + 0.5 * a1*a2*a2);
		c23 = zeta4pi4*invSqrtT * (a2/4.0) * ( 8*(a2 + 3*b4)*(2*b3+a1) + 24*lambda*a1);
		c41 = zeta4pi4*invSqrtT * ( 12*a1*a2*lambda + 2*a2*a2*(b3 + a1) + 24*lambda*lambda*a1 ); 
		c60 = zeta4pi4 * (80*lambda*lambda*lambda - 28 *ytsq*ytsq*ytsq + 1.0/3.0 * a2*a2*a2);
		c06 = zeta4pi4 * (1.0/6.0 * a2*a2*a2 + 9*b4*b4*b4);
		c42 = zeta4pi4 * a2 * (24*lambda*lambda + 2*a2*a2 + a2*(12*lambda + 3*b4));
		c24 = zeta4pi4 * (a2*a2*a2 + 9*a2*b4*b4 + 3*a2*a2*lambda + 6*a2*a2*b4);
	}

	params3D["c05"] = c05;
	params3D["c23"] = c23;
	params3D["c41"] = c41;
	params3D["c60"] = c60;
	params3D["c06"] = c06;
	params3D["c42"] = c42;
	params3D["c24"] = c24;


	return params3D;
}



ParameterMap DimRed::IntegrateSoftModes(const ParameterMap &params3D, 
		const ELoopOrderDR loopOrderDR, const bool bDoDim6, const bool bNLOCubics) 
{

	// TODO (if ever needed?)
	if (bNLOCubics) {
		Die("!!! NNLO cubic terms in DR not implemented yet...", 667);
	}

	bool bIsNLO = (loopOrderDR != ELoopOrderDR::LO);

	// Matching scale
	double scale = GetFromMap(params3D, "RGScale");

	double mD1 = sqrt(GetFromMap(params3D, "mD1sq"));
	double mD2 = sqrt(GetFromMap(params3D, "mD2sq"));
	double mD3 = sqrt(GetFromMap(params3D, "mD3sq"));

	double h3 = GetFromMap(params3D, "h3");
	double hp3 = GetFromMap(params3D, "hp3");
	double hpp3 = GetFromMap(params3D, "hpp3");
	double omega3 = GetFromMap(params3D, "omega3");

	// I only include these in 1-loop mass corrections, drop elsewhere.
	double x3 = GetFromMap(params3D, "x3");
	double xp3 = GetFromMap(params3D, "xp3");
	double y3 = GetFromMap(params3D, "y3");
	double yp3 = GetFromMap(params3D, "yp3");

	double g1sq = GetFromMap(params3D, "g1sq");
	double g2sq = GetFromMap(params3D, "g2sq");
	double lambda = GetFromMap(params3D, "lambda");
	double msqPhi = GetFromMap(params3D, "msqPhi");
	double b1 = GetFromMap(params3D, "b1");
	double b2 = GetFromMap(params3D, "b2");
	double b3 = GetFromMap(params3D, "b3");
	double b4 = GetFromMap(params3D, "b4");
	double a1 = GetFromMap(params3D, "a1");
	double a2 = GetFromMap(params3D, "a2");

	// Resulting parameters at the ultrasoft scale
	ParameterMap paramsUS;
	paramsUS["RGScale"] = scale; // no running done in this function

	/* ===== Matching ===== */
	
	// U(1) gauge coupling squared
	{
		double g1sqUS = g1sq;
		paramsUS["g1sq"] = g1sqUS;
	}
	// SU(2) gauge coupling squared
	{
		double g2sqUS = g2sq * (1.0 - g2sq / (6*4*PI * mD2));
		paramsUS["g2sq"] = g2sqUS;
	}
	// Higgs quartic self interaction 
	{
		double lambdaUS = lambda - 1.0/(2*4*PI) * (3*h3*h3 / mD2 + hp3*hp3 / mD1 + hpp3*hpp3 / (mD1 + mD2));
		paramsUS["lambda"] = lambdaUS;
	}
	// Singlet cubic self interaction
	{
		double b3US = b3; // plus O(x3 y3 / mD)
		paramsUS["b3"] = b3US;
	}
	// Singlet quartic self interaction
	{
		double b4US = b4; // plus O(y3^2 / mD)
		paramsUS["b4"] = b4US;
	}
	// Singlet-Higgs cubic coupling
	{
		double a1US = a1; // plus O(h3 x3 / mD)
		paramsUS["a1"] = a1US;
	}
	// Singlet-Higgs quartic coupling
	{
		double a2US = a2; // plus O(h3 y3 / mD)
		paramsUS["a2"] = a2US;
	}

	/* All 2-loop corrections to singlet b1 and b2 go beyond g^4 accuracy */
	// Singlet tadpole interaction
	{
		double b1US = b1;
		// 1-loop effects from x3, xp3 are NLO because x3 itself is O(g^4)
		if (bIsNLO) {
			b1US += -1.0/(4*PI) * (3*mD2 * x3 + mD1 * xp3);
		} 
		paramsUS["b1"] = b1US;
	}
	// Singlet mass term squared
	{
		double b2US = b2;
		// 1-loop effects from y3, yp3 are NLO because y3 itself is O(g^4)
		if (bIsNLO) {
			b2US += -1.0/(2*PI) * (3*mD2*y3 + mD1*yp3) - 1.0/(4*PI) * (3*x3*x3 / mD2 + xp3*xp3 / mD1);
		} 
		paramsUS["b2"] = b2US;
	}
	// Higgs mass term squared
	{
		double msqPhiUS = msqPhi - 1.0/(4*PI) * (3*h3*mD2 + hp3*mD1);
		
		// 1-loop effects from omega3 is NLO because omega3 itself is O(g^4)
		if (bIsNLO) {
			msqPhiUS += -1.0/(4*PI) * 8*omega3*mD3;
		} 
		// 2 loop, all singlet contributions here go beyond O(g^4)
		if (bIsNLO && loopOrderDR != ELoopOrderDR::NLONo2Loop) {
			msqPhiUS += 1.0 / (4*4*PI*PI) * (3*g2sq*h3 - 3*h3*h3 - hp3*hp3 - 3.0/2.0 * hpp3*hpp3 
						+ log(scale/(2.0*mD2)) * (-3.0/4.0 * g2sq*g2sq + 12*g2sq*h3 - 6*h3*h3)
						- 2*hp3*hp3 * log(scale/(2.0*mD1)) - 3*hpp3*hpp3 * log(scale/(mD1 + mD2))
			);
		}
		paramsUS["msqPhi"] = msqPhiUS;
	}

	/** Dimension 5 and 6 operators **/
	// TODO. for now just copy them. corrections from A0 are subleading anyway.
	if (bDoDim6) {

	} 

	paramsUS["c05"] = GetFromMap(params3D, "c05");
	paramsUS["c23"] = GetFromMap(params3D, "c23");
	paramsUS["c41"] = GetFromMap(params3D, "c41");
	paramsUS["c60"] = GetFromMap(params3D, "c60");
	paramsUS["c06"] = GetFromMap(params3D, "c06");
	paramsUS["c42"] = GetFromMap(params3D, "c42");
	paramsUS["c24"] = GetFromMap(params3D, "c24");

	return paramsUS;
}


ParameterMap DimRed::DoFullDimRed(const double T, const ParameterMap &MSParams, 
        const ELoopOrderDR loopOrderDR, double const finalScale, const bool bDoDim6, const bool bNLOCubics) 
{

	// DR step 1
	ParameterMap paramsDR = DimRed::IntegrateHardModes(T, MSParams, loopOrderDR, finalScale, bDoDim6, bNLOCubics);

	// DR step 2
	paramsDR = DimRed::IntegrateSoftModes(paramsDR, loopOrderDR, bDoDim6, bNLOCubics);

	if (GetFromMap(paramsDR, "RGScale") != finalScale) {
		// TODO evolve using the 3D RG equations
	}
	
	return paramsDR;
}