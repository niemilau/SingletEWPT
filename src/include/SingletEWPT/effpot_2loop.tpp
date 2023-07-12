#include "effpot.h"
#include "integrals3D.h"

/* This file calculates diagrams required for the 2-loop potential. 
Most expressions have been copy-pasted from Mathematica using CForm */

template <typename Complex>
Complex EffPot<Complex>::V2() {

	Complex res{0,0};

	// Figure-8 diagrams
	res += NeutralScalar8();
	res += ChargedScalar8();
	res += MixedScalar8();
	res += PureGauge8();
	res += MixedGaugeScalar8();

	// Sunset diagrams
	res += NeutralScalarSunset();
	res += MixedScalarSunset();
	res += VVS();
	//res += VVStest();
	res += VSS();
	res += VVV();
	res += VGG();

	// In my sign conventions, potential is minus the sum of 2-loop diagrams
	res = -res;

	return res;
}


template <typename Complex>
Complex EffPot<Complex>::NeutralScalar8() {

	// Alias for the loop integral of correct type. Cannot use 'using' because it is a static member function
	const auto I3 = Integrals::I3;

	return ( (((-a2 - 3*(b4 + lambda) + 3*(-a2 + b4 + lambda)*cos(4*theta))*I3(mh1sq)*I3(mh2sq))/4. + 
		I3(mGsq)*I3(mh1sq)*(-2*lambda*(cos(theta)*cos(theta)) - a2*(sin(theta)*sin(theta))) + 
		I3(mGsq)*I3(mh2sq)*(-(a2*(cos(theta)*cos(theta))) - 2*lambda*(sin(theta)*sin(theta))))/
		4. + (-6*lambda*(I3(mGsq)*I3(mGsq)) - 
		6.0*(I3(mh1sq)*I3(mh1sq))*(lambda*(cos(theta)*cos(theta)*cos(theta)*cos(theta)) + 
		a2*(cos(theta)*cos(theta))*(sin(theta)*sin(theta)) + 
		b4*(sin(theta)*sin(theta)*sin(theta)*sin(theta))) - 
		6.0*(I3(mh2sq)*I3(mh2sq))*(b4*(cos(theta)*cos(theta)*cos(theta)*cos(theta)) + 
		a2*(cos(theta)*cos(theta))*(sin(theta)*sin(theta)) + 
		lambda*(sin(theta)*sin(theta)*sin(theta)*sin(theta))))/8.
	);
}

template <typename Complex>
Complex EffPot<Complex>::ChargedScalar8() {

	const auto I3 = Integrals::I3;
	return -2. * lambda * I3(mGpmsq)*I3(mGpmsq);
}

template <typename Complex>
Complex EffPot<Complex>::MixedScalar8() {

	const auto I3 = Integrals::I3;
	return ( (-2*lambda*I3(mGpmsq)*I3(mGsq) + I3(mGpmsq)*I3(mh1sq)*
      (-2*lambda*(cos(theta)*cos(theta)) - a2*(sin(theta)*sin(theta))) + 
     I3(mGpmsq)*I3(mh2sq)*(-(a2*(cos(theta)*cos(theta))) - 2*lambda*(sin(theta)*sin(theta)))) / 2.
	);
}


template <typename Complex>
Complex EffPot<Complex>::PureGauge8() {

	const auto DVV = Integrals::DVV;
	Complex res = -0.5*(g2sq*DVV(mWsq,mWsq)) - (g2sq*g2sq*DVV(mWsq,mZsq))/(g1sq + g2sq);

	// Sensibility check: these diagrams should go to 0 in the symmetric phase
	if (abs(v) < vWarningThreshold && abs(res) > 100.) {
		warnings++;
		std::cout << "!!! Warning in VV diagrams. At v = " << v << ", VV = " << res << "\n";
	}

	return res;
}

template <typename Complex>
Complex EffPot<Complex>::MixedGaugeScalar8() {

	Float d = 3;
	const auto I3 = Integrals::I3;

	return ( -0.5*((-1 + d)*g2sq*I3(mGpmsq)*I3(mWsq)) - 
			((-1 + d)*((-g1sq + g2sq)*(-g1sq + g2sq))*I3(mGpmsq)*I3(mZsq))/(4.*(g1sq + g2sq)) + 
			((-1 + d)*(-0.5*(g2sq*I3(mGsq)*I3(mWsq)) - 
			(g2sq*(cos(theta)*cos(theta))*I3(mh1sq)*I3(mWsq))/2. - 
			(g2sq*I3(mh2sq)*I3(mWsq)*(sin(theta)*sin(theta)))/2.))/2. + 
			((-1 + d)*(((-g1sq - g2sq)*I3(mGsq)*I3(mZsq))/2. - 
			((g1sq + g2sq)*(cos(theta)*cos(theta))*I3(mh1sq)*I3(mZsq))/2. - 
			((g1sq + g2sq)*I3(mh2sq)*I3(mZsq)*(sin(theta)*sin(theta)))/2.))/4.
	);
}

template <typename Complex>
Complex EffPot<Complex>::NeutralScalarSunset() {

	const auto S3 = Integrals::S3;

	return ( (S3(mh2sq,mh2sq,mh2sq,scale)*((-2*(b3 + 3*b4*x)*(cos(theta)*cos(theta)*cos(theta)) + 
		3*a2*v*(cos(theta)*cos(theta))*sin(theta) - 
		(3*(a1 + 2*a2*x)*cos(theta)*(sin(theta)*sin(theta)))/2. + 
		6*lambda*v*(sin(theta)*sin(theta)*sin(theta)))*
		(-2*(b3 + 3*b4*x)*(cos(theta)*cos(theta)*cos(theta)) + 
		3*a2*v*(cos(theta)*cos(theta))*sin(theta) - 
		(3*(a1 + 2*a2*x)*cos(theta)*(sin(theta)*sin(theta)))/2. + 
		6*lambda*v*(sin(theta)*sin(theta)*sin(theta)))) + 
		S3(mh1sq,mh1sq,mh1sq,scale)*((-6*lambda*v*(cos(theta)*cos(theta)*cos(theta)) - 
		(3*(a1 + 2*a2*x)*(cos(theta)*cos(theta))*sin(theta))/2. - 
		3*a2*v*cos(theta)*(sin(theta)*sin(theta)) - 
		2*(b3 + 3*b4*x)*(sin(theta)*sin(theta)*sin(theta)))*
		(-6*lambda*v*(cos(theta)*cos(theta)*cos(theta)) - 
		(3*(a1 + 2*a2*x)*(cos(theta)*cos(theta))*sin(theta))/2. - 
		3*a2*v*cos(theta)*(sin(theta)*sin(theta)) - 
		2*(b3 + 3*b4*x)*(sin(theta)*sin(theta)*sin(theta)))))/12. + 
		(S3(mh2sq,mGsq,mGsq,scale)*((-0.5*((a1 + 2*a2*x)*cos(theta)) + 2*lambda*v*sin(theta))*
		(-0.5*((a1 + 2*a2*x)*cos(theta)) + 2*lambda*v*sin(theta))) + 
		S3(mh1sq,mGsq,mGsq,scale)*((-2*lambda*v*cos(theta) - ((a1 + 2*a2*x)*sin(theta))/2.)*
		(-2*lambda*v*cos(theta) - ((a1 + 2*a2*x)*sin(theta))/2.)) + 
		S3(mh1sq,mh1sq,mh2sq,scale)*((-0.5*((a1 + 2*a2*x)*(cos(theta)*cos(theta)*cos(theta))) - 
		2*(a2 - 3*lambda)*v*(cos(theta)*cos(theta))*sin(theta) + 
		(a1 - 2*b3 + 2*a2*x - 6*b4*x)*cos(theta)*(sin(theta)*sin(theta)) + 
		a2*v*(sin(theta)*sin(theta)*sin(theta)))*
		(-0.5*((a1 + 2*a2*x)*(cos(theta)*cos(theta)*cos(theta))) - 
		2*(a2 - 3*lambda)*v*(cos(theta)*cos(theta))*sin(theta) + 
		(a1 - 2*b3 + 2*a2*x - 6*b4*x)*cos(theta)*(sin(theta)*sin(theta)) + 
		a2*v*(sin(theta)*sin(theta)*sin(theta)))) + 
		S3(mh1sq,mh2sq,mh2sq,scale)*((-(a2*v*(cos(theta)*cos(theta)*cos(theta))) + 
		(a1 - 2*b3 + 2*a2*x - 6*b4*x)*(cos(theta)*cos(theta))*sin(theta) - 
		6*lambda*v*cos(theta)*(sin(theta)*sin(theta)) - 
		((a1 + 2*a2*x)*(sin(theta)*sin(theta)*sin(theta)))/2. + 
		a2*v*sin(theta)*sin(2*theta))*
		(-(a2*v*(cos(theta)*cos(theta)*cos(theta))) + 
		(a1 - 2*b3 + 2*a2*x - 6*b4*x)*(cos(theta)*cos(theta))*sin(theta) - 
		6*lambda*v*cos(theta)*(sin(theta)*sin(theta)) - 
		((a1 + 2*a2*x)*(sin(theta)*sin(theta)*sin(theta)))/2. + 
		a2*v*sin(theta)*sin(2*theta))))/4. 
	);
}


template <typename Complex>
Complex EffPot<Complex>::MixedScalarSunset() {

	const auto S3 = Integrals::S3;

	return ( (S3(mh2sq,mGpmsq,mGpmsq,scale)*((-0.5*((a1 + 2*a2*x)*cos(theta)) + 2*lambda*v*sin(theta))*
		(-0.5*((a1 + 2*a2*x)*cos(theta)) + 2*lambda*v*sin(theta))) + 
		S3(mh1sq,mGpmsq,mGpmsq,scale)*((-2*lambda*v*cos(theta) - ((a1 + 2*a2*x)*sin(theta))/2.)*
		(-2*lambda*v*cos(theta) - ((a1 + 2*a2*x)*sin(theta))/2.)))/2.
	);
}

template <typename Complex>
Complex EffPot<Complex>::VVS() {

	/* VVS diagrams are numerically challenging because the general expression for DVVS integral
	is IR divergent for gauge mass -> 0. In the mass = 0 case we could use the pre-calculated special case,
	but for small nonzero mass it's still necessary to use the general form. The combination couplings * integral 
	is still IR safe but may need large numerical precision for small gauge masses.
	Fundamentally the issue is that we are adding small m^2 to another small m^2 both inside the integrals and 
	then summing up small diagrams, and floating-point errors accumulate here. */ 
	/* Instead of finding (cumbersome) workarounds or using slow multiprecision,
	I note that since VVS -> 0 as v -> 0, it is a good approximation to simply drop these diagrams if v is small enough. 
	I've tested (with 32-digit precision) that already when v = 1e-4, these diagrams add up to 1e-13 in GeV^3 units. 
	This was in our BM3 point, and the test should not be sensitive to details of the scalar potential because the diagrams here
	are proportional to gauge couplings. */

	Complex res{0,0};

	bool bSmallFieldApproximation = true;
	if (bSmallFieldApproximation && abs(v) < vWarningThreshold) {
		// Just return 0, the error should be negligible for small v
		return res; 
	}  

	const auto DVVS = Integrals::DVVS;

	Float v2 = v*v;
	
	res = ( (g1sq*(g2sq*g2sq)*(v2)*DVVS(mGpmsq,mWsq,0,scale))/(4.*(g1sq + g2sq)) + 
		(g1sq*g1sq*g2sq*(v2)*DVVS(mGpmsq,mWsq,mZsq,scale))/(4.*(g1sq + g2sq)) + 
		((g2sq*g2sq*(v2)*(cos(theta)*cos(theta))*DVVS(mh1sq,mWsq,mWsq,scale))/4. + 
		(g2sq*g2sq*(v2)*DVVS(mh2sq,mWsq,mWsq,scale)*(sin(theta)*sin(theta)))/4.)/2. + 
		(((g1sq + g2sq)*(g1sq + g2sq)*(v2)*(cos(theta)*cos(theta))*DVVS(mh1sq,mZsq,mZsq,scale))/
		4. + ((g1sq + g2sq)*(g1sq + g2sq)*(v2)*DVVS(mh2sq,mZsq,mZsq,scale)*
		(sin(theta)*sin(theta)))/4.)/4.
	);

	// Sensibility check: these diagrams should go to 0 in the symmetric phase
	if (abs(v) < vWarningThreshold && abs(res) > 100.) {
		warnings++;
		std::cout << "!!! Warning in VVS diagrams. At v = " << v << ", VVS = " << res << "\n";
	}

	return res;
}


template <typename Complex>
Complex EffPot<Complex>::VSS() {

	const auto DVSS = Integrals::DVSS;

	return ( ((g1sq*g2sq*DVSS(mGpmsq,mGpmsq,0,scale))/(g1sq + g2sq) + 
		((g1sq - g2sq)*(g1sq - g2sq)*DVSS(mGpmsq,mGpmsq,mZsq,scale))/(4.*(g1sq + g2sq)))/2. + 
		(g2sq*DVSS(mGsq,mGpmsq,mWsq,scale))/4. + 
		(g2sq*(cos(theta)*cos(theta))*DVSS(mh1sq,mGpmsq,mWsq,scale))/4. + 
		(g2sq*DVSS(mh2sq,mGpmsq,mWsq,scale)*(sin(theta)*sin(theta)))/4. + 
		(((g1sq + g2sq)*(cos(theta)*cos(theta))*DVSS(mh1sq,mGsq,mZsq,scale))/4. + 
		((g1sq + g2sq)*DVSS(mh2sq,mGsq,mZsq,scale)*(sin(theta)*sin(theta)))/4.)/2.
	);
}

template <typename Complex>
Complex EffPot<Complex>::VVV() {

	const auto DVVV = Integrals::DVVV;
	Complex res = ( ((g1sq*g2sq*DVVV(mWsq,mWsq,0,scale))/(g1sq + g2sq) + 
    	(g2sq*g2sq*DVVV(mWsq,mWsq,mZsq,scale))/(g1sq + g2sq))/2.
	);

	// Sensibility check: these diagrams should go to 0 in the symmetric phase
	if (abs(v) < vWarningThreshold && abs(res) > 100.) {
		warnings++;
		std::cout << "!!! Warning in VVV diagrams. At v = " << v << ", VVV = " << res << "\n";
	}
	return res;
}

template <typename Complex>
Complex EffPot<Complex>::VGG() {

	const auto DVGG = Integrals::DVGG;
	Complex res = ( (-2*g1sq*g2sq*DVGG(mWsq,scale))/(g1sq + g2sq) - 
		(2*(g2sq*g2sq)*DVGG(mWsq,scale))/(g1sq + g2sq) - (g2sq*g2sq*DVGG(mZsq,scale))/(g1sq + g2sq)
	);

	// Sensibility check: these diagrams should go to 0 in the symmetric phase
	if (abs(v) < vWarningThreshold && abs(res) > 100.) {
		warnings++;
		std::cout << "!!! Warning in VGG diagrams. At v = " << v << ", VGG = " << res << "\n";
	}
	return res;
}


template <typename Complex>
void EffPot<Complex>::TestDiagrams(typename EffPot<Complex>::Float v, typename EffPot<Complex>::Float x) {

	SetBackgroundFields(v, x);

	using std::cout;
	cout.precision(12);

	cout << "=== Start diagram test ===\n";
	cout << "NeutralScalar8(): " << NeutralScalar8() << "\n";
	cout << "ChargedScalar8();: " << ChargedScalar8() << "\n";
	cout << "MixedScalar8(): " << MixedScalar8() << "\n";
	cout << "PureGauge8(): " << PureGauge8() << "\n";
	cout << "MixedGaugeScalar8(): " << MixedGaugeScalar8() << "\n";
	cout << "NeutralScalarSunset(): " << NeutralScalarSunset() << "\n";
	cout << "MixedScalarSunset(): " << MixedScalarSunset() << "\n";
	cout << "VVS(): " << VVS() << "\n";
	cout << "VSS(): " << VSS() << "\n";
	cout << "VVV(): " << VVV() << "\n";
	cout << "VGG(): " << VGG() << "\n";
	cout << "=== End diagram test ===\n";
}