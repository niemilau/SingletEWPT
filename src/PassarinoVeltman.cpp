/** Implementation of Passarino-Veltman 1-loop integrals in d = 4-2eps dimensions **/
// The integration measure here is int_p = mu^{2eps}/i * int d^dp / (2pi)^d , with d = 4-2eps. 
// See Appendix A.3 in https://arxiv.org/pdf/2103.07467.pdf and references there for definitions.

#include "renormalization.h"


Renormalization::Integral Renormalization::A0(const Float &m) {
	if (IsZeroMass(m)) {
		return Integral(0., 0.);
	}
	Float div = pi4sqInv * m*m;
	Float fin = pi4sqInv * m*m * (1. + log(scale*scale / (m*m)));
	return Integral(fin, div);
}

Renormalization::Complex Renormalization::Fb(const Float &k, const Float &m1, const Float &m2) {
	Float k2 = k*k;
	Float re = 0.;
	Float im = 0.;

	if (k2 <= (m1 - m2)*(m1 - m2)) {
		Float r1 = sqrt( (m1+m2)*(m1+m2) - k2 );
		Float r2 = sqrt( (m1-m2)*(m1-m2) - k2 );
		re = 1./k2 * r1 * r2 * log( (r1+r2) / (r1-r2) );
	} else if ( (m1 - m2)*(m1 - m2) < k2 && k2 < (m1 + m2)*(m1 + m2) ) {

		Float r1 = sqrt( (m1+m2)*(m1+m2) - k2 );
		Float r2 = sqrt( k2 - (m1-m2)*(m1-m2) );
		re = -2. / k2 * r1*r2 * atan(r2 / r1);
	} else {

		// now k2 >= (m1+m2)^2. Fb develops an imaginary part here
		Float r1 = sqrt( k2 - (m1+m2)*(m1+m2) );
		Float r2 = sqrt( k2 - (m1-m2)*(m1-m2) );
		re = -1. / k2 * r1*r2 * log( (r1+r2) / (r2-r1) );
		im = 1. / k2 * r1*r2 * PI; 
	}

	// In principle, should have special cases here for vanishing p, m1, m2 etc. 
	// But since Fb appears only in B0, it's enough to have special cases for B0 which have clean expressions without Fb.

	return Renormalization::Complex{re, im}; 
}



Renormalization::Integral Renormalization::B0(const Float &k, const Float &M1, const Float &M2) {
	// Use symmetry M1 <-> M2 to reduce the number of special cases we need to write down 
	Float m1 = M1; 
	Float m2 = M2;
	if (IsZeroMass(m1) && !IsZeroMass(m2)) {
		// Swap masses so that m1 != 0, m2 == 0
		Float temp = m2;
		m2 = m1;
		m1 = temp; 
	}

	// Divergent part
	Float div = pi4sqInv;
	// MS-bar scale squared
	const Float musq = scale*scale;
	const Float k2 = k*k;
	Float fin = 0.; 

	// Special cases for the finite part; will drop imaginary parts
	// B0(k, 0, 0)
	if (!IsZeroMass(k) && IsZeroMass(m1) && IsZeroMass(m2)) {
		fin = pi4sqInv * (2. + log(musq / k2)); // plus -pi4sqInv * i*pi, 
	} 
	// B0(m, m, 0)
	else if (!IsZeroMass(k) && IsZeroMass(m2) && IsZeroMass(k - m1)) {
		fin = pi4sqInv * (2. + log(musq / k2));
	}
	// B0(k, m, 0)
	else if (!IsZeroMass(k) && !IsZeroMass(m1) && IsZeroMass(m2)) {
		Float msq = m1*m1;
		// Now there is a complex log if k2 < msq, take real part only
		Float logReal = log(abs(1. - k2 / msq));
		fin = pi4sqInv * (2. + log(musq / msq) + (-1. + msq/k2) * logReal );
	} 
	// B0(m, m, m)
	else if (!IsZeroMass(k) && IsZeroMass(k-m1) && IsZeroMass(k-m2)) {
		fin = pi4sqInv * (2. - PI / sqrt(3.) + log(musq / (m1*m1)));
	}
	// B0(0, 0, 0) 
	else if (IsZeroMass(k) && IsZeroMass(m1) && IsZeroMass(m2)) {
		fin = 0.;
		div = 0.;
	}
	// B0(0, m, m)
	else if (IsZeroMass(k) && !IsZeroMass(m1) && IsZeroMass(m1 - m2)) {
		fin = pi4sqInv * log(musq / (m1*m1));
	} 
	// B0(0, m, 0) 
	else if (IsZeroMass(k) && !IsZeroMass(m1) && IsZeroMass(m2)) {
		fin = pi4sqInv * (1. + log(musq / (m1*m1)));
	} 
	// B0(0, m1, m2)
	else if (IsZeroMass(k) && !IsZeroMass(m1) && !IsZeroMass(m2)) {
		fin = 1./ (m2*m2 - m1*m1) * ( A0(m2).fin() - A0(m1).fin() );
	}
	// General case B0(k, m1, m2)
	else {
		fin = pi4sqInv * ( 0.5 * log(musq / (m1*m1)) + 0.5 * log(musq / (m2*m2)) + 2. 
			+ 0.5 * (m1*m1 - m2*m2)/k2 * log(m2*m2 / (m1*m1)) + Fb(k, m1, m2).real() );
	}

	return Integral(fin, div);
}


Renormalization::Integral Renormalization::B1(const Float &k, const Float &m1, const Float &m2) {

	// B1(0, m1, m2) = 0
	if (IsZeroMass(k)) {
		return Integral(0., 0.);
	}
	Float k2 = k*k;
	return 1. / (2.*k2) * ( (m2*m2 - m1*m1 - k2) * B0(k, m1, m2) + A0(m1) - A0(m2) );
}

Renormalization::Integral Renormalization::B00(const Float &k, const Float &m1, const Float &m2) {

	Float div = 0.;
	Float fin = 0.;
	// B00(0, m1, m2) = 1/d (A0(m2) + m1^2 B0(0, m1, m2))
	// TEST THIS CAREFULLY!! 
	if (IsZeroMass(k)) {
		Integral a0 = A0(m2);
		Integral b0 = B0(0., m1, m2);
		div = 1./4. * (a0.div() + m1*m1 * b0.div());
		fin = 1./8. * (a0.div() + 2.*a0.fin() + m1*m1 * b0.div() + 2.*m1*m1 * b0.fin() ); 
	} else {
		Integral a0 = A0(m2);
		Integral b0 = B0(k, m1, m2);
		Integral b1 = B1(k, m1, m2);
		Float k2 = k*k;
		div = 1./6. * (a0.div() + 2.*m1*m1*b0.div() + (k2 + m1*m1 - m2*m2) * b1.div());
		fin = 1./9. * (a0.div() + 2.*m1*m1*b0.div() + (k2 + m1*m1 - m2*m2) * b1.div())
			+ 1./6. * ( a0.fin() + 2.*m1*m1*b0.fin() + (k2 + m1*m1 - m2*m2) * b1.fin() );
	}

	return Integral(fin, div);
}

// Test function
void Renormalization::TestPaVe() {

	using std::cout;
	auto oldPrecision = cout.precision();
	double oldScale = this->scale;
	this->scale = 100.;
	Float smallNumber = 1e-5;
	Float k = smallNumber;
	Float m1 = 80.;
	Float m2 = 200.;

	cout.precision(15);
	cout << "========= Begin test of Passarino-Veltman functions (div ; fin) =========\n";

	cout << "k = " << k << "\n";
	cout << "B0(k, m1, m2) : " << B0(k, m1, m2).div() << " ; " <<  B0(k, m1, m2).fin() << "\n";
	cout << "B0(0, m1, m2) : " << B0(0., m1, m2).div() << " ; " <<  B0(0., m1, m2).fin() << "\n";
	cout << "B00(k, m1, m2) : " << B00(k, m1, m2).div() << " ; " <<  B00(k, m1, m2).fin() << "\n";
	cout << "B00(0, m1, m2) : " << B00(0., m1, m2).div() << " ; " <<  B00(0., m1, m2).fin() << "\n";
	

	cout << "========= End Passarino-Veltman test =========\n";
	cout.precision(oldPrecision);
	this->scale = oldScale;
}
