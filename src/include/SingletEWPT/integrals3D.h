#ifndef INTEGRALS3D_H
#define INTEGRALS3D_H

#include "common.h"

/* 3D integrals for the effective potential, in MS-bar regularization. Finite parts only.
Cannot template a namespace, so I made this is kind of a "static class". */
template<typename Complex>
class Integrals3D {

public:
	using Float = typename Complex::value_type;

private:

	Integrals3D() {}
    ~Integrals3D() {}

	// High-precision pi
	static const Float pi;

	/* Used to check if some input mass squared is practically zero. 
	Keep in mind that double is reliable only up to 15 digits. */
	static const Float massCutoff;
	static const Float pi4sqInv;
	
	// part of S3 integral that is proportional to 1/eps. Needed in other integrals that can contain (3-2eps)*S3
	static const Float divPartOfS3;
	// TODO check that I've written these parts correctly (d has 2eps, not 1eps)

public:
	static bool IsZeroMsq(Float msq) {
		return abs(msq) < massCutoff;
	}

	// J3(m^2) = 1/2 * int_p ln(p^2 + m^2) = -1/(12pi) * m^3
	static Complex J3(Float msq) {
		return -1.0/(12.0*pi) * msq * sqrt( Complex(msq) );
	}

	// I3(m^2) = int_p 1/(p^2 + m^2)
	static Complex I3(Float msq) {
		return -1.0/(4.0*pi) * sqrt( Complex(msq) );
	} 


	/* ============ 2-loop ============ */

	// S3(msq1, msq2, msq3) = int_pk 1 / [(p^2 + m1sq)(k^2 + m2sq)((p+k)^2 + m3sq) ]
	static Complex S3(Float msq1, Float msq2, Float msq3, Float scale) {

		if (IsZeroMsq(msq1) && IsZeroMsq(msq2) && IsZeroMsq(msq3)) {
			return 0.0;
		}

		// Should I take abs value here?!?
		
		Complex m1 = sqrt(Complex(msq1));
		Complex m2 = sqrt(Complex(msq2));
		Complex m3 = sqrt(Complex(msq3));
		

		// TEST
		//Complex m1 = sqrt(abs(msq1));
		//Complex m2 = sqrt(abs(msq2));
		//Complex m3 = sqrt(abs(msq3));

		return pi4sqInv * ( 0.5 + log(scale / (m1 + m2 + m3)) );
	}


	/* DVVS(msq1, msq2, msq3) = int_pk Dij(p, m2) Dkl(k, m3) dik djl / ((p+k)^2 + msq1).
	Scalar mass is m1. */
	static Complex DVVS(Float msq1, Float msq2, Float msq3, Float scale) {

		Complex res{0,0};
		/* Special cases for vanishing mass(es). Integral is symmetric wrt. msq2 <-> msq3, 
		can use this to reduce the number of special cases we need to write down. */
		if (IsZeroMsq(msq2) && !IsZeroMsq(msq3)) {
			// Swap masses so that m2 != 0, m3 == 0
			// !! Error prone
			Float temp = msq3;
			msq3 = msq2;
			msq2 = temp; 
		}

		if (IsZeroMsq(msq1) && IsZeroMsq(msq2) && IsZeroMsq(msq3)) {
			res = 0.0;

		} else if (!IsZeroMsq(msq1) && IsZeroMsq(msq2) && IsZeroMsq(msq3)) {
			res = 3*(3-1) / 4.0 * S3(msq1, 0, 0, scale) - 5.0/2.0 * divPartOfS3;
			
		} else if (IsZeroMsq(msq1) && !IsZeroMsq(msq2) && IsZeroMsq(msq3)) {
			res = 3.0 * 3*(3-1) / 4.0 * S3(msq2, 0, 0, scale) - 3*5.0/2.0 * divPartOfS3;
			
		} else if (!IsZeroMsq(msq1) && !IsZeroMsq(msq2) && IsZeroMsq(msq3)) {
			
			res = -(3.0-1.0) / (4.0*msq2) * ( (msq1 - 3*msq2) * S3(msq1, msq2, 0, scale) 
				- msq1 * S3(msq1, 0, 0, scale) + I3(msq1)*I3(msq2))
				+ 1.0 / (2*msq2) * divPartOfS3 * (-3*msq2);

		} else {

			// General result. Here the very last term is a finite contribution from (3-2eps) * S3
			res = 1.0/(4*msq2*msq3) * ( -msq3 * I3(msq1) * I3(msq2) 
				+ I3(msq3) * (-msq2 * I3(msq1) + (msq3 + msq2 - msq1)*I3(msq2))
				+ msq1*msq1 * S3(0, 0, msq1, scale) - (msq2 - msq1)*(msq2 - msq1) * S3(0, msq2, msq1, scale)
				- (msq3 - msq1)*(msq3 - msq1) * S3(msq3, 0, msq1, scale) 
				+ S3(msq3, msq2, msq1, scale) * ((msq2-msq1)*(msq2-msq1) + msq3*msq3 - 2*msq1*msq3 + 6*msq3*msq2) 
				- 8.0*msq3*msq2 * divPartOfS3
			);
		}

		return res;

	}

	/* DVSS(msq1, msq2, msq3) = int_pk (2p+k)_i (2p+k)_j Dij(k, m3) / [(p^2 + msq1)(p+k)^2 + msq2)].
	Gauge mass is m3 */
	static Complex DVSS(Float msq1, Float msq2, Float msq3, Float scale) {

		Complex res = 0.0;
		// Special cases for vanishing mass(es)
		if (IsZeroMsq(msq1) && IsZeroMsq(msq2) && IsZeroMsq(msq3)) {
			res = 0.0;
		} else if (IsZeroMsq(msq3)) {
			res = -(3.0-1.0) *  ( S3(msq1, msq2, 0, scale) * (msq1 + msq2) + I3(msq1)*I3(msq2))
				+ 2.0 * (msq1 + msq2) * divPartOfS3;
		} else {
			// General result
			res = 1.0 / msq3 * ( (msq3 + msq2 - msq1) * I3(msq2) * I3(msq3) 
				+ I3(msq1) * (-msq3 * I3(msq2) + I3(msq3) * (msq1 - msq2 + msq3))
				- (msq1 - msq2)*(msq1 - msq2) * S3(msq1, 0, msq2, scale) 
				+ (msq1*msq1 + msq2*msq2 + msq3*msq3 - 2*msq1*msq2 - 2*msq1*msq3 - 2*msq2*msq3) 
				* S3(msq1, msq3, msq2, scale)
			);
		}

		return res;
	}

	/* DVVV(msq1, msq2, msq3) is 4 times the integral in eq. (119) of hep-ph/9404201 */
	static Complex DVVV(Float msq1, Float msq2, Float msq3, Float scale) {
		Complex res = 0.0;

		// Special cases for zero or equal masses (in practice these are all that are needed)
		if (IsZeroMsq(msq1) && IsZeroMsq(msq2) && IsZeroMsq(msq3)) {
			res = 0.0;

		} else if (IsZeroMsq(msq1 - msq2) && IsZeroMsq(msq3)) {
			// msq1 == msq2 and msq3 == 0
			/* Need to check for complex logs here => I write this in an 'obscure' form using S3
			* to let the S3 function do the checking. */
			res = msq1/(12.0*pi*pi) * (89.0/16.0 + log(64.0)) - 10*msq1 * S3(msq1, 0, 0, scale);

		} else if (IsZeroMsq(msq1 - msq2) && !IsZeroMsq(msq3)) {
			// msq1 == msq2 and msq3 != 0
			res = -I3(msq1)*I3(msq3) / (6.0*msq1*msq3) * (3*msq1*msq1 - 11*msq1*msq3 - 15*msq3*msq3)
				+ 1.0/(192*pi*pi*msq1) * (56*msq1*msq1 - 30*msq1*msq3 - 3*msq3*msq3)
				- (msq1-msq3)*(msq1-msq3)*(msq1*msq1 + 6*msq1*msq3 + msq3*msq3) / (2*msq1*msq1*msq3) * S3(msq1, msq3, 0, scale) 
				+ 4 *(msq1-msq3)*(msq1-msq3) / msq1 * divPartOfS3
				- (4*msq1-msq3) * (8*msq1*msq1 + 12*msq1*msq3 + msq3*msq3) / (4*msq1*msq1) * S3(msq1, msq1, msq3, scale)
				+ 2.0/msq1 * (4*msq1 - msq3)*(msq1 + 2*msq3) * divPartOfS3
				+ msq3*msq3*msq3 / (4*msq1*msq1) * S3(msq3, 0, 0, scale) + msq1*msq1/(2*msq3) * S3(msq1, 0, 0, scale);

		} else {
			// Generic result: copypasted from Mathematica with CForm
			res = (3*msq1)/(32.*(pi*pi)) + (3*msq2)/(32.*(pi*pi)) + (3*msq3)/(32.*(pi*pi)) + 
				(13.*I3(msq1)*I3(msq2))/6. + (5*msq1*I3(msq1)*I3(msq2))/(4.*msq2) + 
				(5*msq2*I3(msq1)*I3(msq2))/(4.*msq1) - (5*msq3*I3(msq1)*I3(msq2))/(4.*msq1) - 
				(5*msq3*I3(msq1)*I3(msq2))/(4.*msq2) - 
				(msq3*msq3*I3(msq1)*I3(msq2))/(4.*msq1*msq2) + (13.*I3(msq1)*I3(msq3))/6. - 
				(5*msq2*I3(msq1)*I3(msq3))/(4.*msq1) + (5*msq1*I3(msq1)*I3(msq3))/(4.*msq3) - 
				(5*msq2*I3(msq1)*I3(msq3))/(4.*msq3) - 
				(msq2*msq2*I3(msq1)*I3(msq3))/(4.*msq1*msq3) + 
				(5*msq3*I3(msq1)*I3(msq3))/(4.*msq1) + (13.*I3(msq2)*I3(msq3))/6. - 
				(5*msq1*I3(msq2)*I3(msq3))/(4.*msq2) - (5.*msq1*I3(msq2)*I3(msq3))/(4.*msq3) - 
				(msq1*msq1*I3(msq2)*I3(msq3))/(4.*msq2*msq3) + 
				(5*msq2*I3(msq2)*I3(msq3))/(4.*msq3) + (5.*msq3*I3(msq2)*I3(msq3))/(4.*msq2) + 
				(msq1*msq1*msq1*S3(msq1,0,0,scale))/(4.*msq2*msq3) - 
				(msq1*msq1*S3(msq1,msq2,0,scale))/msq3 - 
				(msq1*msq1*msq1*S3(msq1,msq2,0,scale))/(4.*msq2*msq3) + 
				(5*msq1*msq2*S3(msq1,msq2,0,scale))/(2.*msq3) - 
				(msq2*msq2*S3(msq1,msq2,0,scale))/msq3 - 
				(msq2*msq2*msq2*S3(msq1,msq2,0,scale))/(4.*msq1*msq3) - 
				5*msq1*S3(msq1,msq2,msq3,scale) + (msq1*msq1*S3(msq1,msq2,msq3,scale))/msq2 - 
				5*msq2*S3(msq1,msq2,msq3,scale) + (msq2*msq2*S3(msq1,msq2,msq3,scale))/msq1 + 
				(msq1*msq1*S3(msq1,msq2,msq3,scale))/msq3 + 
				(msq1*msq1*msq1*S3(msq1,msq2,msq3,scale))/(4.*msq2*msq3) - 
				(5*msq1*msq2*S3(msq1,msq2,msq3,scale))/(2.*msq3) + 
				(msq2*msq2*S3(msq1,msq2,msq3,scale))/msq3 + 
				(msq2*msq2*msq2*S3(msq1,msq2,msq3,scale))/(4.*msq1*msq3) - 
				5*msq3*S3(msq1,msq2,msq3,scale) - 
				(5*msq1*msq3*S3(msq1,msq2,msq3,scale))/(2.*msq2) - 
				(5*msq2*msq3*S3(msq1,msq2,msq3,scale))/(2.*msq1) + 
				(msq3*msq3*S3(msq1,msq2,msq3,scale))/msq1 + 
				(msq3*msq3*S3(msq1,msq2,msq3,scale))/msq2 + 
				(msq3*msq3*msq3*S3(msq1,msq2,msq3,scale))/(4.*msq1*msq2) - 
				(msq1*msq1*S3(msq1,msq3,0,scale))/msq2 - 
				(msq1*msq1*msq1*S3(msq1,msq3,0,scale))/(4.*msq2*msq3) + 
				(5*msq1*msq3*S3(msq1,msq3,0,scale))/(2.*msq2) - 
				(msq3*msq3*S3(msq1,msq3,0,scale))/msq2 - 
				(msq3*msq3*msq3*S3(msq1,msq3,0,scale))/(4.*msq1*msq2) + 
				(msq2*msq2*msq2*S3(msq2,0,0,scale))/(4.*msq1*msq3) - 
				(msq2*msq2*S3(msq2,msq3,0,scale))/msq1 - 
				(msq2*msq2*msq2*S3(msq2,msq3,0,scale))/(4.*msq1*msq3) + 
				(5*msq2*msq3*S3(msq2,msq3,0,scale))/(2.*msq1) - 
				(msq3*msq3*S3(msq2,msq3,0,scale))/msq1 - 
				(msq3*msq3*msq3*S3(msq2,msq3,0,scale))/(4.*msq1*msq2) + 
				(msq3*msq3*msq3*S3(msq3,0,0,scale))/(4.*msq1*msq2);
		}

		return res;
	} // end DVVV

	/* DVV(msq1, msq2), figure-eight integral with gauge propagators */
	static Complex DVV(Float msq1, Float msq2) {
		return 8.0/3.0 * I3(msq1)*I3(msq2);
	}

	/* DVVG(msq), gauge-ghost-ghost vacuum bubble */
	static Complex DVGG(Float msq, Float scale) {
		return 1.0/4.0 * msq * S3(msq, 0, 0, scale);
	}


	// Test integral results
	static void TestIntegrals() {

		Float msq1 = (Float) pow(142.85421, 2);
		Float msq2 = (Float) pow(94.11231, 2);
		Float msq3 = (Float) pow(112.95122, 2);
		Float scale = (Float) 541.51223;

		Float testSmallNumber = (Float)1e-11;

		using std::cout;

		cout.precision(12);
		cout << "I3(msq1): " << I3(msq1) << "\n";

		cout << "S3(msq1, msq2, msq3): " << S3(msq1, msq2, msq3, scale) << "\n";
		cout << "S3(0, 0, 0): " << S3(0, 0, 0, scale) << "\n";
		cout << "S3(1e-11, 1-e11, 1e-11): " << S3(testSmallNumber, testSmallNumber, testSmallNumber, scale) << "\n";

		cout << "DVVS(msq1, msq2, msq3): " << DVVS(msq1, msq2, msq3, scale) << "\n";
		cout << "DVVS(msq1, msq2, 0): " << DVVS(msq1, msq2, 0, scale) << "\n";
		cout << "DVVS(msq1, 0, 0): " << DVVS(msq1, 0, 0, scale) << "\n";
		cout << "DVVS(msq1, 0, msq3): " << DVVS(msq1, 0, msq3, scale) << "\n";
		cout << "DVVS(0, 0, msq3): " << DVVS(0, 0, msq3, scale) << "\n";
		cout << "DVVS(0, msq2, msq3): " << DVVS(0, msq2, msq3, scale) << "\n";
		cout << "DVVS(0, msq2, 0): " << DVVS(0, msq2, 0, scale) << "\n";
		cout << "DVVS(0, 0, 0): " << DVVS(0, 0, 0, scale) << "\n";
		cout << "DVVS(1e-11, 1e-11, 1e-11): " << DVVS(testSmallNumber, testSmallNumber, testSmallNumber, scale) << "\n";

		cout << "DVSS(msq1, msq2, msq3): " << DVSS(msq1, msq2, msq3, scale) << "\n";
		cout << "DVSS(msq1, msq2, 0): " << DVSS(msq1, msq2, 0, scale) << "\n";
		cout << "DVSS(msq1, 0, msq3): " << DVSS(msq1, 0, msq3, scale) << "\n";
		cout << "DVSS(0, msq2, msq3): " << DVSS(0, msq2, msq3, scale) << "\n";
		cout << "DVSS(msq1, 0, 0): " << DVSS(msq1, 0, 0, scale) << "\n";
		cout << "DVSS(0, msq2, 0): " << DVSS(0, msq2, 0, scale) << "\n";
		cout << "DVSS(0, 0, msq3): " << DVSS(0, 0, msq3, scale) << "\n";
		cout << "DVSS(0, 0, 0): " << DVSS(0, 0, 0, scale) << "\n";
		cout << "DVSS(1e-11, 1e-11, 1e-11): " << DVSS(testSmallNumber, testSmallNumber, testSmallNumber, scale) << "\n";

		cout << "DVVV(msq1, msq2, msq3): " << DVVV(msq1, msq2, msq3, scale) << "\n";
		cout << "DVVV(msq1, msq1, msq3): " << DVVV(msq1, msq1, msq3, scale) << "\n";
		cout << "DVVV(msq1, msq1, 0): " << DVVV(msq1, msq1, 0, scale) << "\n";
		cout << "DVVV(0, 0, 0): " << DVVV(0, 0, 0, scale) << "\n";
		cout << "DVVV(1e-11, 1e-11, 1e-11): " << DVVV(testSmallNumber, testSmallNumber, testSmallNumber, scale) << "\n";
	};

};

/* Set the constants. Need funny typename syntax here because the Float type is defined within the template class. */
template <typename Complex>
const typename Integrals3D<Complex>::Float Integrals3D<Complex>::pi = PI; // TODO precision
template <typename Complex>
const typename Integrals3D<Complex>::Float Integrals3D<Complex>::pi4sqInv = 1.0/(16.0 * pi*pi);
template <typename Complex>
const typename Integrals3D<Complex>::Float Integrals3D<Complex>::divPartOfS3 = 1.0 / (64.0*pi*pi);
template <typename Complex>
const typename Integrals3D<Complex>::Float Integrals3D<Complex>::massCutoff = 1e-10;


#endif