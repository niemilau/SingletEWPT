#ifndef EFFPOT_T0_H
#define EFFPOT_T0_H

#include "common.h"
#include "effpot.h"

// dlib library provides methods for finding function minima 
#include <dlib/optimization.h>

/* Calculate T=0 effective potential 3D SM + singlet, in Rxi Landau gauge.
* The eff potential is a function of two background fields, Veff = Veff(v,x). 
* These are defined by shifts phi -> phi + 1/sqrt(2) (0, v) and S -> S + x. */ 

// Loop order for the potential. Using a separate typename in case this gets extended later 
using ELoopOrderVeff = ELoopOrder;

// Much of this is copy-pasted from effpot.h. Would be better to use eg. interfaces, but making it work nicely would require
// a lot of code refactorization... 


/* Class for computing the potential using user-specified complex number type (double, float, multiprecision etc) */
template <typename Complex>
class EffPotT0 {

public:
	using Float = typename Complex::value_type;

private:
	// Background fields; these variables are used in all internal calculations
	Float v, x; 
	// Field-dependent mixing angle
	Float theta;
	
	Float ytsq;
	Float g1sq, g2sq, lambda, msqPhi;
	Float b1, b2, b3, b4, a1, a2;
	Float scale;

	Float mh1sq, mh2sq, mGsq, mGpmsq;
	Float mWsq, mZsq;
	// Top quark mass (field dependent)
	Float mtsq;
	
	bool bIsZ2Symmetric;

public:

	// Keep track of possible issues in computation of the potential
	int warnings = 0;

	/* Constructor that sets the couplings. (ParameterMap only has doubles, so precision loss can occur!) */
	EffPotT0(const ParameterMap &params3D) {
		SetParams(params3D);
		v = 0;
		x = 0;
		theta = 0;
	}


	void SetParams(const ParameterMap &params) {
		scale = GetFromMap(params, "RGScale");
		ytsq = GetFromMap(params, "ytsq");
		g1sq = GetFromMap(params, "g1sq");
		g2sq = GetFromMap(params, "g2sq");
		lambda = GetFromMap(params, "lambda");
		msqPhi = GetFromMap(params, "msqPhi");

		b1 = GetFromMap(params, "b1");
		b2 = GetFromMap(params, "b2");
		b3 = GetFromMap(params, "b3");
		b4 = GetFromMap(params, "b4");
		a1 = GetFromMap(params, "a1");
		a2 = GetFromMap(params, "a2");

		// Z2 symmetric model?
		double smallNumber = 1e-6;
		bIsZ2Symmetric = (abs(b1) < smallNumber && abs(b3) < smallNumber && abs(a1) < smallNumber); 
	}

	// Fix background fields and calculate mass eigenvalues and theta
	void SetBackgroundFields(Float v, Float x) {
		this->v = v;
		this->x = x;
		CalcMixingAngle();
		CalcMassEigenvalues();
	}

	// Calculate field-dependent mixing angle and store in 'theta'
	void CalcMixingAngle();

	/* Calculate eigenvalues of scalar and gauge mass matrices (squared masses)
	and store them in member variables (with obvious names) */
	void CalcMassEigenvalues();

	// Return the presently stored mixing angle
	inline Float GetMixingAngle() { return theta; }

	// Return scalar and gauge masses squared
	ParameterMap GetMassEigenvalues();

	/* Tree level potential in conventions of 2103.07467. */
	Float V0();

	// 1-loop correction to Veff. Can be complex outside minima.
	Complex V1();

	/* Evaluate Veff(v, x) to a given loop accuracy. */
	Complex EvaluatePotential(Float v, Float x, const ELoopOrderVeff loopOrder);

	
	/* Evaluate the potential with double valued fields instead and return a complex double. 
	Intermediate computations still use Floats. Needed for dlib minimization routines.
	I could not use overloading here because if Float = double, the compiler reports a conflicting overload. */
	std::complex<double> EvaluatePotentialAsDouble(double v, double x, const ELoopOrderVeff loopOrder) {
		Complex Veff = EvaluatePotential((Float)v, (Float)x, loopOrder);
		return (std::complex<double>)Veff;
	}

public:

	/* Analytically finds all extrema of the tree-level potential so that v,x are real and v >= 0. 
	Does not distinguish between minima/maxima but does not include saddle points. */
	std::vector<std::array<double, 2>> TreeLevelMinima() const;

	/* Global minimization. Looks for several local minima (based on intuitive guesses) and takes the deepest of those. 
	Works with doubles, which get converted to Complex for internal computations */
	ParameterMap FindGlobalMinimum(const ELoopOrderVeff loopOrder);
	
	// Find a local minimum with initial guess (v0, x0). Returns doubles.
	ParameterMap FindLocalMinimum(const ELoopOrderVeff loopOrder, double v0, double x0, const MinimizationParams &minParams);

	// Call FindLocalMinimum with default MinimizationParams struct
	inline ParameterMap FindLocalMinimum(const ELoopOrderVeff loopOrder, double v0, double x0) {
		return FindLocalMinimum(loopOrder, v0, x0, MinimizationParams());	
	}

private:
	/* Return set of (v,x) pairs to use as starting points for FindGlobalMinimum(). 
	This includes at least the tree-level extrema of the potential 
	and possibly some hand-picked points corresponding to symmetric, Higgs, singlet phases (if not included in the former) */ 
	std::vector<std::array<double, 2>> InitialSearchPoints() const;


	// 4D integral for 1-loop Veff: J4(m^2) = 1/2 \int \ln(p^2 + m^2). Finite part only
	Complex J4(Float msq) {
		if (abs(msq) < 1e-8) return 0;
		
		return -msq*msq / (64.*PI*PI) * ( 3./2. + log(scale*scale / Complex(msq) ) );
	}

};

#include "effpot_T0.tpp"

#endif