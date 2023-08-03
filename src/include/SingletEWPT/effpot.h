#ifndef EFFPOT_H
#define EFFPOT_H

#include "common.h"
#include "integrals3D.h"

// dlib library provides methods for finding function minima 
#include <dlib/optimization.h>


/* Calculate effective potential for 3D SM + singlet, in Rxi Landau gauge.
* The eff potential is a function of two background fields, Veff = Veff(v,x). 
* These are defined by shifts phi -> phi + 1/sqrt(2) (0, v) and S -> S + x. */ 

// Loop order for the 3D potential. Using a separate typename in case this gets extended later 
using ELoopOrderVeff = ELoopOrder;


// parameters for BOBYQA minimization
struct MinimizationParams {

	// These values seem to work pretty well in generic situations
	double initialTrustRadius = 10;
	double stoppingTrustRadius = 1e-4;
	double maxFunctionEvaluations = 1000;
};


/* Class for computing the potential using user-specified complex number type (double, float, multiprecision etc) */
template <typename Complex>
class EffPot {

public:
	using Float = typename Complex::value_type;

private:
	// Background fields; these variables are used in all internal calculations
	Float v, x; 
	// Field-dependent mixing angle
	Float theta;
	
	Float g1sq, g2sq, lambda, msqPhi;
	Float b1, b2, b3, b4, a1, a2;
	Float scale;

	Float mh1sq, mh2sq, mGsq, mGpmsq;
	Float mWsq, mZsq;

	// dim5 and 6
	Float c05, c23, c41, c60, c06, c42, c24;

	// Flag for checking if our theory has Z2 symmetry. Set in SetParams()
	bool bIsZ2Symmetric = false;

public:

	// Keep track of possible issues in computation of the potential
	int warnings = 0;
	// Start looking for issues when v is smaller than this
	Float vWarningThreshold = 1e-5;

	// Class-specific typedef for integrals of correct type
	typedef Integrals3D<Complex> Integrals;

	/* Constructor that sets the couplings. (ParameterMap only has doubles, so precision loss can occur!) */
	EffPot(const ParameterMap &params3D) {
		SetParams(params3D);
		v = 0;
		x = 0;
		theta = 0;
	}


	void SetParams(const ParameterMap &params3D) {
		scale = GetFromMap(params3D, "RGScale");
		g1sq = GetFromMap(params3D, "g1sq");
		g2sq = GetFromMap(params3D, "g2sq");
		lambda = GetFromMap(params3D, "lambda");
		msqPhi = GetFromMap(params3D, "msqPhi");

		b1 = GetFromMap(params3D, "b1");
		b2 = GetFromMap(params3D, "b2");
		b3 = GetFromMap(params3D, "b3");
		b4 = GetFromMap(params3D, "b4");
		a1 = GetFromMap(params3D, "a1");
		a2 = GetFromMap(params3D, "a2");

		c05 = GetFromMap(params3D, "c05");
		c23 = GetFromMap(params3D, "c23");
		c41 = GetFromMap(params3D, "c41");
		c60 = GetFromMap(params3D, "c60");
		c06 = GetFromMap(params3D, "c06");
		c42 = GetFromMap(params3D, "c42");
		c24 = GetFromMap(params3D, "c24");

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

	/* Tree level potential in conventions of 2103.07467.
	If bDoDim6 == true, includes dimension 5 and 6 operators. */
	Float V0(bool bDoDim6);

	// 1-loop correction to Veff. Can be complex outside minima.
	Complex V1();

	// 2-loop correction to Veff. Defined in effpot_2loop.cpp.
	Complex V2();

	/* Evaluate Veff(v, x) to a given loop accuracy. */
	Complex EvaluatePotential(Float v, Float x, const ELoopOrderVeff loopOrder, bool bDoDim6);

	
	/* Evaluate the potential with double valued fields instead and return a complex double. 
	Intermediate computations still use Floats. Needed for dlib minimization routines.
	I could not use overloading here because if Float = double, the compiler reports a conflicting overload. */
	std::complex<double> EvaluatePotentialAsDouble(double v, double x, const ELoopOrderVeff loopOrder, bool bDoDim6) {
		Complex Veff = EvaluatePotential((Float)v, (Float)x, loopOrder, bDoDim6);
		return (std::complex<double>)Veff;
	}
	

private:
	// Struct for passing data to wrapper functions
	struct DataWrapper {
		ELoopOrderVeff loopOrder;
		bool bDoDim6 = false;
		EffPot<Complex> *ptrToObject;
	};

public:
	/* Static wrapper function for NLOPT minimization routines. x contains field values, grad contains gradient info (not used ATM)
	and data will be cast to */
	static double EvaluatePotentialWrapper(const std::vector<double>& x, std::vector<double>& grad, void* data) {
		(void)grad;
		DataWrapper* dataWrapper = static_cast<DataWrapper*>(data);
		return real( dataWrapper->ptrToObject->EvaluatePotentialAsDouble(x[0], x[1], dataWrapper->loopOrder, dataWrapper->bDoDim6) );
	}

	/* Calculate dv/v, dx/x, where dv is the shift to tree-level VEV caused by dimension 5,6 operators.
	The inputs should be tree-level values for the minima.
	If v=0 or x=0 (meaning small), sets the shift to 0. First element in the returned vector is dv/v and second is dx/x. */
	std::vector<double> FieldShiftsDim6(Float v, Float x);


	/* Analytically finds all minima of the tree-level potential so that v,x are real and v >= 0. */
	std::vector<std::array<double, 2>> TreeLevelMinima() const;

	/* Global minimization. Looks for several local minima (based on intuitive guesses) and takes the deepest of those. 
	Works with doubles, which get converted to Complex for internal computations */
	ParameterMap FindGlobalMinimum(const ELoopOrderVeff loopOrder, bool bDoDim6);
	
	// Find a local minimum with initial guess (v0, x0). Returns doubles.
	ParameterMap FindLocalMinimum(const ELoopOrderVeff loopOrder, bool bDoDim6, double v0, double x0, const MinimizationParams &minParams);

	// Call FindLocalMinimum with default MinimizationParams struct
	inline ParameterMap FindLocalMinimum(const ELoopOrderVeff loopOrder, bool bDoDim6, double v0, double x0) {
		return FindLocalMinimum(loopOrder, bDoDim6, v0, x0, MinimizationParams());	
	}

private:
	/* Return set of (v,x) pairs to use as starting points for FindGlobalMinimum(). 
	This includes at least the tree-level extrema of the potential 
	and possibly some hand-picked points corresponding to symmetric, Higgs, singlet phases (if not included in the former) */ 
	std::vector<std::array<double, 2>> InitialSearchPoints() const;

public:

	/* 2-loop diagrams. These are defined in effpot_2loop.tpp */

	// Non-charged scalar figure-8 diagrams
	Complex NeutralScalar8();
	// Charged scalar figure-8 diagrams
	Complex ChargedScalar8();
	// Mixed charged and neutral scalar figure-8 diagrams
	Complex MixedScalar8();
	// Pure gauge figure-8 diagrams
	Complex PureGauge8();
	// Mixed gauge-scalar figure-8 diagrams
	Complex MixedGaugeScalar8();

	// Neutral scalar sunset diagrams
	Complex NeutralScalarSunset();
	// Mixed neutral-charged scalar sunset diagrams
	Complex MixedScalarSunset();
	// Sunset diagrams with 2 gauge, 1 scalar lines
	Complex VVS();
	// Sunset diagrams with 1 gauge, 2 scalar lines
	Complex VSS();
	// Pure gauge sunset diagrams
	Complex VVV();
	// Ghost-gauge sunset diagrams
	Complex VGG();

	// Sunset diagrams with 2 gauge, 1 scalar lines
	Complex VVStest();

	// Test function: print out numerical results for individual diagram types
	void TestDiagrams(Float v, Float x);

	// Test v->0 limit. Prints abs of (Veff(v) - Veff(0))/Veff(0) as function of small v to outfile.
	void TestSmallLimit(std::string outfile);
};


// Include headers containing class implementation (generic template implementations need to be in headers)
#include "effpot.tpp"
#include "effpot_2loop.tpp"

#endif