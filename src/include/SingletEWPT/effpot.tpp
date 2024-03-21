#include <gsl/gsl_poly.h>




template <typename Complex>
Complex EffPot<Complex>::EvaluatePotential(Float v, Float x, const ELoopOrderVeff loopOrder, bool bDoDim6)
{
	Complex Veff{0, 0};

	this->SetBackgroundFields(v, x);
	// Now mass eigenvalues have been calculated

	int order = static_cast<int>(loopOrder);

	Veff += this->V0(bDoDim6);
	
	if (order > 0) {
		Veff += this->V1();
	}
	if (order > 1) {
		Veff += this->V2();
	}
	
	return Veff;
}

template <typename Complex>
void EffPot<Complex>::CalcMixingAngle() {

	/* Mixing angle is defined as the rotation angle that diagonalizes mass matrix
	of neutral Higgs modes. If the mass matrix is M = (A, C // C, B) and 
	P = ( cos(theta), sin(theta) // -sin(theta) cos(theta) ) is a rotation matrix, then
	P.M.P^T should be diagonal. 
	More specifically, M is chosen so that the quadratic tree potential is 
	V = 1/2 (h, S).M.(h // S) where h is the SU(2) doublet component that got shifted 
	and S is the singlet field. The diagonal eigenstates are (h1 // h2) = P.(h // S)
	==> V = 1/2 (h1, h2) P.M.P^T (h1 // h2).
	*/

	/* Here I always choose the first branch of arctan (to be consistent with the vertex
	rules copy-pasted from Mathematica). A concequence is that for A > B, h1 is actually the
	heavier eigenstate and h2 is lighter (see CalcMassEigenstates()). 
	If we wanted, we could shift theta -> theta + n*pi/2 with n = odd to interchange
	the mass eigenstates, but this would require rewriting of the vertex rules.
	*/

	/* NB! Equation (B3) in 2103.07467 corresponds to mixing angle with the opposite sign
	of what I use here (choice in the paper matches that used in T=0 calculation). 
	Similarly the masses in eqs. (B1)-(B2) and the list of couplings after (B13) have 
	opposite sign for theta than what I have in this code. */

	Float v2 = v*v;
	Float x2 = x*x;

	Float A = msqPhi + 3*lambda*v2 + a1*x / 2. + a2*x2 / 2.;
	
	Float B = b2 + 3*b4*x2 + 2*b3*x + a2*v2 / 2.;

	Float C = a1*v/2. + a2*v*x;

	//Float smallNumber = 1e-12;

	Float d = C / (A-B);

	theta = 0.5 * atan(2.0*d);
	/*
	if (A - B == 0) {
		if (C > 0) theta =  PI/4.;
		else theta = -PI/4.;
	} else {
		
		theta = 0.5 * atan(2*C / (A - B));	
	}
	*/
}

template <typename Complex>
void EffPot<Complex>::CalcMassEigenvalues() 
{

	Float v2 = v*v;
	Float x2 = x*x;


	// Gauge bosons
	this->mWsq = 1.0/4.0 * g2sq * v2;
	this->mZsq = 1.0/4.0 * (g1sq + g2sq) * v2;


	// Neutral Goldstone mode
	Float mGsq = msqPhi + lambda*v2 + 0.5*a1*x + 0.5*a2*x2;
	this->mGsq = mGsq;
	
	// Charged Goldstone mode
	Float mGpmsq = mGsq;
	this->mGpmsq = mGpmsq;


	/* For neutral Higgs modes we need to diagonalize mass matrix of form (A, C // C, B).
	Furthermore, for the 2-loop calculation to work we need the eigenvalues mh1 and mh2
	to match the eigenstates h1, h2 used in vertex rules. h1 and h2 are obtained by rotating
	the original Higgs and singlet modes by a rotation given by the mixing angle 
	==> use this rotation here to get correct mh1, mh2, 
	instead of assuming that we'd always have mh2 >= mh1. */

	Float A = msqPhi + 3*v2*lambda + 0.5*a1*x + 0.5*a2*x2;
	Float B = b2 + 3*b4*x2 + 2*b3*x + 0.5*a2*v2;
	Float C = 0.5*a1*v + a2*v*x;

	Float st = sin(theta);
	Float ct = cos(theta);

	Float mh1sq = A*ct*ct + B*st*st + 2.* C *st*ct;
	Float mh2sq = A*st*st + B*ct*ct - 2.* C *st*ct;

	// Can check that for A > B one gets mh1sq > mh2sq.


	/* // Old calculation that doesn't care about mixing angle => always mh2 >= mh1, 
	   // which does not work at 2-loop because of how the vertex rules are set up!! 
	double A = msqPhi + 3*v2*lambda + 0.5*a1*x + 0.5*a2*x2;
	double B = b2 + 3*b4*x2 + 2*b3*x + 0.5*a2*v2;
	double C = 0.5*a1*v + a2*v*x;

	// This is always >= 0 => no issues with taking square root
	double DD = A*A - 2*A*B + B*B + 4*C*C;

	// Lighter Higgs mode
	double mh1sq = 0.5 * (A + B - sqrt(DD));
	// Heavier Higgs mode
	double mh2sq = 0.5 * (A + B + sqrt(DD));
	*/

	this->mh1sq = mh1sq;
	this->mh2sq = mh2sq;

/*
	// Need special cutoff for mh1sq, mh2sq to avoid log divergence?!
	Float smallMassSq = 100;
	if (abs(mh1sq) < smallMassSq )
		this->mh1sq  = smallMassSq;

	if (abs(mh2sq) < smallMassSq )
		this->mh2sq = smallMassSq;
*/

}

// need funny typename syntax here because the return type is defined within the template class
template <typename Complex>
typename EffPot<Complex>::Float EffPot<Complex>::V0(bool bDoDim6) 
{
	Float v2 = v*v;
	Float v4 = v2*v2;
	Float x2 = x*x;
	Float x3 = x*x*x;
	Float x4 = x*x*x*x;

	Float res = 0.5*msqPhi*v2 + 0.25*lambda*v4 + 0.25*a1*v2*x + 0.25*a2 * v2*x2
				+ b1*x + 0.5*b2 * x2 + 1.0/3.0*b3 * x3 + 0.25*b4 *x4; 

	if (bDoDim6) {

		Float x5 = x*x*x*x*x;
		Float x6 = x5*x;

		res += c05 * x5 + 0.5*c23 * v2*x3 + 0.25*c41 * v4*x + 1.0/8.0 *c60 * v2*v2*v2
				+ c06 * x6 + 0.25*c42 *v4*x2 + 0.5*c24 * v2*x4;
	}

	return res;
}

template <typename Complex>
Complex EffPot<Complex>::V1()
{
	// Alias for the loop integral of correct type. Cannot use 'using' because it is a static member function
	const auto J3 = Integrals::J3;

	Float d = 3;
	Complex res{0,0};

	// Gauge + ghost loops
	res += 2.0*(d-1.0) * J3(mWsq) + (d-1.0) * J3(mZsq); 
	// Scalar loops
	res += J3(mh1sq) + J3(mh2sq) + J3(mGsq) + 2.0*J3(mGpmsq);

	return res;
} 


template <typename Complex>
std::vector<std::array<double, 2>> EffPot<Complex>::TreeLevelMinima() const {

	// Need to solve dV/dv = 0, dV/dx = 0 so that d^2V/dv^2 > 0, hessian > 0
	std::vector<std::array<double, 2>> minima;
	minima.reserve(6);
	// Will solve 'manually' so that v is first solved from dV/dv = 0 
	// => dV/dx = 0 becomes cubic eq: 	a x^3 + b x^2 + c x + d = 0. 
	// Note that gsl_poly_solve_cubic() solves x^3 + b x^2 + c x + d = 0, so need to scale with coeff. of x^3.

	// NB: v > 0 case only makes sense for lambda > 0 but our loop-corrected lambda may actually be slightly negative.
	// To allow for this case I just take abs of lambda here.
	double lambda_abs = std::abs(lambda);	
	
	// Finds real roots of a x^3 + b x^2 + c x + d = 0. Has either 1 or 3 real roots
	auto PolySolve3 = [](double a, double b, double c, double d) {
		double x[3];
		std::vector<double> res;
		res.reserve(3);
		int nroots = 0;

		if (std::abs(a) < 1e-14) {
			// now essentially a = 0 => is actually a quadratic eq
			nroots = gsl_poly_solve_quadratic(b, c, d, &x[0], &x[1]);

		} else {
			// Now we have a proper cubic eq and can use gsl_poly_solve_cubic
			nroots = gsl_poly_solve_cubic(b/a, c/a, d/a, &x[0], &x[1], &x[2]);
		}

		for (int i=0; i<nroots; ++i) {
			res.push_back(x[i]);
		}

		return res;
	};

	// Hessian determinant
	auto Hessian = [&](double v, double x) {
		
		return b2*msqPhi - (a1*a1*(v*v))/4. + (a2*msqPhi*(v*v))/2. + (a1*b2*x)/2. + 
				2*b3*msqPhi*x - (3*a1*a2*(v*v)*x)/4. + (a2*b2*(x*x))/2. + 
				a1*b3*(x*x) + 3*b4*msqPhi*(x*x) - (3*(a2*a2)*(v*v)*(x*x))/4. + 
				a2*b3*(x*x*x) + (3*a1*b4*(x*x*x))/2. + 
				(3*a2*b4*(x*x*x*x))/2. + 3*b2*(v*v)*lambda_abs + 
				(3*a2*(v*v*v*v)*lambda_abs)/2. + 6*b3*(v*v)*x*lambda_abs + 9*b4*(v*v)*(x*x)*lambda_abs;
	};
	// d^2V / dv^2
	auto d2Vdv2 = [&](double v, double x) {
		return msqPhi + 3.0*v*v * lambda_abs + 0.5*a1*x + 0.5*a2*x*x; 
	};
	
	// Solve cubic eq for x
	double a, b, c, d;
	// v = 0 case
	a = b4;
	b = b3;
	c = b2;
	d = b1;
	std::vector<double> xsolns = PolySolve3(a, b, c, d);
	for (double x : xsolns) {
		double v = 0;
		if (Hessian(v, x) > 0 && d2Vdv2(v, x) > 0) {
			minima.push_back({v, x});
		}
	}

	// v > 0 case. 
	a = b4 - a2*a2/(4.0*lambda_abs);
	b = b3 - 3.0*a1*a2 / (8.0*lambda_abs);
	c = b2 - a1*a1 / (8.0*lambda_abs) - a2*msqPhi / (2.0*lambda_abs);
	d = b1 - a1*msqPhi / (4.0*lambda_abs);

	xsolns = PolySolve3(a, b, c, d);
	for (double x : xsolns) {

		// Check if real soln for v can exist
		double kappa = -2.0*msqPhi - a1*x - a2*x*x;
		if (kappa < 0 && lambda > 0) continue;

		double v = std::sqrt( kappa / (2.0*lambda_abs));
		if (Hessian(v, x) > 0 && d2Vdv2(v, x) > 0) {
			minima.push_back({v, x});
		}
	}

	return minima;
}

template <typename Complex>
std::vector<std::array<double, 2>> EffPot<Complex>::InitialSearchPoints() const {

	std::vector<std::array<double, 2>> searchPoints = TreeLevelMinima();

	bool bHasHiggsPhase = false;
	bool bHasSingletPhase = false;
	bool bHasSymmetricPhase = false;

	if (bIsZ2Symmetric) {
		// In Z2 symmetric limit we just need x >= 0
		for (std::size_t i = 0; i < searchPoints.size(); ++i) {
			auto point = searchPoints[i];
			double x = point[1];

			if (x < 0) {
				searchPoints.erase(searchPoints.begin() + i);
			}
		}
	}

	// Check that we have all sensible points
	double smallFieldValue = 1e-3;
	for (const auto& point: searchPoints) {
		double v = point[0];
		double x = point[1];

		if (std::abs(v) > smallFieldValue) bHasHiggsPhase = true;
		if (std::abs(x) > smallFieldValue) bHasSingletPhase = true;
		if (std::abs(v) <= smallFieldValue && std::abs(x) <= smallFieldValue) bHasSymmetricPhase = true;
		
	}
	
	if (!bHasSymmetricPhase) {
		searchPoints.push_back({smallFieldValue, smallFieldValue});
	}
	if (!bHasHiggsPhase) {
		searchPoints.push_back({50.0, smallFieldValue});
	}
	if (!bHasSingletPhase) {
		searchPoints.push_back({smallFieldValue, 50.0});
	}

	return searchPoints;
}


// Find a local minimum with initial guess (v0, x0). Returns doubles.
template <typename Complex>
ParameterMap EffPot<Complex>::FindLocalMinimum(const ELoopOrderVeff loopOrder, bool bDoDim6, double v0, double x0, const MinimizationParams &minParams) {

	using ColumnVector = dlib::matrix<double, 0, 1>;

	// Lambda for passing real part of Veff to the minimization routine
	auto VeffValue = [&](ColumnVector vevs) {
		auto v = vevs(0);
		auto x = vevs(1); 
		//return real( EvaluatePotential(v, x, loopOrder, bDoDim6) );
		return real( EvaluatePotentialAsDouble(v, x, loopOrder, bDoDim6) );
	};

	// ---- BOBYQA minimization from dlib

	double singletLowerBound = -1e10;
	if (bIsZ2Symmetric) 
		singletLowerBound = 0.0;

	//ColumnVector lowerBound{5e-5, singletLowerBound};
	ColumnVector lowerBound{0, singletLowerBound};
	ColumnVector upperBound{1e10, 1e10};
	ColumnVector minimum = {v0, x0};
	try { 
		dlib::find_min_bobyqa(VeffValue, 
			minimum, 
			6,    // number of interpolation points
			lowerBound,  // lower bound constraint
			upperBound,  // upper bound constraint
			minParams.initialTrustRadius,    // initial trust region radius
			minParams.stoppingTrustRadius,  // stopping trust region radius
			minParams.maxFunctionEvaluations    // max number of objective function evaluations
		);
	} catch (dlib::bobyqa_failure& exc) {
		// Only print the exceptions in debug mode, otherwise there's too much spam. Warning flag is always toggled though 
		DEBUG("!!! Exception thrown in FindLocalMinimum():");
		DEBUG(exc.what());
		DEBUG("Reached point (v, x) = (" << minimum(0) << ", " << minimum(1) << ")\n\n");
		warnings++;
	}
	// Result goes into 'minimum'
	ParameterMap res;
	double v = minimum(0);
	double x = minimum(1);
	res["v"] = v;
	res["x"] = x;

	std::complex<double> val = this->EvaluatePotentialAsDouble(v, x, loopOrder, bDoDim6);
	res["Veff.re"] = real(val);
	res["Veff.im"] = imag(val);
	
	return res;
}



/* Global minimization. Looks for several local minima (based on intuitive guesses) and takes the deepest of those. 
Works with doubles, which get converted to Complex for internal computations */
template <typename Complex>
ParameterMap EffPot<Complex>::FindGlobalMinimum(const ELoopOrderVeff loopOrder, bool bDoDim6) {


	std::vector<std::array<double, 2>> startingPoints = InitialSearchPoints();

	// Our current result for the global minimum
	ParameterMap globalMinimum;

	MinimizationParams minParams;

	for (unsigned int i=0; i<startingPoints.size(); i++) {

		double v = startingPoints[i][0];
		double x = startingPoints[i][1];

		// singlet field can take large values if b4 is very small
		minParams.initialTrustRadius = std::max(10.0, std::abs(0.6*x)); 

		ParameterMap minimum = this->FindLocalMinimum(loopOrder, bDoDim6, v, x, minParams);

/* 		std::cout << "\nMinimum at" << "\n";
		PrintMap(minimum);
		std::cout << "\n"; */

		// Check if the present minimum is deeper than what we found earlier
		double val = GetFromMap(minimum, "Veff.re");
		if (i == 0 || val < GetFromMap(globalMinimum, "Veff.re")) {
			globalMinimum = minimum;
		}
	}

	// Potential should be real in its minima. Print warning if there is a relatively large imag part
	double re = GetFromMap(globalMinimum, "Veff.re");
	double im = GetFromMap(globalMinimum, "Veff.im");
	if (abs(im/re) > 1e-5) {
		// std::cout << "! Imaginary part in free energy\n";
		warnings++;
	}

	return globalMinimum;
}

template <typename Complex>
ParameterMap EffPot<Complex>::GetMassEigenvalues() {

	ParameterMap masses;
	masses["mh1sq"] = (double)mh1sq;
	masses["mh2sq"] = (double)mh2sq;
	masses["mGsq"] = (double)mGsq;
	masses["mGpmsq"] = (double)mGpmsq;
	masses["mWsq"] = (double)mWsq;
	masses["mZsq"] = (double)mZsq;
	return masses;
}



/* Calculate relative shifts to the VEVs due to dim-6 operators (simple tree-level estimate).
The logic here is not fully airtight: The analytical formulae used by the function assume calculate relative shifts 
to the TREE LEVEL minimum (v0, x0), yet here I plug in the minimum of the loop-corrected potential. So what this calculation
really gives is the relative shifts to the minimum location at this temperature, IF the tree level minimum would be at 
(v0, x0) = (v, x). So we're calculating the error for hypothetical, hand-picked field values like in the paper hep-ph/9508379.
*/
// UPDATE 12-04-2023: This is not used in the actual scans, see main.cpp for alternative (numerical minimization with dim-6). 
template <typename Complex>
std::vector<double> EffPot<Complex>::FieldShiftsDim6(Float v, Float x) {
	Float vRelativeShift = 0.;
	Float xRelativeShift = 0.;
	// Dim-6 error estimate (analytical expression based on tree-level minimization)
	// Taken from Mathematica
	Float msq = msqPhi;

	if (abs(v) > 10.*vWarningThreshold) {
		vRelativeShift = ( (-(v*(a1 + 2*a2*x)*(c41*(v*v*v*v) + 
			2*x*(c42*(v*v*v*v) + 
			x*(3*c23*(v*v) + 4*c24*(v*v)*x + 10*c05*(x*x) + 
			12*c06*(x*x*x))))) + 
			(2*b2 + a2*(v*v) + 4*b3*x + 6*b4*(x*x))*
			(3*c60*(v*v*v*v*v) + 
			4*v*x*(c41*(v*v) + x*(c42*(v*v) + x*(c23 + c24*x)))))/
			(8.*v*(((a1*v)/2. + a2*v*x)*((a1*v)/2. + a2*v*x) - 
			((2*msq + 6*lambda*(v*v) + a1*x + a2*(x*x))*
			(2*b2 + a2*(v*v) + 4*b3*x + 6*b4*(x*x)))/4.))
		);
	}
	if (abs(x) > 10.*vWarningThreshold) {
		xRelativeShift = ( ((2*msq + 6*lambda*(v*v) + a1*x + a2*(x*x))*
			(c41*(v*v*v*v) + 2*x*
			(c42*(v*v*v*v) + 
			x*(3*c23*(v*v) + 4*c24*(v*v)*x + 10*c05*(x*x) + 
			12*c06*(x*x*x)))) - 
			v*(a1 + 2*a2*x)*(3*c60*(v*v*v*v*v) + 
			4*v*x*(c41*(v*v) + x*(c42*(v*v) + x*(c23 + c24*x)))))/
			(8.*x*(((a1*v)/2. + a2*v*x)*((a1*v)/2. + a2*v*x) - 
			((2*msq + 6*lambda*(v*v) + a1*x + a2*(x*x))*
			(2*b2 + a2*(v*v) + 4*b3*x + 6*b4*(x*x)))/4.))
		);
	}
	std::vector<double> res{(double) vRelativeShift, (double) xRelativeShift};
	return res;
}



template <typename Complex>
void EffPot<Complex>::TestSmallLimit(std::string outfile) {

	using std::ofstream;
	ofstream file(outfile);
	if (file.is_open()) {

		file.precision(24);
		typename EffPot<Complex>::Float v = 0;

		// For reference, evaluate exactly at v=0
		Float veff0 = EvaluatePotential(0, this->x, ELoopOrderVeff::loop2, false).real();

		for (v=1e-16; v < 1e-3; v *= 1.25) {
			Float val = EvaluatePotential(v, this->x, ELoopOrderVeff::loop2, false).real();
			file << v << " " << abs((val - veff0)/veff0) << "\n";
		}
		
		file.close();
	}
  	else std::cout << "! Unable to open output file for TestSmallLimit(). (in effpot.tpp)\n";

}