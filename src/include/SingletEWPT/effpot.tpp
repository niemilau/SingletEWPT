
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
std::vector<double> EffPot<Complex>::TreeLevelMinima() const {

	// Need to solve dV/dv = 0, dV/dx = 0 so that d^2V/dv^2 > 0 and hessian > 0
	std::vector<double> res(1, 0.0);
	return res;
}


// Find a local minimum with initial guess (v0, x0). Returns doubles.
template <typename Complex>
ParameterMap EffPot<Complex>::FindLocalMinimum(const ELoopOrderVeff loopOrder, bool bDoDim6, double v0, double x0) {

	using ColumnVector = dlib::matrix<double, 0, 1>;

	// Calculate real part of the potential in a lambda expression and pass it to dlib for minimization
	auto VeffValue = [&](ColumnVector vevs) {
		auto v = vevs(0);
		auto x = vevs(1); 
		//return real( EvaluatePotential(v, x, loopOrder, bDoDim6) );
		return real( EvaluatePotentialAsDouble(v, x, loopOrder, bDoDim6) );
	};

	// BOBYQA minimization from dlib (TODO tune the parameters here)
	double singletLowerBound = -10000;
	if (bIsZ2Symmetric) 
		singletLowerBound = 0.0;

	ColumnVector lowerBound{5e-5, singletLowerBound};
	ColumnVector upperBound{10000, 10000};
	ColumnVector minimum = {v0, x0};
	try { 
		dlib::find_min_bobyqa(VeffValue, 
			minimum, 
			6,    // number of interpolation points
			lowerBound,  // lower bound constraint
			upperBound,  // upper bound constraint
			5,    // initial trust region radius
			1e-4,  // stopping trust region radius
			1000    // max number of objective function evaluations
		);
	} catch (dlib::bobyqa_failure& exc) {
		std::cout << "!!! Exception thrown in FindLocalMinimum():\n";
		std::cout << exc.what() << "\n";
		std::cout << "Reached point (v, x) = (" << minimum(0) << ", " << minimum(1) << ")\n\n";
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

	using ColumnVector = std::vector<double>;


	// Initial guesses for minima
	int NPOINTS = 5;
	std::vector<ColumnVector> startingPoint(NPOINTS);

/*
	double singletLowerBound = -500;
	if (bIsZ2Symmetric) {
		singletLowerBound = 0.0;
	}
	
	for (int n=0; n<NPOINTS; n++) {
		// TODO seed!
		double d1 = drand48() / ((double) RAND_MAX + 1);
		double d2 = drand48() / ((double) RAND_MAX + 1);
		startingPoint[n] = {1e-6 + d1 * (1000.0 - 1e-6), singletLowerBound + d2 * (1000 - singletLowerBound)};
	}
*/

	startingPoint[0] = {1e-1, 0.0};
	startingPoint[1] = {20, 0.0};
	// In Z2 symmetric model, search only singlet values >= 0
	if (bIsZ2Symmetric) {
		startingPoint[2] = {1e-1, 20};
		startingPoint[3] = {1e-2, 5};
	}
	else {
		
		// If b4 is very small, singlet VEV can get large so set the initial search point accordingly
		// This is relevant in small mass region (mh2 << 100 GeV)
		if (abs(this->b4) < 0.01) {
			startingPoint[2] = {10, -500.0};
			startingPoint[3] = {1e-2, 500};
		} else {
			startingPoint[2] = {10, -4.0};
			startingPoint[3] = {1e-2, 5.0};
		}

	}
	// One point at very large field values
	startingPoint[4] = {500, 500};


	// Our current result for the global minimum
	ParameterMap globalMinimum;

	for (int i=0; i<NPOINTS; i++) {

		double v = startingPoint[i][0];
		double x = startingPoint[i][1];
		ParameterMap minimum = this->FindLocalMinimum(loopOrder, bDoDim6, v, x);

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

	// DEBUG: print mass eigenvalues at the global minimum
	//this->SetBackgroundFields(GetFromMap(globalMinimum, "v"), GetFromMap(globalMinimum, "x"));
	//std::cout << scale << " v=" << v << " x=" << x << " theta=" << theta << " " << mh1sq << " " << mh2sq << " " << mGsq << " " << mGpmsq << " " << mWsq << "\n";
	//std::cout << "re = " << re << " im = " << im << ". Z2 symmetric? "<< bIsZ2Symmetric << "\n";

	//TestDiagrams(GetFromMap(globalMinimum, "v"), GetFromMap(globalMinimum, "x"));
	//std::cout << "\n";

/*
	double xmin = 2.43746;
	for (double vmin=0.02; vmin<0.04; vmin+=0.002) {
		this->SetBackgroundFields(vmin, xmin);
		std::cout << scale << " v=" << vmin << " x=" << xmin << " theta=" << theta << " " << mh1sq << " " << mh2sq << " " << mGsq << " " << mGpmsq << " " << mWsq << "\n";
		TestDiagrams(vmin, xmin);
		std::cout << "\n";
	}
*/
	
	/*
	// Try global minimization? 

	auto VeffValue = [&](dlib::matrix<double, 0, 1> vevs) {
		auto v = vevs(0);
		auto x = vevs(1); 
		return real( EvaluatePotentialAsDouble(v, x, loopOrder, bDoDim6) );
	};


	Float singletLowerBound = -1000;
	if (bIsZ2Symmetric) {
		singletLowerBound = 0.0;
	}

	auto result = dlib::find_min_global(VeffValue, 
                                  {1e-6, singletLowerBound}, // lower bounds
                                  {1000, 1000}, // upper bounds
                                  std::chrono::milliseconds(5000) // run this long
                                  );
	double v = result.x(0);
	double x = result.x(1);
	// Calculate Veff once more in the minimum to get possible imag part
	Complex Veff = EvaluatePotential(v, x, loopOrder, bDoDim6);
	double re = real(Veff);
	double im = imag(Veff);

	ParameterMap globalMinimum;
	globalMinimum["v"] = v;
	globalMinimum["x"] = x;
	globalMinimum["Veff.re"] = re;
	globalMinimum["Veff.im"] = im;
	*/

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