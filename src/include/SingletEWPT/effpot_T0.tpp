
template <typename Complex>
Complex EffPotT0<Complex>::EvaluatePotential(Float v, Float x, const ELoopOrderVeff loopOrder)
{
	Complex Veff{0, 0};

	this->SetBackgroundFields(v, x);
	// Now mass eigenvalues have been calculated

	int order = static_cast<int>(loopOrder);

	Veff += this->V0();
	
	if (order > 0) {
		Veff += this->V1();
	}

	return Veff;
}

template <typename Complex>
void EffPotT0<Complex>::CalcMixingAngle() {

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

	/* NB! Equation B3 in 2103.07467 corresponds to mixing angle with the opposite sign
	of what I use here (choice in the paper matches that used in T=0 calculation). 
	Similarly the masses in eqs. (B1)-(B2) and the list of couplings after (B13) have 
	opposite sign for theta than what I have in this code. */

	Float v2 = v*v;
	Float x2 = x*x;

	Float A = msqPhi + 3*lambda*v2 + a1*x / 2. + a2*x2 / 2.;
	
	Float B = b2 + 3*b4*x2 + 2*b3*x + a2*v2/2.;

	Float C = a1*v/2. + a2*v*x;

	if (A - B == 0) {
		if (C > 0) theta =  PI/4.;
		else theta = -PI/4.;
	} else {
		
		theta = 0.5 * atan(2*C / (A - B));	
	}
}

template <typename Complex>
void EffPotT0<Complex>::CalcMassEigenvalues() 
{

	Float v2 = v*v;
	Float x2 = x*x;

	// Gauge bosons
	this->mWsq = 1.0/4.0 * g2sq * v2;
	this->mZsq = 1.0/4.0 * (g1sq + g2sq) * v2;

	// Top quark
	this->mtsq = 1./2. * ytsq * v2;

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
}

// need funny typename syntax here because the return type is defined within the template class
template <typename Complex>
typename EffPotT0<Complex>::Float EffPotT0<Complex>::V0() 
{
	Float v2 = v*v;
	Float v4 = v2*v2;
	Float x2 = x*x;
	Float x3 = x*x*x;
	Float x4 = x*x*x*x;

	Float res = 0.5*msqPhi*v2 + 0.25*lambda*v4 + 0.25*a1*v2*x + 0.25*a2 * v2*x2
				+ b1*x + 0.5*b2 * x2 + 1.0/3.0*b3 * x3 + 0.25*b4 *x4; 

	return res;
}

template <typename Complex>
Complex EffPotT0<Complex>::V1()
{

	Float d = 4;
	Complex res{0,0};


	// Top quark loop
	res += -4. * 3. * J4(mtsq);

	// Gauge + ghost loops
	res += 2.0*(d-1.0) * J4(mWsq) + (d-1.0) * J4(mZsq); 
	// Scalar loops
	res += J4(mh1sq) + J4(mh2sq) + J4(mGsq) + 2.0*J4(mGpmsq);

	return res;
} 


// Find a local minimum with initial guess (v0, x0). Returns doubles.
template <typename Complex>
ParameterMap EffPotT0<Complex>::FindLocalMinimum(const ELoopOrderVeff loopOrder, double v0, double x0) {

	using ColumnVector = dlib::matrix<double, 0, 1>;

	// Calculate real part of the potential in a lambda expression and pass it to dlib for minimization
	auto VeffValue = [&](ColumnVector vevs) {
		double v = vevs(0);
		double x = vevs(1); 
		return real( EvaluatePotentialAsDouble(v, x, loopOrder) );
	};

	// BOBYQA minimization from dlib
	ColumnVector lowerBound{1e-6, -1000};
	ColumnVector upperBound{1000, 1000};
	ColumnVector minimum = {v0, x0};
	try { 
		dlib::find_min_bobyqa(VeffValue, 
			minimum, 
			6,    // number of interpolation points
			lowerBound,  // lower bound constraint
			upperBound,  // upper bound constraint
			1,    // initial trust region radius
			1e-5,  // stopping trust region radius
			200    // max number of objective function evaluations
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

	std::complex<double> val = this->EvaluatePotentialAsDouble(v, x, loopOrder);
	res["Veff.re"] = real(val);
	res["Veff.im"] = imag(val);
	return res;
}



/* Global minimization. Looks for several local minima (based on intuitive guesses) and takes the deepest of those. 
Works with doubles, which get converted to Complex for internal computations */
template <typename Complex>
ParameterMap EffPotT0<Complex>::FindGlobalMinimum(const ELoopOrderVeff loopOrder) {

	using ColumnVector = std::vector<double>;

	int NPOINTS = 3;

	// Initial guesses for minima
	std::vector<ColumnVector> startingPoint(NPOINTS);
	startingPoint[0] = {1e-1, 0.0};
	startingPoint[1] = {300, 0.0};
	startingPoint[2] = {1e-1, -50.0};

	// Our current result for the global minimum
	ParameterMap globalMinimum;

	for (int i=0; i<NPOINTS; i++) {

		double v = startingPoint[i][0];
		double x = startingPoint[i][1];
		ParameterMap minimum = this->FindLocalMinimum(loopOrder, v, x);

		// Check if the present minimum is deeper than what we found earlier
		double val = GetFromMap(minimum, "Veff.re");
		if (i == 0 || val < GetFromMap(globalMinimum, "Veff.re")) {
			globalMinimum = minimum;
		}
	}

	return globalMinimum;
}

template <typename Complex>
ParameterMap EffPotT0<Complex>::GetMassEigenvalues() {

	ParameterMap masses;
	masses["mtsq"] = (double)mtsq;
	masses["mh1sq"] = (double)mh1sq;
	masses["mh2sq"] = (double)mh2sq;
	masses["mGsq"] = (double)mGsq;
	masses["mGpmsq"] = (double)mGpmsq;
	masses["mWsq"] = (double)mWsq;
	masses["mZsq"] = (double)mZsq;
	return masses;
}