#include "renormalization.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

/*
// Old implementation that I used with Boost odeint
void BetaFunctions(const std::vector<double> &x, std::vector<double> &dxdmu, double scale) {
	
	// This routine calculates dx / dmu = ... 
	// where x = coupling and mu = MS-bar scale. Remember that beta function is mu * dx / dmu

	// Common prefactor: dx/dmu = 1/mu 1/(8pi^2) * (...)
	const double pf = 1.0 / (8.0*PI*PI * scale);
	// Number of fermion __generations__
	const int Nf = 3;

	double g1sq = x[0];
	double g2sq = x[1];
	double g3sq = x[2];
	double ytsq = x[3];
	double msqPhi = x[4];
	double lambda = x[5];
	//double b1 = x[6]; // this does not appear in the beta functions
	double b2 = x[7];
	double b3 = x[8];
	double b4 = x[9];
	double a1 = x[10];
	double a2 = x[11];
	
	// U(1) coupling squared
	dxdmu[0] = pf * g1sq*g1sq * (1./6. + 20.*Nf / 9.);
	// SU(2) coupling squared
	dxdmu[1] = -pf * g2sq*g2sq * (43./6. - 4.*Nf / 3.);
	// SU(3) coupling squared. Does not run in our approximation!  
	dxdmu[2] = 0.;
	// Top yukawa coupling squared
	dxdmu[3] = pf * (9./2. * ytsq*ytsq - 8*g3sq*ytsq - 9./4. * g2sq*ytsq - 17./12. * g1sq*ytsq);
	// Higgs mass parameter m^2
	dxdmu[4] = pf * (msqPhi * (-9./4.*g2sq - 3./4.*g1sq + 6*lambda + 3*ytsq) + a1*a1 / 4. + 0.5 * a2*b2);
	// Higgs lambda
	dxdmu[5] = pf * (12*lambda*lambda + 9./16.*g2sq*g2sq + 3./8.*g2sq*g1sq + 3./16.*g1sq*g1sq - 3*ytsq*ytsq
				- 3./2.*lambda * (3*g2sq + g1sq) + 6*lambda*ytsq + 1./4. * a2*a2
			);
	// Singlet b1
	dxdmu[6] = pf * (b3*b2 + msqPhi*a1);
	// Singlet b2
	dxdmu[7] = pf * (2*b3*b3 + 0.5*a1*a1 + 3*b4*b2 + 2*a2*msqPhi);
	// Singlet b3
	dxdmu[8] = 3.*pf * (3*b3*b4 + 0.5*a1*a2);
	// Singlet b4
	dxdmu[9] = pf * (9*b4*b4 + a2*a2);
	// a1
	dxdmu[10] = pf * (2*a2*b3 + a1 * (-9./4.*g2sq - 3./4.*g1sq + 6*lambda + 2*a2 + 3*ytsq));
	// a2
	dxdmu[11] = pf * a2 * (-9./4.*g2sq - 3./4.*g1sq + 3*ytsq + 2*a2 + 6*lambda + 3*b4);

	// Derivatives done.
}
*/

namespace gslBetaFunctions {

	int BetaFunctions(double scale, const double x[], double dxdmu[], void *params) {
		// suppress unused parameter warnings
		(void)params;

		// This routine calculates dx / dmu = ... 
		// where x = coupling and mu = MS-bar scale. Remember that beta function is mu * dx / dmu

		// Common prefactor: dx/dmu = 1/mu 1/(8pi^2) * (...)
		const double pf = 1.0 / (8.0*PI*PI * scale);
		// Number of fermion __generations__
		const int Nf = 3;

		double g1sq = x[0];
		double g2sq = x[1];
		double g3sq = x[2];
		double ytsq = x[3];
		double msqPhi = x[4];
		double lambda = x[5];
		//double b1 = x[6]; // this does not appear in the beta functions
		double b2 = x[7];
		double b3 = x[8];
		double b4 = x[9];
		double a1 = x[10];
		double a2 = x[11];
		
		// U(1) coupling squared
		dxdmu[0] = pf * g1sq*g1sq * (1./6. + 20.*Nf / 9.);
		// SU(2) coupling squared
		dxdmu[1] = -pf * g2sq*g2sq * (43./6. - 4.*Nf / 3.);
		// SU(3) coupling squared. Does not run in our approximation!  
		dxdmu[2] = 0.;
		// Top yukawa coupling squared
		dxdmu[3] = pf * (9./2. * ytsq*ytsq - 8*g3sq*ytsq - 9./4. * g2sq*ytsq - 17./12. * g1sq*ytsq);
		// Higgs mass parameter m^2
		dxdmu[4] = pf * (msqPhi * (-9./4.*g2sq - 3./4.*g1sq + 6*lambda + 3*ytsq) + a1*a1 / 4. + 0.5 * a2*b2);
		// Higgs lambda
		dxdmu[5] = pf * (12*lambda*lambda + 9./16.*g2sq*g2sq + 3./8.*g2sq*g1sq + 3./16.*g1sq*g1sq - 3*ytsq*ytsq
					- 3./2.*lambda * (3*g2sq + g1sq) + 6*lambda*ytsq + 1./4. * a2*a2
				);
		// Singlet b1
		dxdmu[6] = pf * (b3*b2 + msqPhi*a1);
		// Singlet b2
		dxdmu[7] = pf * (2*b3*b3 + 0.5*a1*a1 + 3*b4*b2 + 2*a2*msqPhi);
		// Singlet b3
		dxdmu[8] = 3.*pf * (3*b3*b4 + 0.5*a1*a2);
		// Singlet b4
		dxdmu[9] = pf * (9*b4*b4 + a2*a2);
		// a1
		dxdmu[10] = pf * (2*a2*b3 + a1 * (-9./4.*g2sq - 3./4.*g1sq + 6*lambda + 2*a2 + 3*ytsq));
		// a2
		dxdmu[11] = pf * a2 * (-9./4.*g2sq - 3./4.*g1sq + 3*ytsq + 2*a2 + 6*lambda + 3*b4);

		return GSL_SUCCESS;
	}

}


ParameterMap Renormalization::RunToScale(double scaleOut, const ParameterMap &params) {

	// Integration range
	double scaleInit = GetFromMap(params, "RGScale");
	double scaleGoal = scaleOut;
	// Step size
	double delta = 0.05;

	// If we run from high to low scale, need to flip the sign of delta
	if (scaleGoal < scaleInit) 
		delta *= -1;

	// Put initial couplings into an array. Ordering needs to be the same as in BetaFunctions()
	const size_t NPARAMS = params.size() - 1; // don't count the RGScale
	std::vector<double> x(NPARAMS);
	x[0] = GetFromMap(params, "g1sq");
	x[1] = GetFromMap(params, "g2sq");
	x[2] = GetFromMap(params, "g3sq");
	x[3] = GetFromMap(params, "ytsq");
	x[4] = GetFromMap(params, "msqPhi");
	x[5] = GetFromMap(params, "lambda");
	x[6] = GetFromMap(params, "b1");
	x[7] = GetFromMap(params, "b2");
	x[8] = GetFromMap(params, "b3");
	x[9] = GetFromMap(params, "b4");
	x[10] = GetFromMap(params, "a1");
	x[11] = GetFromMap(params, "a2");

	// GSL ODE system and driver
	gsl_odeiv2_system sys = {gslBetaFunctions::BetaFunctions, nullptr, NPARAMS, nullptr};
	// last 3 arguments here are irrelevant when using fixed step algorithm
	gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, delta, 1e-6, 0.0);

	// Now integrate. For fixed-step evolution we need to manually calculate the number of steps
	// this many full steps
	const unsigned long nsteps = (scaleGoal - scaleInit) / delta; 
	// leftover delta for final step, so that we end up exactly at scaleGoal
	double deltaFinal = scaleGoal - (scaleInit + nsteps*delta);

	double scale = scaleInit;

	// now do the N full steps
	int status1 = gsl_odeiv2_driver_apply_fixed_step(driver, &scale, delta, nsteps, &x[0]);
	// and final adjustment step
	int status2 = gsl_odeiv2_driver_apply_fixed_step(driver, &scale, deltaFinal, 1, &x[0]);

	if (status1 != GSL_SUCCESS || status2 != GSL_SUCCESS) {
		std::cerr << "!! Error when integrating beta functions. Reached scale: " << scale << "\n";
    }

	// Put the resulting couplings back to a ParameterMap and return it 
	ParameterMap newParams;
	newParams["RGScale"] = scaleGoal;
	newParams["g1sq"] = x[0];
	newParams["g2sq"] = x[1];
	newParams["g3sq"] = x[2];
	newParams["ytsq"] = x[3];
	newParams["msqPhi"] = x[4];
	newParams["lambda"] = x[5];
	newParams["b1"] = x[6];
	newParams["b2"] = x[7];
	newParams["b3"] = x[8];
	newParams["b4"] = x[9];
	newParams["a1"] = x[10];
	newParams["a2"] = x[11];

	gsl_odeiv2_driver_free(driver);

	/* !!! Too much manual copying here. I should probably have a whole class for describing the SM+singlet model,
	* with a vector container for all parameters in some fixed ordering. Then have a member hashmap for accessing those by name.
	* Could then use the same class for all different EFTs etc? */

	return newParams;
}