
#include "renormalization.h"


ParameterMap Renormalization::CalcMS(ELoopOrderMS loopOrder) {

	ParameterMap MSParams;

	if (loopOrder != ELoopOrderMS::tree && loopOrder != ELoopOrderMS::loop1) {
		Die("!!! Invalid loop order in CalcMS()\n", 100);
	}

	if (loopOrder == ELoopOrderMS::loop1) 
		CalculateSelfEnergies();

	MSParams["RGScale"] = scale;
	/* These are direct inputs in our scheme, no loop corrections */
	MSParams["a2"] = a2;
	MSParams["b3"] = b3;
	MSParams["b4"] = b4;
	// QCD coupling does not run in this approximation
	MSParams["g3sq"] = g3sq;

	// Then parameters that get loop corrections (if loopOrder > 0)
	MSParams["ytsq"] = CalcYtsq(loopOrder);
	MSParams["g2sq"] = Calcg2sq(loopOrder);
	MSParams["g1sq"] = Renormalization::Calcg1sq(loopOrder);
	MSParams["lambda"] = Renormalization::CalcLambda(loopOrder);
	MSParams["msqPhi"] = Renormalization::CalcMsqPhi(loopOrder);
	MSParams["a1"] = Renormalization::Calca1(loopOrder);
	MSParams["b2"] = Renormalization::Calcb2(loopOrder);
	MSParams["b1"] = Renormalization::Calcb1(loopOrder);

	return MSParams;
}

/** Implement conversion "physical input" -> MS-bar parameters. From appendix A of 2103.07467 **/

double Renormalization::CalcYtsq(ELoopOrderMS loopOrder) {

	double res = 0.0;
	if (loopOrder == ELoopOrderMS::tree) {
		// Use tree level relation
		res = g0sq * Mt*Mt / (2.0 * MW*MW);

	} else if (loopOrder == ELoopOrderMS::loop1) {
		res = g0sq * Mt*Mt / (2.0 * MW*MW) * (1.0 + gsqDelta - selfEnergyW/(MW*MW) + 2.0*selfEnergyTop);
	}

	return res;
}


double Renormalization::Calcg2sq(ELoopOrderMS loopOrder) {

	double res = 0.0;
	if (loopOrder == ELoopOrderMS::tree) {
		res = g0sq;

	} else if (loopOrder == ELoopOrderMS::loop1) {
		res = g0sq * ( 1.0 + gsqDelta);
	} 

	return res;
}

double Renormalization::Calcg1sq(ELoopOrderMS loopOrder) {

	double res = 0.0;
	if (loopOrder == ELoopOrderMS::tree) {
		// Use tree level relation
		res = g0sq * (MZ*MZ - MW*MW) / (MW*MW);

	} else if (loopOrder == ELoopOrderMS::loop1) {
		double MW2 = MW*MW;
		double MZ2 = MZ*MZ;
		res = g0sq / MW2 * ( (MZ2 - MW2) * (1.0 + gsqDelta) - MZ2*selfEnergyW / MW2 + selfEnergyZ );
	} 

	return res;
}

double Renormalization::CalcLambda(ELoopOrderMS loopOrder) {

	double res = 0.0;

	if (loopOrder == ELoopOrderMS::tree) {
		// Use tree level relation
		res = (g0sq / (16.0*MW*MW)) * ( Mh1*Mh1 + Mh2*Mh2 + (Mh1 - Mh2)*(Mh1 + Mh2) * cos(2*theta) );

	} else if (loopOrder == ELoopOrderMS::loop1) {
		res = (g0sq / (16.0*MW*MW)) * ( 
			( Mh1*Mh1 + Mh2*Mh2 + (Mh1 - Mh2)*(Mh1 + Mh2) * cos(2*theta) ) * ( 1.0 + gsqDelta - selfEnergyW / (MW*MW))
			+ 2.0*ct*ct*selfEnergyH1 + 2.0*st*st*selfEnergyH2
		); 
	}

	return res;
}


double Renormalization::CalcMsqPhi(ELoopOrderMS loopOrder) {


	double res = 0.0;
	if (loopOrder == ELoopOrderMS::tree) {
		// Use tree level relation
		res = (1.0/4.0) * (-Mh1*Mh1 - Mh2*Mh2 + (-Mh1*Mh1 + Mh2*Mh2) * cos(2*theta) );

	} else if (loopOrder == ELoopOrderMS::loop1) {
		res = -1.0/4.0 * (Mh1*Mh1 + Mh2*Mh2 + (Mh1*Mh1 - Mh2*Mh2)*cos(2.0*theta)
			+ 2.0*ct*ct*selfEnergyH1 + 2.0*st*st*selfEnergyH2
		);
	} 

	return res;
}


double Renormalization::Calca1(ELoopOrderMS loopOrder) {


	double res = 0.0;
	if (loopOrder == ELoopOrderMS::tree) {
		// Use tree level relation
		res = (sqrt(g0sq) / MW) * (-Mh1*Mh1 + Mh2*Mh2) * cos(theta)*sin(theta);

	} else if (loopOrder == ELoopOrderMS::loop1) {
		res = sqrt(g0sq)/MW * (Mh2*Mh2 - Mh1*Mh1)*ct*st 
			* ( 1.0 + 0.5*gsqDelta + (selfEnergyH2 - selfEnergyH1)/(Mh2*Mh2 - Mh1*Mh1) - 0.5*selfEnergyW/(MW*MW) );
	} 

	return res;
}

double Renormalization::Calcb2(ELoopOrderMS loopOrder) {

	double res = 0.0;

	if (loopOrder == ELoopOrderMS::tree) {
		// Use tree level relation
		res = 0.5 * (Mh1*Mh1 + Mh2*Mh2 + (Mh2*Mh2 - Mh1*Mh1)*cos(2*theta) 
				- (4 * MW*MW * a2 / g0sq)  
			);
	} else if (loopOrder == ELoopOrderMS::loop1) {
		res = 0.5 * (Mh1*Mh1 + Mh2*Mh2 + (Mh2*Mh2 - Mh1*Mh1)*cos(2*theta) 
						+ 2.0*selfEnergyH1*st*st + 2.0*selfEnergyH2*ct*ct)
			- 2.0*MW*MW*a2/g0sq * (1.0 - gsqDelta + selfEnergyW / (MW*MW));
	} 

	return res;
}

double Renormalization::Calcb1(ELoopOrderMS loopOrder) {

	double res = 0.0;
	
	if (loopOrder == ELoopOrderMS::tree) {
		// Use tree level relation
		res = (MW / sqrt(g0sq)) * (Mh1*Mh1 - Mh2*Mh2) * cos(theta)*sin(theta);
	} else if (loopOrder == ELoopOrderMS::loop1) {
		res = -MW/sqrt(g0sq) * (Mh2*Mh2 - Mh1*Mh1)*ct*st 
			* ( 1.0 - 0.5*gsqDelta + (selfEnergyH2 - selfEnergyH1)/(Mh2*Mh2 - Mh1*Mh1) + 0.5*selfEnergyW/(MW*MW) );
	} 

	return res;
}

bool Renormalization::CheckPerturbativity(const ParameterMap &params) {

	/* Naive perturbativity bounds, from 0909.0520, below eq (25), and converted to my conventions */
	double lambda = GetFromMap(params, "lambda");
	double a2 = GetFromMap(params, "a2");
	double b4 = GetFromMap(params, "b4");

	return ( (lambda < 2.*PI / 3.) && (a2 < 4.*PI) && (b4 < 2.*PI / 3.) );
}


void Renormalization::CalculateSelfEnergies() {

	// Evaluate self energies at external momentum p^2 = M^2 (the corresponding pole mass)
	// See selfEnergy.cpp for the actual formulas (long!)

	// W boson
	this->selfEnergyW = CalcSelfEnergyW(MW);
	// Z boson
	this->selfEnergyZ = CalcSelfEnergyZ(MZ);
	// H1 boson
	this->selfEnergyH1 = CalcSelfEnergyH1(Mh1);
	// H2 boson
	this->selfEnergyH2 = CalcSelfEnergyH2(Mh2);
	// Top quark, scalar + vector parts
	this->selfEnergyTop = CalcSelfEnergyTop(Mt);

	// Correction to gauge coupling (match muon decay)
	this->gsqDelta = CalcgsqCorrection();

	/*
	std::cout.precision(16);
	std::cout << "Pi_W : " << selfEnergyW << "\n";
	std::cout << "Pi_z : " << selfEnergyZ << "\n";
	std::cout << "Pi_H1 : " << selfEnergyH1 << "\n";
	std::cout << "Pi_H2 : " << selfEnergyH2 << "\n";
	std::cout << "Pi_t : " << selfEnergyTop << "\n";
	std::cout << "Delta gsq : " << gsqDelta << "\n";
	*/
}

