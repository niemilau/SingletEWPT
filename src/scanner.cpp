#include "scanner.h"

std::vector<double> Scanner::MakeLinearGrid(double min, double max, double delta) {
	std::vector<double> values;
	values.reserve(std::abs((max-min) / delta));

	// Using some trickery for upper bound here to dodge rounding errors
	for (double x = min; x <= max + 0.01*delta; x += delta) {
		values.push_back(x);
	}
	return values;
}

std::vector<double> Scanner::ReadNumbersFromFile(std::string fname) {
	Die("TODO!!", 42131);
    return std::vector<double>();
}

void Scanner::ReadScannerParams(std::string fname)
{

    std::fstream paramFile;
	paramFile.open(fname);

	if (!paramFile.is_open()) {
		Die("!!! Cannot read, paramFile is not open...\n", 999);
	}

	std::string line;

	/* Need to interpret key,value pairs. I do this in two steps to make checks on input as automatic as possible:
	1) Read all numbers and strings into separate maps
	2) Using the safe-access method GetFromMap(), separate 'option flags' from 'scanning range' parameters */
	ParameterMap nums;
	std::map<std::string, std::string> strings;

	while( std::getline(paramFile, line) ) {

		// Skip comments
		if (line.empty() || line[0] == '#') continue;

		std::istringstream is_line(line);
        std::string key;
		std::string value;

		// Get each word in the string, this accounts for both tab and space delimiters
		// For now, assume only 2 words per line. TODO make this able to read eg. T		60 180 5
		is_line >> key;
		is_line >> value;
		// Step 1)
		if (IsNumber(value)) {
			nums[key] = std::stod(value);
		} else {
			strings[key] = value;
		}
	}

	paramFile.close();

	// Now step 2) with hard-coded key names

	// DR loop order needs manual interpreting (should do this more elegantly with a map, or by turning the enum into a proper class)
	{
		std::string orderAsString = GetFromMap(strings, "loopOrderDR");
		if (orderAsString == "LO") 
			loopOrderDR = ELoopOrderDR::LO;
		else if (orderAsString == "NLO")
			loopOrderDR = ELoopOrderDR::NLO;
		else if (orderAsString == "NLONo2Loop")
			loopOrderDR = ELoopOrderDR::NLONo2Loop;
		else 
			Die("!!! Invalid loopOrderDR: Choose from 'LO', 'NLO', 'NLONo2Loop'.\n", 62);
	}

	// Other enums, convert double -> int -> enum
	loopOrderMS = static_cast<ELoopOrder>( (int) GetFromMap(nums, "loopOrderMS") );
	nums.erase("loopOrderMS");
	loopOrderVeff = static_cast<ELoopOrder>( (int) GetFromMap(nums, "loopOrderVeff") );
	nums.erase("loopOrderVeff");
	loopOrderVeffT0 = static_cast<ELoopOrder>( (int) GetFromMap(nums, "loopOrderVeffT0") );
	nums.erase("loopOrderVeffT0");

	// Boolean flags
	bCalculateCondensates = !!GetFromMap(nums, "calculateCondensates");
	nums.erase("calculateCondensates");
	bStopAtSymmetricPhase = !!GetFromMap(nums, "stopAtSymmetricPhase");
	nums.erase("stopAtSymmetricPhase");
	bWriteAtEachTemperature = !!GetFromMap(nums, "writeAtEachTemperature");
	nums.erase("writeAtEachTemperature");
	bSolveBetas = !!GetFromMap(nums, "solveBetas");
	nums.erase("solveBetas");
	bOnlySearchEWPT = !!GetFromMap(nums, "onlySearchEWPT");
	nums.erase("onlySearchEWPT");

	// Miscellaneous
	jumpThreshold = GetFromMap(nums, "jumpThreshold");
	nums.erase("jumpThreshold");
	dim6TemperatureFraction = GetFromMap(nums, "dim6TemperatureFraction");
	nums.erase("dim6TemperatureFraction");

	matchingScaleOverride = GetFromMap(nums, "matchingScaleOverride");
	nums.erase("matchingScaleOverride");
	scale3DOverride = GetFromMap(nums, "scale3DOverride");
	nums.erase("scale3DOverride");

	if (matchingScaleOverride > 0.0) {
		bUseDefaultMatchingScale = false;
	}
	if (scale3DOverride > 0.0) {
		bUseDefault3DScale = false;
	}

	// This leaves only the parameter scanning ranges
	std::vector<double> range;

	auto SetParameterRange = [&](const std::string &name) {
		// Try to read the values from file first (name eg. range_a2)
		std::string fileName = "range_" + name;
		std::vector<double> values;
		if (FileExists(fileName)) {
			std::cout << "Reading " << name << " range from file " << fileName << "\n";
			values = ReadNumbersFromFile(fileName);
		} else {
			// Just use a linear grid based on number given in parameters file
			values = MakeLinearGrid(GetFromMap(nums, name+"_min"), GetFromMap(nums, name+"_max"), GetFromMap(nums, name+"_delta"));
		}
		this->scanningRange.insert( { name, values });

	};

	SetParameterRange("mh2");
	SetParameterRange("a2");
	SetParameterRange("b3");
	SetParameterRange("b4");
	SetParameterRange("sinTheta");
	SetParameterRange("T");

}

void Scanner::FindTransitionPoints() {

	// Safety checking: if minimization is not accurate, we may jump back and forth near a transition
    bool bPreviousTemperatureHasJump = false;

	std::vector<double> previousFields(2);

	for (int i=0; i < static_cast<int>(resultsForT.size()); i++) {
		std::vector<double> fields(2);
		fields[0] = GetFromMap(resultsForT[i], "vByT");
		fields[1] = GetFromMap(resultsForT[i], "xByT");

		// Compare field/T values at current and previous temperatures to see if there was a transition 
		if (i==0) {
			previousFields = fields;
			continue;
		} 

		// field/T jumps: Usually previous > current if going from low-to-high T 
		double vJump = previousFields[0] - fields[0];
		double xJump = previousFields[1] - fields[1];

		bool bFoundTransition = false;
		// Now check if there is a phase transition; the condition for this depends on whether we only search for 
		// "genuine" EWPT where Higgs jumps from v = 0 to nonzero (or vice versa), OR if we just base the check on
		// how much v/T or x/T changed in one temperature step. The latter method needs a minimum threshold value
		// and may result in a lot of messy data at low temperatures, where v/T may change rapidly even without a transition.
		// A warning flag will be raised in this case though

		// TODO make this less messy: fields[0] is v/T and fields[1] x/T. I know I'll get these wrong if I need to add more stuff here

		if (this->bOnlySearchEWPT) {
			// v/T smaller than this is regarded as being in symmetric phase v=0 
			double symmPhaseThreshold = 1e-4; 
			if ( (abs(previousFields[0]) < symmPhaseThreshold && abs(fields[0]) > symmPhaseThreshold)
				|| (abs(previousFields[0]) > symmPhaseThreshold && abs(fields[0]) < symmPhaseThreshold) ) {
					bFoundTransition = true;
				}
		} else {
			// Just check if any field had a large enough jump
			if (abs(vJump) >= jumpThreshold || abs(xJump) >= jumpThreshold) {
				bFoundTransition = true;
			}
		}

		if (bFoundTransition) {
			// Found transition, now calculate discontinuities and latent heat
			ParameterMap transitionValues;

			// Keep track of possible issues. Add up warnings from nearby temperatures that are used for calculation here
			int warningsMinimization = 0;
			int warningsDerivatives = 0; // latent heat, condensates etc
			// Two transitions in a row?!
			if (bPreviousTemperatureHasJump) 
				warningsDerivatives++; // could be a minimization issue too tbh

			ParameterMap p1 = resultsForT[i-1];
			ParameterMap p2 = resultsForT[i];
			warningsMinimization += GetFromMap(p1, "warningsMinimization") + GetFromMap(p2, "warningsMinimization");
			warningsDerivatives += GetFromMap(p1, "warningsDerivatives") + GetFromMap(p2, "warningsDerivatives");

			// Perturbativity is considered lost if any nearby T has nonperturbative couplings
			bool bIsPerturbative = GetFromMap(p1, "isPerturbative") && GetFromMap(p2, "isPerturbative");

			// Critical temperature. Take from p1 (lower temperature if using low-to-high scanning)
			double Tc = GetFromMap(p1, "T");

			// Effective potential value
			double VeffRe = GetFromMap(p1, "Veff.re");
			double VeffIm = GetFromMap(p1, "Veff.im");

			// Higgs condensate discontinuity v_phys / T = sqrt( 2\Delta <phi^+phi> ) / T. Take abs because of the sqrt
			double vphysByT = GetFromMap(p1, "phisqByT2") - GetFromMap(p2, "phisqByT2");
			vphysByT = sqrt(2. * abs(vphysByT));

			double LByT4 = 0.0;

			// Latent heat. Need dV/dT on both sides of the transition
			if (i-2 < 0 || i+1 >= static_cast<int>(resultsForT.size())) {
				warningsDerivatives++;
				DEBUG("!!! Not enough data points for calculating latent heat at Tc = " << Tc);
			} else {
				// Low-T derivative
				ParameterMap pOther = resultsForT[i-2];
				double T1 = GetFromMap(pOther, "T");
				double T2 = GetFromMap(p1, "T");
				double val1 = GetFromMap(pOther, "Veff.re");
				double val2 = GetFromMap(p1, "Veff.re");
				double dVdT1 = (val2 - val1) / (T2 - T1);
				warningsMinimization += GetFromMap(pOther, "warningsMinimization");
				warningsDerivatives += GetFromMap(pOther, "warningsDerivatives");
				bIsPerturbative *= GetFromMap(pOther, "isPerturbative");

				// High-T derivative
				pOther = resultsForT[i+1];
				T1 = GetFromMap(pOther, "T");
				T2 = GetFromMap(p2, "T");
				val1 = GetFromMap(pOther, "Veff.re");
				val2 = GetFromMap(p2, "Veff.re");
				double dVdT2 = (val2 - val1) / (T2 - T1);
				warningsMinimization += GetFromMap(pOther, "warningsMinimization");
				warningsDerivatives += GetFromMap(pOther, "warningsDerivatives");
				bIsPerturbative *= GetFromMap(pOther, "isPerturbative");
				
				LByT4 = abs(dVdT2 - dVdT1) / (Tc*Tc);
			}


			// Calculate dim-6 error at LOWER temperature than Tc. This is necessary because addition of 
			// dim-6 operators can shift Tc so much that eg. the minimum changes completely. 
			// Also, we ultimately need to know the error at T < Tc because of supercooling 

			// Lower T by some percents
			double dim6T = dim6TemperatureFraction * Tc;
			double dT = GetFromMap(p2, "T") - Tc;
			// Index of T = dim6T 
			int ind_Tc = i-1;
			int ind_new = ind_Tc - (Tc - dim6T) / dT;
			
			double vShiftDim6 = 0.0; // GetFromMap(p1, "vShiftDim6");
			double xShiftDim6 = 0.0; // GetFromMap(p1, "xShiftDim6");

			if (ind_new < 0 || ind_new >= static_cast<int>(resultsForT.size())) {
				// index out of bounds, can't calculate dim-6 error
				warningsDerivatives++;
				std::cout << "! Can't calculate dim-6 error; Tc was " << Tc << "\n";
			} else {
				ParameterMap params_for_dim6 = resultsForT[ind_new];
				vShiftDim6 = GetFromMap(params_for_dim6, "vShiftDim6");
				xShiftDim6 = GetFromMap(params_for_dim6, "xShiftDim6");
				//std::cout << "Tc = " << Tc << ", dim-6 error at T = " << GetFromMap(params_for_dim6, "T") << "\n";
			}

			// Put these in a vector and write to file
			const int NCOLUMNS = 18;
			std::vector<double> res(NCOLUMNS);
			res[0] = GetFromMap(currentInput, "Mh1");
			res[1] = GetFromMap(currentInput, "Mh2");
			res[2] = GetFromMap(currentInput, "a2");
			res[3] = GetFromMap(currentInput, "b4");
			res[4] = GetFromMap(currentInput, "sinTheta");
			res[5] = GetFromMap(currentInput, "b3");
			res[6] = Tc;
			res[7] = vphysByT;
			res[8] = vJump;
			res[9] = xJump;
			res[10] = LByT4;
			// Dim-6 error estimate dv/v and dx/x at the transition point. Take always from low-T phase; this avoids issues with v = 0
			res[11] = vShiftDim6;
			res[12] = xShiftDim6;
			// Store Veff / T^4, but our Veff is in 3D units => actually Veff3D/T^3
			res[13] = VeffRe / (Tc*Tc*Tc);
			res[14] = VeffIm / (Tc*Tc*Tc);
			res[15] = (int) bIsPerturbative;
			res[16] = warningsMinimization;
			res[17] = warningsDerivatives;

			AppendToFile(this->transitionsFileName, res);

			bPreviousTemperatureHasJump = true;				
		} else {
			bPreviousTemperatureHasJump = false;	
		}
		
		previousFields = fields;
	}

}

void Scanner::WriteTemperatureData() {

	// Write results separately for each temperature. 
	// Use same ordering as in WriteDataLabels()

	const int NCOLUMNS = 18;
	std::vector<double> res(NCOLUMNS);

	for (ParameterMap map : resultsForT) {

		double T = GetFromMap(map, "T");
		res[0] = GetFromMap(currentInput, "Mh1");
		res[1] = GetFromMap(currentInput, "Mh2");
		res[2] = GetFromMap(currentInput, "a2");
		res[3] = GetFromMap(currentInput, "b4");
		res[4] = GetFromMap(currentInput, "sinTheta");
		res[5] = GetFromMap(currentInput, "b3");
		res[6] = GetFromMap(map, "T"); // this is now obviously T, not Tc
		res[7] = GetFromMap(map, "phisqByT2"); // NOT v/T, but square of that
		res[8] = GetFromMap(map, "vByT");
		res[9] = GetFromMap(map, "xByT");
		res[10] = 0.0; // dummy value for "latent heat" 
		res[11] = GetFromMap(map, "vShiftDim6");
		res[12] = GetFromMap(map, "xShiftDim6");
		// Make Veff values dimensionless like in the other output
		res[13] = GetFromMap(map, "Veff.re") / (T*T*T);
		res[14] = GetFromMap(map, "Veff.im") / (T*T*T);
		res[15] = GetFromMap(map, "isPerturbative");
		res[16] = GetFromMap(map, "warningsMinimization");
		res[17] = GetFromMap(map, "warningsDerivatives");
		AppendToFile(this->temperatureDataFileName, res);
	}

}

void Scanner::WriteDataLabels() {
	// Phase transition labels 
	std::ofstream f;
	f.open("labels");
	f << "1 mh1\n";
	f << "2 mh2\n";
	f << "3 a2\n";
	f << "4 b4\n";
	f << "5 sinTheta\n";
	f << "6 b3\n";
	f << "7 Tc\n";
	f << "8 sqrt(2 Delta <phi^+phi>) / Tc\n";
	f << "9 v/Tc jump\n";
	f << "10 x/Tc jump\n";
	f << "11 L/Tc^4\n"; 
	f << "12 Dim-6 error estimate: delta v / v (abs value)\n";
	f << "13 Dim-6 error estimate: delta x / x (abs value)\n";
	f << "14 Re Veff / Tc^4\n";
	f << "15 Im Veff / Tc^4\n";
	f << "16 perturbative\n";
	f << "17 warnings (minimization)\n";
	f << "18 warnings (derivatives)\n";
	f.close();

	// T=0 labels
	f.open("labels_T0");
	f << "1 mh1\n";
	f << "2 mh2\n";
	f << "3 a2\n";
	f << "4 b4\n";
	f << "5 sinTheta\n";
	f << "6 b3\n";
	f << "7 stability\n";
	f << "8 perturbativity\n";
	f.close();

}


bool Scanner::CheckT0Stability(const ParameterMap &MSParams) {

	EffPotT0<std::complex<double>> effPotT0(MSParams);
	ParameterMap minimum = effPotT0.FindGlobalMinimum(this->loopOrderVeffT0);
	bool bGlobalMinimumIsEW = ( abs(GetFromMap(minimum, "v")) > 1 );
	bool bIsPerturbative = Renormalization::CheckPerturbativity(MSParams);

	// Write to T0 output file
	const int NCOLUMNS = 8;
	std::vector<double> res(NCOLUMNS);
	res[0] = GetFromMap(currentInput, "Mh1");
	res[1] = GetFromMap(currentInput, "Mh2");
	res[2] = GetFromMap(currentInput, "a2");
	res[3] = GetFromMap(currentInput, "b4");
	res[4] = GetFromMap(currentInput, "sinTheta");
	res[5] = GetFromMap(currentInput, "b3");
	res[6] = (int) bGlobalMinimumIsEW;
	res[7] = (int) bIsPerturbative;
	AppendToFile("data_T0.dat", res);
	return bGlobalMinimumIsEW;
}

void Scanner::AppendToFile(std::string fname, std::vector<double> data) {
	std::ofstream f;
	f.open(fname, std::ios_base::app);
	for (int i=0; i<static_cast<int>(data.size()); i++) {
		if (i > 0) 
			f << " "; 
		f << data[i];
	}
	f << "\n";
	f.close();
}


double Scanner::GetMatchingScale() const {
	if (bUseDefaultMatchingScale) {
		return 4.0*PI * currentTemperature * exp(-EULERGAMMA);
	} else {
		return matchingScaleOverride;
	}
 }

double Scanner::Get3DScale() const {
	if (bUseDefault3DScale) {
		return currentTemperature;
	} else {
		return scale3DOverride;
	}
 }