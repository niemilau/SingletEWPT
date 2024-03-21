
#include "common.h"
#include "scanner.h"
#include "renormalization.h"
#include "dimred.h"
#include "effpot.h"
#include "effpot_T0.h"

// Execution time evaluation
#include <chrono>
#include <utility>
typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define Duration(a) std::chrono::duration_cast<std::chrono::seconds>(a).count()
#define TimeNow() std::chrono::high_resolution_clock::now()

// Multiprecision is too slow to use in real scans... double is 15 digits and running a test scan over T (151 values) takes 0.120s. 
// With multiprecision 16 digits it takes 6.255s, which is 52 times more! So use standard double whenever possible.

//int main(int argc, char *argv[]) {
int main() {

	using Complex = std::complex<double>;

	/* Read parameter ranges and other options */
	Scanner scanner("parameters");

	std::string statusFileName = "status";

    std::cout << "====== Scanner options =====\n";
	scanner.PrintScanner();
	std::cout << "============================\n\n"; 
	std::cout << std::flush;

	scanner.WriteDataLabels();

	// Other options (TODO move these inside the Scanner class)
	bool bDoDim6 = true; 
	bool bNLOCubics = false;

	TimeVar startTime = TimeNow();

	/****** Parameter loops ******/

	// Keep track of how many parameter points we've visited. Used for evaluation control
	long pointCount = 0;
	long checkpointInterval = 10000;

	for (double const &mh2 : GetFromMap(scanner.scanningRange, "mh2") ) 
	for (double const &sinTheta : GetFromMap(scanner.scanningRange, "sinTheta") ) 
	for (double const &a2 : GetFromMap(scanner.scanningRange, "a2") )  
	for (double const &b3 : GetFromMap(scanner.scanningRange, "b3") ) 
	for (double const &b4 : GetFromMap(scanner.scanningRange, "b4") ) 
	{
		scanner.currentInput = {
			{"Mh1", ExperimentalInput::MH},
			{"Mh2", mh2},
			{"a2", a2},
			{"sinTheta", sinTheta},
			{"b3", b3},
			{"b4", b4}
		};

		pointCount++;
		if (pointCount % checkpointInterval == 0) {
			double seconds = Duration(TimeNow() - startTime);
			std::cout << pointCount << " points done, time taken so far: " << seconds << "s. Current parameters:\n";
			PrintMap(scanner.currentInput);
			std::cout << std::endl;
		}

		// Calculate MS-bar parameters at Z pole
		double inputScale = ExperimentalInput::MZ;
		Renormalization renorm(scanner.currentInput, inputScale);

		ParameterMap MSParams = renorm.CalcMS(scanner.loopOrderMS);

		// Store these for optimizing RG running inside the T-loop 
		ParameterMap paramsAtPreviousT = MSParams;

		/** Minimization of T=0 potential, ie. check if v != 0 minimum is the global one at our starting scale **/
		bool bGlobalMinimumIsEW = scanner.CheckT0Stability(MSParams);
		if (!bGlobalMinimumIsEW) {
			// Don't look for phase transitions if already the T=0 vacuum is wrong
			continue;
		}

		/******** Now the T-loop ********/
		scanner.StartTemperatureLoop();

		/* if bStopAtSymmetricPhase == true, the scanner will stop once v becomes small in the global minimum.
		But it can happen that at small temperature our high-T approx breaks, and the minimum shifts to v=0.
		This is clearly not physical since we assumed v != 0 at T = 0. To avoid stopping the scanner,
		introduce a check that we've been in broken phase at low T before applying the bStopAtSymmetricPhase check. */
		bool bFoundBrokenPhase = false;

		// For optimization at v = 0
		int symmPhasePointsMax = 4;
		int symmPhasePoints = symmPhasePointsMax;
		bool bInSymmetricPhase = false;

		for (double T : GetFromMap(scanner.scanningRange, "T") ) 
		{

			// RG running will not work if you forget to set this here!
			scanner.SetTemperature(T);

			// Run parameters to matching scale by integrating beta functions. 
			// If the matching scale is T-dependent, can optimize this by getting the params at previous T and running from there.
			// This is faster than repeatedly integrating eg. from MZ to 7T at least if the T-spacing is small and we are not doing binary search
			ParameterMap paramsForDR;
			if (scanner.bSolveBetas) {
				paramsForDR = Renormalization::RunToScale(scanner.GetMatchingScale(), paramsAtPreviousT);
			} else {
				paramsForDR = MSParams;
			}
			paramsAtPreviousT = paramsForDR;

			bool bIsPerturbative = Renormalization::CheckPerturbativity(paramsForDR);

			// Dimensional reduction to 3D SM + singlet (incl. running to 3D scale)
			double scale3D = scanner.Get3DScale();
			ParameterMap params3D = DimRed::DoFullDimRed(T, paramsForDR, scanner.loopOrderDR, scale3D, bDoDim6, bNLOCubics);

			// Construct and minimize effective potential. This operates in 3D units: [v] = [x] = 1/2, [Veff] = 3. Convert these to 4D units below
			EffPot<Complex> effPot(params3D);
			ParameterMap minimum = effPot.FindGlobalMinimum(scanner.loopOrderVeff, false);

			// Keep track of possible issues. Separate these into "minimization", ie those reported by eff. pot. evaluation/minimization,
			// and "derivatives", ie issues encountered when computing latent heat or condensates etc
			int warningsMinimization = effPot.warnings;
			int warningsDerivatives = 0;

			double v = GetFromMap(minimum, "v");
			double x = GetFromMap(minimum, "x");
			// Convert minimum location to 4D units and normalize by T (=> v/T)
			double vByT = v / sqrt(T);
			double xByT = x / sqrt(T);


			// Relative VEV shifts due to dim-6 operators 
			std::vector<double> fieldShiftsDim6{0.0, 0.0};

			// OLD: don't use this. See below
			//fieldShiftsDim6 = effPot.FieldShiftsDim6(v, x);


			// Calculate quadratic Higgs condensate <phi^+ phi> = dV/dmsq. This is relatively expensive because we need to minimize Veff again 
			double phisq = 0.0;
			if (scanner.bCalculateCondensates) {
				ParameterMap newParams3D(params3D); // copy-constructor
				double msqPhi = GetFromMap(params3D, "msqPhi");
				// Some care required here: need to change msqPhi, but we don't want the system to end up in a different minimum.
				// Can help this by choosing smaller msq in the broken phase, and larger msq in the symmetric phase
				double diff = 0.015*msqPhi; // 0.001 seems too small here. at 2 loop it gives messy results, 1loop seems fine. TODO test carefully
				if (abs(v) > 1e-4) 
					diff *= -1.;
					
				double msqPhi_new = msqPhi + diff;
				newParams3D["msqPhi"] = msqPhi_new;

				// Don't modify the original effPot here, would be too error prone
				EffPot<Complex> effPotNew(newParams3D);
				// Find LOCAL minimum near to the original one
				// Better results by reducing initial trust radius?
				MinimizationParams minParams;
				minParams.initialTrustRadius = 5;
				ParameterMap newMinimum = effPotNew.FindLocalMinimum(scanner.loopOrderVeff, false, v, x, minParams);

				// Sensibility check
				double sensitivity = 0.5;
				double vNew = GetFromMap(newMinimum, "v");
				double xNew = GetFromMap(newMinimum, "x");
				if ( abs( vNew / sqrt(T) - vByT) > sensitivity) 
				{
					warningsDerivatives++;
					(void)xNew;
					DEBUG("!!! Warning: Higgs condensate at T = " << T << ". New minimum is at (v, x) = (" << vNew << ", " << xNew 
								<< "), used to be (" << v << ", " << x << ")");
				}

				phisq = ( GetFromMap(newMinimum, "Veff.re") - GetFromMap(minimum, "Veff.re") ) / (msqPhi_new - msqPhi);	
			}

			// Minimize again with dim-6 included to find relative shift in v/T. Since we do this at all temperatures, 
			// can later in scanner.FindTransitionPoints() choose any T <= Tc as a reference T at which we report the dim-6 error. 
			// Generally T = Tc is NOT good for this estimate since at Tc, dim-6 effect may be large enough that the broken minimum vanishes
			// and we can't really calculate relative shift in v/T
			{
				EffPot<Complex> effPotNew(params3D);

				bool bDim6 = true;
				// Find LOCAL minimum near to the original one
				ParameterMap newMinimum = effPotNew.FindLocalMinimum(scanner.loopOrderVeff, bDim6, v, x);

				double vNew = GetFromMap(newMinimum, "v");
				double xNew = GetFromMap(newMinimum, "x");
				//std::cout << "(" << vByT << "," << xByT << "); new = (" << vNew/sqrt(T) << "," << xNew/sqrt(T) << ")\n";
 
				fieldShiftsDim6[0] = (vByT - vNew/sqrt(T)) / vByT;
				fieldShiftsDim6[1] = (xByT - xNew/sqrt(T)) / xByT;
			}


			// Results at this temperature
			ParameterMap results = {
				{"T", T},
				// Convert phisq to 4D units and normalize with T. This gives ~(v/T)^2 / 2
				{"phisqByT2", phisq / T},
				{"vByT", vByT},
				{"xByT", xByT},
				{"vShiftDim6", fieldShiftsDim6[0]},
				{"xShiftDim6", fieldShiftsDim6[1]},
				{"Veff.re", GetFromMap(minimum, "Veff.re")},
				{"Veff.im", GetFromMap(minimum, "Veff.im")},
				{"isPerturbative", (int) bIsPerturbative},
				{"warningsMinimization", warningsMinimization},
				{"warningsDerivatives", warningsDerivatives}
			};	
			scanner.resultsForT.push_back(results);

			// ---- Some optimization if we don't care about the v = 0 phase

			// Minimum v/T value before we consider it to be symmetric phase
			double symmPhaseThreshold = 1e-3;

			if (abs(vByT) > symmPhaseThreshold) bFoundBrokenPhase = true;

			// Don't apply any v = 0 optimizations if we haven't found the broken phase yet (avoid spurious effects at small T)
			if (!bFoundBrokenPhase) continue;

			if (bInSymmetricPhase && abs(vByT) > symmPhaseThreshold) {
				// Now bInSymmetricPhase was set to true at previous T
				// but current T is no longer in symmetric phase
				// ==> very weird behavior since we scan from low-T to high-T!
				// This probably means that the minimization failed in the previous point, so reset counter
				symmPhasePoints = symmPhasePointsMax;
				bInSymmetricPhase = false;

			} else if (abs(vByT) < symmPhaseThreshold) {

				bInSymmetricPhase = true;
				if ( scanner.bStopAtSymmetricPhase) {
					symmPhasePoints--;
					if (symmPhasePoints < 0) break;
				}
			}
			// End v = 0 phase optimization
		
		} // end T loop 


		if (scanner.bWriteAtEachTemperature) {
			scanner.WriteTemperatureData();
		}

		scanner.FindTransitionPoints();

	} // end parameter scan loops

	double seconds = Duration(TimeNow() - startTime);
	std::cout << "Scan complete, did " << pointCount << " points total. Time taken: " << seconds << "s.\n";

	return 0;
}
