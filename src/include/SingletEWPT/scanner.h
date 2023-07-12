#ifndef SCANNER_H
#define SCANNER_H

#include "common.h"
#include "effpot_T0.h"
#include "renormalization.h"


/* Class ParameterGrid - The grid over which a parameter is varied in scanning. [min, max], inclusive on both ends.
Can be log-scale (base 10). */
class ParameterGrid {

private:
    bool bLogScale = false;
    // Uniform grid spacing
    double uniformDelta;
    // For log scale: now the grid points A_i are so that log10 A_(i+1) - log10 A_i = const = b. Store here 10^b.
    double logScaleFactor = 10;

    // How many points in the range (so there will be len-1 bins)
    int len;

public:

    double min, max;

    inline int Length() const { return this->len; }

    ParameterGrid() {}

    /* Construct a grid with uniform spacing delta. If logScalePoints > 1, will instead use a log-scale grid with this many points. 
    For log-grids the delta is not used. */ 
    ParameterGrid(double minIn, double maxIn, double delta, int logScalePoints) {
        this->min = minIn;
        this->max = maxIn;
        this->uniformDelta = delta;
        if (delta <= 0) {
            Die("Cannot have delta <= 0!!\n", 862);
        }

        if (logScalePoints > 1 && abs(max - min) > 1e-14) {
            this->bLogScale = true;
            this->len = logScalePoints;
        } else {
            this->bLogScale = false;
            this->len = 1 + (max - min) / uniformDelta;
        }
        

        if (bLogScale && (min <= 0 || max <= 0)) {
            Die("Negative start value in parameter range, cannot use log-scale!!\n", 863);
        }

        if (bLogScale) {
            // Calculate the log-scale grid spacing 
            double b = (log10(max) - log10(min)) / len;
            logScaleFactor = pow(10.0, b); 
        }

    } // End constructor
    
    // Return next point in the grid
    double GetNextPoint(double current) const {
        
        double next;
        if (bLogScale) {
            // A_(i+1) = 10^b A_i
            next = current * logScaleFactor;
        } else {
            next = current + uniformDelta;
        }

        // Dodge floating point rounding
        if (abs(next - max) < 1e-14)
            next = max;

        return next;
    }

    friend std::ostream& operator<<(std::ostream& out, const ParameterGrid &grid) {
        std::string header;
        double spacing;
        if (grid.bLogScale) {
            header = "Log10 scale grid";
            spacing = grid.logScaleFactor;
        }
        else {
            header = "Uniform grid";
            spacing = grid.uniformDelta;
        }
        out << "{ " << header << ": min = " << grid.min << " , max = " << grid.max 
                    << " , delta = " << spacing << " , points = " << grid.len;
        out << " }\n";
        return out;
    }
    
};


/* Class Scanner - Contains control parameters for scanning, eg. parameter ranges and various options */
class Scanner {
private: 
    // Minimum discontinuity in a field (v/T) until the jump is considered a phase transition
    double jumpThreshold = 0.3;

    // How much to reduce temperature from Tc when calculating dim-6 error estimate
    double dim6TemperatureFraction = 0.8; 

    // RG scale to use in 4D -> 3D matching? 0 means we use the optimal choice of approx. 7T 
    double matchingScaleOverride = 0;

    // RG scale to use in 3D EFT? 0 means default value of T 
    double scale3DOverride = 0;

    bool bUseDefaultMatchingScale = true;
    bool bUseDefault3DScale = true;

    double currentTemperature = 0.0;

    // Only report phase transitions where Higgs field jumps from v=0 to nonzero (or vice versa)?
    bool bOnlySearchEWPT = false;


public:
    Scanner(std::string configFileName) {
        ReadScannerParams(configFileName);
    }


    inline void SetTemperature(double T) { currentTemperature = T; };

    // Read parameters from file
    void ReadScannerParams(std::string fname);

    void PrintScanner() {
        // Enums
        std::cout << "loopOrderMS : " << static_cast<int>(loopOrderMS) << "\n";
        std::cout << "loopOrderDR : " << static_cast<int>(loopOrderDR) << "\n";
        std::cout << "loopOrderVeff : " << static_cast<int>(loopOrderVeff) << "\n";
        std::cout << "loopOrderVeffT0 : " << static_cast<int>(loopOrderVeffT0) << "\n";
        // Booleans
        std::cout << "bCalculateCondensates : " << bCalculateCondensates << "\n";
        std::cout << "bStopAtSymmetricPhase : " << bStopAtSymmetricPhase << "\n";
        std::cout << "bSolveBetas : " << bSolveBetas << "\n";
        std::cout << "bOnlySearchEWPT : " << bOnlySearchEWPT << "\n";


        // Miscellaneous
        std::cout << "jumpThreshold : " << jumpThreshold << "\n";
        std::cout << "dim6TemperatureFraction : " << dim6TemperatureFraction << "\n";
        std::cout << "matchingScaleOverride : " << matchingScaleOverride << "\n";
        std::cout << "scale3DOverride : " << scale3DOverride << "\n";

        // Scanning ranges. Also count how many parameter points we have
        long points = 1;
        for (auto const& x : scanningRange) {
            std::cout << "Range " << x.first << ":\n";
            std::cout << x.second;
            if (x.first != "T") {
                points *= x.second.Length(); 
            }
        }
        std::cout << points << " parameter points in total (not counting the T-loop)\n";

        // Print info about what renormalization group scales are used
        if (bUseDefaultMatchingScale) {
            std::cout << "RG scale for 4D -> 3D matching: ~ 7.055T\n";
        } else {
            std::cout << "RG scale for 4D -> 3D matching: " << matchingScaleOverride << "\n";
        }
        if (bUseDefault3DScale) {
            std::cout << "RG scale for 3D calculations: T\n";
        } else {
            std::cout << "RG scale for 3D calculations: " << scale3DOverride << "\n";
        }
    }

    void StartTemperatureLoop() {
        resultsForT.clear();
        // Initialize the T-results vector to correct size
        /*
        const int NTEMPERATURE = 1 + ceil( 
            ( GetFromMap(scanningRange, "T_max") - GetFromMap(scanningRange, "T_min") ) / GetFromMap(scanningRange, "T_delta")
        );
        resultsForT.resize(NTEMPERATURE);
        */
    }

    /* Minimize the T=0 effective potential and check if EW minimum is the global one.
    In practice, just checks that v > 0 in the minimum. Will write the outcome to a separate file. */
    bool CheckT0Stability(const ParameterMap &MSParams);

    // Scan resultsForT vector for transitions and calculate discontinuities etc across the transition
    void FindTransitionPoints();

    // Write labels for transition and T=0 data into files
    void WriteDataLabels();

    // Write a vector to file in a column-by-column format
    void AppendToFile(std::string fname, std::vector<double> data);

    // Writes resultsForT vector to file with a sensible ordering
    void WriteTemperatureData();

    // Current input parameters for the scanner
    ParameterMap currentInput;

    // Container for storing results at each temperature (Veff value, minimum location etc)
    std::vector<ParameterMap> resultsForT;

    std::map<std::string, ParameterGrid> scanningRange;
    ELoopOrder loopOrderMS;
    ELoopOrder loopOrderVeff;
    ELoopOrder loopOrderVeffT0;
    ELoopOrderDR loopOrderDR;
    bool bCalculateCondensates;
    // Break T-loop once (v,x) == (0, 0) phase is found?
    bool bStopAtSymmetricPhase;
    /* Store measurements of Veff, global minimum, condensates etc at each temperature?
    Will use a separate file for these */ 
    bool bWriteAtEachTemperature = false;
    // Do 1-loop running to matching scale or not?
    bool bSolveBetas = false;

    /* Returns RG scale for 4D -> 3D matching based on scanner options.
    If matchingScaleOverride == 0, uses optimal choice of approx 7T, where T is the current scanner temperature. 
    Otherwise uses value in matchingScaleOverride. */
    double GetMatchingScale() const;

    /* Returns RG scale for 3D matching based on scanner options.
    If matchingScaleOverride == 0, uses optimal choice of approx 7T, where T is the current scanner temperature. 
    Otherwise uses value in matchingScaleOverride. */
    double Get3DScale() const;

    std::string transitionsFileName = "transitions.dat";
    std::string temperatureDataFileName = "temperature_data.dat";
};

#endif
