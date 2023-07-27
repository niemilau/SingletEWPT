#ifndef RENORM_H
#define RENORM_H

#include "common.h"

// Loop order for MS-bar renormalization. Using a separate typename in case this gets extended later 
using ELoopOrderMS = ELoopOrder;


/* Experimental values for Standard Model observables */
namespace ExperimentalInput {

    // Particle masses in GeV
    const double Mt = 172.76;
    const double MW = 80.379;
    const double MZ = 91.1876;
    const double MH = 125.10;

    // Strong fine-structure constant in MS-bar at Z-pole
    const double alphaQCD = 0.1181;
    // QCD coupling constant squared
    const double gs2 = 4.0*PI*alphaQCD;
    // Fermi constant (GeV^-2)
    const double Gf = 1.1663787e-5;
    // Classical fine-structure constant in the Thomson limit 
    const double alphaEM = 1.0 / 137.035999679;
    // Hadronic contribution to Delta alpha at Z-pole, light quarks only (not used)
    // const double DeltaAlphaHadronic = 276.26e-4;

    // 'Tree-level' SU(2) gauge coupling from Fermi constant and W mass
    const double g0sq = 4.0 * sqrt(2.0) * Gf * MW * MW;
};


/* Class Renormalization. The purpose of this class is to convert SM + singlet inputs {Mh1, Mh2, a2, b3, b4, sinTheta}
and experimental values for MW, MZ etc to MS-bar renormalized parameters {ytsq, g1sq, g2sq, g3sq, msqPhi, lambda, b1, b2, b3, b4, a1, a2}.
The conversion can be done at tree level or at 1-loop level. 
For the latter case the class calculates 1-loop corrected propagators and matches the pole masses at external momentum p = M. */
class Renormalization {

private:
    /****** SM + singlet specific input parameters ******/

    // Pole masses of scalar excitations
    double Mh1, Mh2;
    // Couplings etc (in MS-bar)
    double a2, b3, b4, sinTheta;
    // Mixing angle, solved from sinTheta
    double theta;
    // shorthands for cos and sin of theta
    double ct, st;
    // MS-bar scale at which the above couplings are given
    double scale;

    /** Constants **/
    const double Mt = ExperimentalInput::Mt;
    const double MW = ExperimentalInput::MW;
    const double MZ = ExperimentalInput::MZ;
    
    // Alias for g3sq, SU(3) coupling
    const double gs2 = ExperimentalInput::gs2;
    // SU(3) coupling
    const double g3sq = gs2;
    // Fermi constant (GeV^-2)
    const double Gf = ExperimentalInput::Gf;
    // 'Tree-level' SU(2) gauge coupling from Fermi constant and W mass
    const double g0sq = 4.0 * sqrt(2.0) * Gf * MW * MW;
    
    /***** Self-energies *****/
    double selfEnergyH1, selfEnergyH2, selfEnergyW, selfEnergyZ, selfEnergyTop;
    // SU(2) gauge coupling correction: g^2 = g0^2 ( 1 + Delta g^2) 
    double gsqDelta;

public:

    // Constructor that just calls SetParameters
    Renormalization(ParameterMap singletInputs, double RGScale) {
        SetParameters(singletInputs, RGScale);
    } 

    // Fix SM + singlet inputs
    void SetParameters(ParameterMap singletInputs, double RGScale) {
        double mass1 = GetFromMap(singletInputs, "Mh1");
        double mass2 = GetFromMap(singletInputs, "Mh2");
        if (mass1 < 0. || mass2 < 0.) {
            std::cout << "!!! Negative mass input, got negative: M1 = " << mass1 << " ; M2 = " << mass2 << "\n";
            Die("Exiting...\n", 89);
        }
        
        Mh1 = mass1;
        Mh2 = mass2;

        a2 = GetFromMap(singletInputs, "a2");
        b3 = GetFromMap(singletInputs, "b3");
        b4 = GetFromMap(singletInputs, "b4");
        sinTheta = GetFromMap(singletInputs, "sinTheta");
        this->scale = RGScale; 
        
        // Mixing angle
        theta = GetTheta(sinTheta);
        this->st = sin(theta);
        this->ct = cos(theta);
    }

    // Convert Sin(theta) to the angle theta
    inline double GetTheta(double sinTheta) {
        // Needs to be within |sinTheta| < 1/sqrt(2); this does not restrict the available parameter space
        if (std::abs(sinTheta) > 1./std::sqrt(2.)) {
            std::cout << "!!! Invalid sinTheta = " << sinTheta << " ; needs to be |sinTheta| < 1/sqrt(2)\n";
            Die("Exiting...\n", 90);
        }
        return asin(sinTheta);
    }

    /* Calculate 1-loop self energies. These are complicated functions of the renormalized parameters, 
    but inside loop corrections I replace renormalized masses with the pole masses (known from input). The error is of higher order. */
    void CalculateSelfEnergies();

    /* Get MS-bar parameters corresponding to given input */
    ParameterMap CalcMS(ELoopOrderMS loopOrder);

    // Calculate top Yukawa coupling squared
    double CalcYtsq(ELoopOrderMS loopOrder);

    // Calculate SU(2) gauge coupling squared (Fermi EFT method)
    double Calcg2sq(ELoopOrderMS loopOrder);

    // Calculate U(1) gauge coupling squared
    double Calcg1sq(ELoopOrderMS loopOrder);


    /* Scalar potential parameters. Here phisq = phi^+ phi and potential is
    V = msq_phi phisq + lambda phisq^2 + 1/2 a1 S phisq + 1/2 a2 S^2 phisq 
        + b1 S + 1/2 b2 S^2 + 1/3 b3 S^3 + 1/4 b4 S^4.
    */ 

    // Calculate Higgs quartic self-interaction
    double CalcLambda(ELoopOrderMS loopOrder);

    // Calculate Higgs quadratic mass parameter
    double CalcMsqPhi(ELoopOrderMS loopOrder);

    // Calculate singlet quadratic mass parameter
    double Calcb2(ELoopOrderMS loopOrder);

    // Calculate Higgs-singlet cubic interaction
    double Calca1(ELoopOrderMS loopOrder);

    // Calculate Singlet linear coupling
    double Calcb1(ELoopOrderMS loopOrder);


    /*** Some static functions. Would probably be better to have these elsewhere! ***/

    /* Run MS-bar parameters to a different scale using 1-loop beta functions. 
    Uses a Runge-Kutta4 integrator from boost/odeint */ 
    static ParameterMap RunToScale(double scaleOut, const ParameterMap &params);

    // Checks if the couplings are within NAIVE perturbativity bounds 
    static bool CheckPerturbativity(const ParameterMap &params);

    // Test function for Passarino-Veltman functions
    void TestPaVe();

    // Self energy computations, at external momentum k
    double CalcSelfEnergyW(double k);
    double CalcSelfEnergyZ(double k);
    double CalcSelfEnergyH1(double k);
    double CalcSelfEnergyH2(double k);
    // Top quark self energy: this is actually scalar part + vector part and scaled dimensionless (Minkowskian signature)
    double CalcSelfEnergyTop(double k);

    // SU(2) gauge coupling g^2 = g0^2 (1 + Delta gsq); this calculates Delta gsq. Eq. (A28) in 2103.07467
    inline double CalcgsqCorrection() {
        double MW2 = MW*MW;
        double MZ2 = MZ*MZ;
        double delta_r_remainder = g0sq / (16.0*PI*PI) * ( 4.0 * log(scale*scale / MW2) + (7.0/2.0*MZ2 / (MZ2 - MW2) - 2) * log(MW2/MZ2) + 6.0 );
        return ( CalcSelfEnergyW(MW) - CalcSelfEnergyW(0.0) ) / MW2 - delta_r_remainder;
    }

private:

    /******** Passarino-Veltman integrals + some extra. See PassarinoVeltman.cpp for implementations. ********/

    // class DivergentIntegral: used to represent divergent integrals in MS-bar scheme. The divergent part is everything proportional to 1/eps 
    template <typename Float>
    class DivergentIntegral {
    private:
        Float finitePart, divergentPart;
    public:
        DivergentIntegral(const Float &fin, const Float &div) {
            finitePart = fin;
            divergentPart = div;
        }
        Float fin() const { return finitePart; }
        Float div() const { return divergentPart; }
        // Operator overloads
        DivergentIntegral operator+(const DivergentIntegral &obj) const {
            return DivergentIntegral(finitePart + obj.fin(), divergentPart + obj.div());
        }
        // Remember that compiler automatically assumes this object on the LHS and the other object (argument) on the RHS
        DivergentIntegral operator-(const DivergentIntegral &obj) const {
            return DivergentIntegral(finitePart - obj.fin(), divergentPart - obj.div());
        }
        // DivergentIntegral * Float
        DivergentIntegral operator*(const Float &x) const {
            return DivergentIntegral(finitePart * x, divergentPart * x);
        }
        // Float * DivergentIntegral
        friend DivergentIntegral operator*(const Float &x, const DivergentIntegral &obj) { 
            return obj * x;
        }
        // DivergentIntegral / Float
        DivergentIntegral operator/(const Float &x) const {
            return DivergentIntegral(finitePart / x, divergentPart / x);
        }
    };

    // Lazy typedefs
    using Float = double;
    using Integral = DivergentIntegral<Float>;
    using Complex = std::complex<Float>;

    // Small-mass cutoff, to prevent NaN behavior due to m==0
    const Float smallMassThreshold = 1e-12;
    inline bool IsZeroMass(const Float &m) {
        return (abs(m) < smallMassThreshold);
    }
    // 1/(4pi)^2
    const Float pi4sqInv = 1.0 / (16.0*PI*PI);


    // A0(m) = int_p 1 / (p^2 - m^2)
    Integral A0(const Float &m);

    /* Function F(k, m1, m2) from Appendix B of Böhm, Spiesberger, Hollik 1986. It is finite, but can have an imaginary part. 
    This appears in B0 only. To avoid issues at m1 = m2, it is convenient to define Fb(k, m1, m2), which is F(k, m1, m2) minus first row of their eq. B.1. 
    Here I use Fb(k, m1, m2) exclusively, explicitly: Fb(k, m1, m2) = F(k, m1, m2) - 1 - ( (m1^2 - m2^2) / k^2 - (m1^2 + m2^2) / (m1^2 - m2^2) )*ln(m2/m1)  */
    Complex Fb(const Float &k, const Float &m1, const Float &m2);

    /* B0(k, m1, m2) = int_p 1 / [(p^2 - m1^2)((p+k)^2 - m2^2)].
    This develops an imaginary part for k^2 > (m1+m2)^2, which is dropped here! (we only need real parts of the self energies)  */
    Integral B0(const Float &k, const Float &M1, const Float &M2);

    // B1 and B00. See Appendix A.3 in https://arxiv.org/pdf/2103.07467.pdf for definitions
    Integral B1(const Float &k, const Float &m1, const Float &m2);
    Integral B00(const Float &k, const Float &m1, const Float &m2);

    // K(k, m1, m2, a, b, c) = int_p (ap^2 + b k^mu p_mu + c) / [(p^2 - m1^2)((p+k)^2 - m2^2)]
    Integral K(const Float &k, const Float &m1, const Float &m2, const Float &a, const Float &b, const Float &c) {
        return a * A0(m2) + b*k*k * B1(k, m1, m2) + (a*m1*m1 + c) * B0(k, m1, m2);
    }
};



#endif
