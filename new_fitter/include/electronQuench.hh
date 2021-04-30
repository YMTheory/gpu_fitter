#ifndef _ELECTRONQUENCH_H
#define _ELECTRONQUENCH_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>

using namespace std;

class electronQuench
{
    public:
        electronQuench();
        ~electronQuench();
    
    public:
        // parameter configuration
        static void setBirk1(double birk1)             { m_birk1 = birk1;}
        static double getBirk1()                       {return m_birk1;}
        static void setkA(double val)                  { m_kA = val; }
        static double getkA()                          {return m_kA;}
        static void setEnergyScale(double scale)       { m_scale = scale; }
        static double getEnergyScale()                 {return m_scale;}

        static void LoadStopPowData();
        static double Integral_BirkLaw(double E);

        static void LoadNLData();
        static void LoadScintPE();
        static void Update();

        static double ScintillatorNL (double eTrue);

        static double ScintillatorPE (double eTrue);

        static void Plot ();


    private:
        static double m_birk1;
        static constexpr double m_birkHigh = 75e-4;
        static constexpr double m_birkLow  = 55e-4;
        static double m_kA;
		static const unsigned int m_nKb          = 100;
		static const unsigned int m_nSamples     = 2000;
		static constexpr double m_samplingResol  = 0.01;   // bining: 50keV/bin
		static double  m_kBResid;
        static double p0 ;
        static double p1 ;
        static double p2 ;
        static double p3 ;

        static double m_scale;

        static double m_edep[1000];
        static double m_nonl[1000];

        static double m_quenching_energy[m_nKb][m_nSamples];
        static double* m_quenching_energy_low;
        static double* m_quenching_energy_high;
		static double m_quenchingShape1 [m_nKb][m_nSamples];
		static double* m_quenchingShape1_lowKb;
		static double* m_quenchingShape1_higKb;

        static bool m_loadStopPowData;
        static bool m_loadNLData;
        static bool m_loadScintPE;

        static std::vector<double> m_Etrue;
        static std::vector<double> m_StopPow;

        static double m_simEtrue[900];
        static double m_simScintPE[900];

        static double ScintillatorShape    (double eTrue);
        static double SimulationNLShape    (double eTrue);
        static double IntegralNLShape      (double eTrue);
        static double EmpiricalNLShape     (double eTrue);
        static double SimulationNLCalcShape(double eTrue);

        static TGraph* gNPE_elec;
};

#endif
