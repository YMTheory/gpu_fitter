#ifndef junoNLChiFunction_h
#define junoNLChiFunction_h

#include "gammaData.hh"
#include "junoSpectrum.hh"
#include "gammaResponse.hh"
#include "junoB12.hh"
#include "junoB12_simplified.hh"

#include <TMinuit.h>

class junoNLChiFunction
{
    public:
        junoNLChiFunction();
        ~junoNLChiFunction();

        void LoadData();
        double GetChiSquare           ( double maxChi2 = 100000 );
        static void SetParameters     ( double *par );
        static double GetChi2         ( double maxChi2 = 100000 );  
    
        static void GammaPlot         ();
        static void GammaPEPlot       ();
        static void Plot              ();

    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

        TMinuit* junoNLMinuit;

        static double m_chi2;
        static double m_chi2Min;
        static double m_chi2Gam;
        static double m_chi2B12;
        static double m_chi2C11;
        static double m_chi2C10;
        static int m_nParameter;
        static double m_bestFit[20];
        static double m_bestFitError[20];
        static bool m_DoFit;
        

    private:
        static gammaData* Cs137Data;
        static gammaData* Mn54Data;
        static gammaData* nHData;
        static gammaData* K40Data;
        static gammaData* Co60Data;
        static gammaData* Tl208Data;
        static gammaData* nC12Data;
        static gammaData* O16Data;
        static gammaData* nFe56Data;
        static gammaData* gamma4440Data;
        static gammaData* gamma511Data;

        static junoSpectrum* junoB12data;

        static junoB12_simplified* b12data;
        //static junoB12* b12data;

        static int m_nData;
        static int m_nGam;
        static std::string source_name[20];

        //static gammaData* gammaData_array[20];

        static gammaResponse* gammaData_array[10];
        static gammaResponse* Cs137;
        static gammaResponse* Mn54;
        static gammaResponse* K40;
        static gammaResponse* Ge68;
        static gammaResponse* Co60;
        static gammaResponse* nH;
        static gammaResponse* AmBe;
        static gammaResponse* nC12;
        static gammaResponse* AmC;

        static bool m_doGamFit;
        static bool m_doB12Fit;

    private:
        static string m_nonlMode;

    
};

#endif
