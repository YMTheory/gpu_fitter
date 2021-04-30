#ifndef gammaResponse_h
#define gammaResponse_h

#include <string>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

using namespace std;

class gammaResponse
{
    public:
        gammaResponse(string name, int nBins, double peMin, double peMax);
        ~gammaResponse();

    public:
        void LoadData();
        void LoadPrmBeta();
        
        void preCalculation();
        void Prediction();
        double SampleGamEnergy(int index);
        void calcGamResponse();

        double GetChi2();

        void SaveHist();
        
    
    public:
        string GetName()                      {return m_name;}
        double GetEtrue()                     {return m_Etrue;}
        double GetPEData()                    {return m_totpeData;}
        double GetNonlData()                  {return m_nonlData;}
        double GetNonlCalc()                  {return m_nonlCalc;}
        double GetNonlErr()                   {return m_nonlDataErr;}
        double GetResData()                   {return m_resData;}
        double GetResCalc()                   {return m_resCalc;}

        double GetAmp()                       {return m_amp;}
        void SetAmp(double amp)               {m_amp = amp;}

    private:
        string m_name;
        int m_nBins;
        double m_peMin;
        double m_peMax;
        double m_Etrue;
        double m_totpeData;
        double m_totpeCalc;
        double m_totpeSigmaCalc;
        double m_nonlData;
        double m_nonlDataErr;
        double m_nonlCalc;
        double m_Evis;
        double m_resData;
        double m_resDataErr;
        double m_resCalc;

        bool m_loadData;
        bool m_loadPrm;

        double m_amp; // TF1 gaussian amplitude

        const int m_nEvents  = 5000;
        const int m_nSamples = 5000;
        TH2D* hPrmElec;
        TH2D* hPrmPosi;

        TH1D* hCalc;
        TH1D* hData;
        TF1*  func;

        double m_pemean[5000];
        double m_pesigma[5000];

    private:
        bool m_doSpecFit;
};


#endif
