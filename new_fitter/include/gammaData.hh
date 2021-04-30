#ifndef _gammaData_h
#define _gammaData_h

#include <iostream>
#include <string>
#include <vector>
#include <TH2D.h>
#include <TH1D.h>

using namespace std;

class gammaData
{
    public:
        gammaData(std::string name,
                  double minPE,
                  double maxPE,
                  int nbins
                 );

        ~gammaData();

    
    public:
        double GetEtrue         () {return m_Etrue;}
        double GetEvis          () {return m_Evis;}
        double GetPEData        () {return m_totpeData;}
        double GetPECalc        () {return m_totpeCalc;}
        double GetCerPECalc     () {return m_cerPE;}
        double GetSctPECalc     () {return m_sctPE;}
        double GetCerPEData     () {return m_cerPEData;}
        double GetSctPEData     () {return m_sctPEData;}
        double GetNonlPred      () {return m_nonlCalc;}
        double GetNonlData      () {return m_nonlData;}
        double GetNonlDataErr   () {return m_nonlDataErr;}
        double GetResolPred     () {return m_resCalc;}
        double GetResolData     () {return m_resData;}
        double GetResolDataErr  () {return m_resDataErr;}

        double GetNonlPred1     () {return m_nonlCalc1;}

        void LoadData           ();
        double GetChi2          ();

        void SaveHist           ();

        void LoadGammaData ();
        void LoadGammaPEComp();
        void LoadPrimElecDist();
        void LoadPrimElecSamples();

        void calcGammaResponse();

    
    private:

        std::string m_name;
        double m_minPE;
        double m_maxPE;
        int m_nbins;

        bool m_loadData;

        double m_Etrue;
        double m_Evis;
        double m_totpeCalc;
        double m_totpeData;
        double m_nonlCalc;
        double m_nonlData;
        double m_nonlDataErr;
        double m_resCalc;
        double m_resData;
        double m_resDataErr;

        // back-up
        double m_nonlCalc1;
        double m_sctPE;
        double m_cerPE;
        double m_sctPEData;
        double m_cerPEData;

    private:
        double m_scale;
        double m_nuGamma;

    private:
        static const unsigned int m_nMaxPdf = 1000;
        double m_pdf_eTrue[m_nMaxPdf];
        double m_pdf_prob[m_nMaxPdf];
        double m_max_eTrue;
        double m_pdf_eTrue_posi[m_nMaxPdf];
        double m_pdf_prob_posi[m_nMaxPdf];

    private:
        TH2D* elec_hist;
        static const unsigned int m_nSamples = 5000;
        double m_mean[m_nSamples];
        double m_meanpe[m_nSamples];
        //double m_sigmape[m_nSamples];
        TH1D* hEmean;
        TH1D* hPEmean;
        TH1D* hCerPEmean;
        TH1D* hSctPEmean;

    private:
        static std::string m_calcOption;  // prmelec ; twolayer
        static std::string m_nonlMode;  // histogram; analytic

    public:
        std::string getNonlMode()  {return m_nonlMode;}

};

#endif
