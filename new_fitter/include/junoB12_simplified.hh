#ifndef junoB12_simplified_h
#define junoB12_simplified_h

#include <TF1.h>

using namespace std;

class junoB12_simplified
{
    
    public:
        junoB12_simplified(int nBinsData, double fitMinE, double fitMaxE);
        ~junoB12_simplified();

    public:
        
        void Initialize();
        void LoadDataSpec();
        void LoadTheoSpec();

        void ApplyScintillatorNL();
        void ApplyResolutionNL();

        void Normalize();
        double GetChi2();

        void Plot();

    private:
        bool m_loadData;
        bool m_loadTheo;

        int m_nBin;
        int m_nBinData;
        double m_fitMinPE;
        double m_fitMaxPE;
        double m_eMin;
        double m_eMax;
        double m_peMin;
        double m_peMax;


        double m_binWidth;
        double m_peBinWidth;
        double m_binCenter[1500];
        double m_peBinCenter[1500];
        double m_eData[100];
        double m_eDataErr[100];
        double m_eTru[1500];
        double m_eVis[1500];
        double m_eTheo[100];

    private:
        TF1* gaus;

};



#endif
