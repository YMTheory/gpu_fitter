#ifndef _gammaFitter_h
#define _gammaFitter_h

#include "gammaData.hh"

#include <string>

class gammaFitter
{

    public:
        gammaFitter();
        ~gammaFitter();

        void Initialize();

        static double GetChi2(const double *par);

        static int Minimization();

        static void Plot();

    private:
        static bool m_DoFit;
        static double m_chi2Min;
        static int m_nParameters;
        static double m_bestFit[20];
        static double m_fitError[20];

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

        static int m_nData;
        static std::string source_name[20];

        static gammaData* gammaData_array[20];
};

#endif
