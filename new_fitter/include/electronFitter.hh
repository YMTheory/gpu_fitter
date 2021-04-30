#ifndef _electronFitter_h
#define _electronFitter_h

#include <TGraphErrors.h>

class electronFitter
{
    public:
        electronFitter();
        ~electronFitter();

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

        static int m_nData;

        static TGraphErrors* gData;
};


#endif
