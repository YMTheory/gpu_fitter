#ifndef _spectrumFitter_h
#define _spectrumFitter_h

#include "junoSpectrum.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

class spectrumFitter
{
    public:
        spectrumFitter();
        ~spectrumFitter();

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
        static junoSpectrum* junoB12;
};

#endif
