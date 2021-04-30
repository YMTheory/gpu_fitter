#include "spectrumFitter.hh"
#include "junoParameters.hh"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TSystem.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

bool spectrumFitter::m_DoFit = false;
double spectrumFitter::m_chi2Min;
int spectrumFitter::m_nParameters;
double spectrumFitter::m_bestFit[20];
double spectrumFitter::m_fitError[20];

junoSpectrum* spectrumFitter::junoB12;

spectrumFitter::spectrumFitter()
{;}

spectrumFitter::~spectrumFitter()
{
    delete junoB12;
}

void spectrumFitter::Initialize()
{
    junoB12 = new junoSpectrum(14400, 100,
                               3, 1,
                               0, 15,
                               1, 14,
                               junoParameters::m_nonlMode,
                               "B12");
    junoB12->LoadData();
}

double spectrumFitter::GetChi2(const double *par)
{
    double chi2 = 0;

    electronQuench::setkA(par[0]);
    electronQuench::setEnergyScale(par[2]);
    electronCerenkov::setkC(par[1]);
    electronCerenkov::setEnergyScale(par[2]);

    chi2 += junoB12->GetChi2();
    return 1;
}

int spectrumFitter::Minimization()
{
    ROOT::Math::Minimizer* minimum = 
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(1);

    // function wrapper for Minimizer_new
    ROOT::Math::Functor f(&GetChi2, 4);
    double step[4];
    step[0] = 0.001;
    step[1] = 1e-6;
    step[2] = 0.001;
    step[3] = 0.1;

    // start point :
    double variable[4];
    variable[0] = 0.96;
    variable[1] = 1.0;
    variable[2] = 3300.371/2.223;

    minimum->SetFunction(f);

    minimum->SetVariable(0, "kA", variable[0], step[0]);
    minimum->SetVariable(1, "kC", variable[1], step[1]);
    minimum->SetVariable(2, "energyScale", variable[2], step[2]);

    minimum->SetVariableLimits(0, 0.9, 1.1);
    minimum->SetVariableLimits(1, 0.9, 1.2);
    minimum->SetVariableLimits(2, 1400, 1600);

    //minimum->FixVariable(0);
    //minimum->FixVariable(1);
    //minimum->FixVariable(2);

    minimum->Minimize();

    m_chi2Min = minimum->MinValue();
    std::cout << "=====> Minimized chi2 = " << m_chi2Min << std::endl;

    m_nParameters = minimum->NDim();
    std::cout << "=====> Total fitting parameters number : " << m_nParameters << std::endl;

    const double *X = minimum->X();
    const double *Err = minimum->Errors();
    for (int iPar=0; iPar<m_nParameters; iPar++) {
        m_bestFit[iPar] = X[iPar];
        m_fitError[iPar] = Err[iPar];
        cout << "Parameter " << iPar << " " << m_bestFit[iPar] << " " << m_fitError[iPar] << std::endl;
    }

    m_DoFit = true;
    return 1;
}

void spectrumFitter::Plot()
{;}
