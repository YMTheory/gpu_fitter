#include "electronFitter.hh"
#include "electronResponse.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TSystem.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>

int electronFitter::m_nData;
bool electronFitter::m_DoFit = false;
double electronFitter::m_chi2Min;
int electronFitter::m_nParameters;
double electronFitter::m_bestFit[20];
double electronFitter::m_fitError[20];

TGraphErrors* electronFitter::gData;

electronFitter::electronFitter()
{
}

electronFitter::~electronFitter()
{
    delete gData;
}


void electronFitter::Initialize()
{
    std::cout << "+++++++++++++++++++ >> Begin of electronFitter Initialization ++++++++++++ << " << std::endl;

    electronResponse::loadMinSimElecNonl();
    gData = electronResponse::gMinElecNonl;
}


double electronFitter::GetChi2(const double *par) {
    double chi2 = 0;

    //if (junoParameters::scintillatorParameterization == kSimulationCalc) {
        electronQuench::setkA(par[0]);
        electronQuench::setEnergyScale(par[2]);
        electronCerenkov::setkC(par[1]);
        electronCerenkov::setEnergyScale(par[2]);
    //}
    
    //if (junoParameters::scintillatorParameterization == kSimulation) {
    //    electronQuench::setkA(par[0]);
    //    electronQuench::setBirk1(par[1]);
    //    electronCerenkov::setkC(par[2]);
    //    electronCerenkov::setEnergyScale(par[3]);
    //}

    for (int i=0; i<gData->GetN(); i++) {
        double E = gData->GetPointX(i);
        double data = gData->GetPointY(i);
        double err = gData->GetEY()[i];
        double tmp_nonl = electronResponse::getElecNonl(E);
        
        chi2 += (tmp_nonl - data) * (tmp_nonl - data) / err /err;
    }

    return chi2;
}


int electronFitter::Minimization()
{
    ROOT::Math::Minimizer* minimum = 
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(1);

    // function wrapper for Minimizer_new
    ROOT::Math::Functor f(&GetChi2, 3);

    //if (junoParameters::scintillatorParameterization == kSimulation) {
    //    double step[4];
    //    step[0] = 0.001;
    //    step[1] = 1e-5;
    //    step[2] = 0.001;
    //    step[3] = 0.1;

    //    // start point :
    //    double variable[4];
    //    variable[0] = 0.96;
    //    variable[1] = 6.5e-3;
    //    variable[2] = 1.0;
    //    variable[3] = 3350/2.22;

    //    minimum->SetFunction(f);

    //    minimum->SetVariable(0, "kA", variable[0], step[0]);
    //    minimum->SetVariable(1, "kB", variable[1], step[1]);
    //    minimum->SetVariable(2, "kC", variable[2], step[2]);
    //    minimum->SetVariable(3, "energyScale", variable[3], step[3]);

    //    minimum->SetVariableLimits(0, 0.5, 1.5);
    //    minimum->SetVariableLimits(1, 5.0e-3, 7.5e-3);
    //    minimum->SetVariableLimits(2, 0.5, 1.5);
    //    minimum->SetVariableLimits(3, 1400, 1600);


    //}

    //if (junoParameters::scintillatorParameterization == kSimulationCalc) {
        double step[3];
        step[0] = 0.001;
        step[1] = 0.001;
        step[2] = 0.1;

        // start point :
        double variable[3];
        variable[0] = 1.0;
        variable[1] = 1.0;
        variable[2] = 3350/2.22;

        minimum->SetFunction(f);

        minimum->SetVariable(0, "kA", variable[0], step[0]);
        minimum->SetVariable(1, "kC", variable[1], step[1]);
        minimum->SetVariable(2, "energyScale", variable[2], step[2]);

        minimum->SetVariableLimits(0, 0.5, 1.5);
        minimum->SetVariableLimits(1, 0.5, 1.5);
        minimum->SetVariableLimits(2, 1400, 1600);


    //}


    //minimum->FixVariable(0);
    //minimum->FixVariable(1);
    //minimum->FixVariable(2);
    //minimum->FixVariable(3);

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


void electronFitter::Plot()
{
    if (not m_DoFit) Minimization();

    TGraphErrors* gNonlCalc = new TGraphErrors();
    gNonlCalc->SetName("gNonlCalc");

    int index = 0;
    for(int iData=0; iData<gData->GetN(); iData++) {
        double tmp_E = gData->GetPointX(iData);
        double tmp_pred = electronResponse::getElecNonl(tmp_E);
        gNonlCalc->SetPoint(index, tmp_E, tmp_pred);

        index++;
    }


    gData->SetMarkerStyle(20);
    gData->SetMarkerColor(kBlue+1);
    gData->SetLineColor(kBlue+1);
    gData->SetLineWidth(2);
    gData->SetMarkerSize(1.2);
    gNonlCalc->SetMarkerStyle(21);
    gNonlCalc->SetMarkerColor(kRed+1);
    gNonlCalc->SetMarkerSize(1.2);
    gNonlCalc->SetLineColor(kRed+1);
    gNonlCalc->SetLineWidth(2);
    TCanvas* c1 = new TCanvas("Nonlinearity", "Nonlinearity");
    c1->cd(); c1->SetGrid();
    gData->SetTitle("Nonlinearity Fitting; Etrue/MeV; Nonlinearity");
    //gNonlData->GetYaxis()->SetRangeUser(0.01,0.045);
    gData->Draw("APL");
    gNonlCalc->Draw("P SAME");
    TLegend* led = new TLegend();
    led->SetFillColor(kWhite);
    led->AddEntry(gData, "data", "PL");
    led->AddEntry(gNonlCalc, "calc", "PL");
    led->Draw("SAME");

    c1->SaveAs("ElecNLFakeFit.root");    
}
























