#include "gammaFitter.hh"
#include "gammaData.hh"
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

using namespace std;

bool gammaFitter::m_DoFit = false;
double gammaFitter::m_chi2Min;
int gammaFitter::m_nParameters;
double gammaFitter::m_bestFit[20];
double gammaFitter::m_fitError[20];

gammaData* gammaFitter::Cs137Data;
gammaData* gammaFitter::Mn54Data;
gammaData* gammaFitter::nHData;
gammaData* gammaFitter::K40Data;
gammaData* gammaFitter::Co60Data;
gammaData* gammaFitter::Tl208Data;
gammaData* gammaFitter::nC12Data;
gammaData* gammaFitter::O16Data;
gammaData* gammaFitter::nFe56Data;
int gammaFitter::m_nData;
std::string gammaFitter::source_name[20];
gammaData* gammaFitter::gammaData_array[20];

gammaFitter::gammaFitter()
{}

gammaFitter::~gammaFitter()
{
    delete Cs137Data;
    delete Mn54Data;
    delete nHData;
    delete K40Data;
    delete Co60Data;
    delete Tl208Data;
    delete nC12Data;
    delete O16Data;
    delete nFe56Data;
}

void gammaFitter::Initialize()
{

    cout << endl;
    cout << " ++++++++++++++++++++++++++++++++++++++++++++ >> Begin of Initialization << ++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << endl;

    ProcInfo_t info;;
    gSystem->GetProcInfo(&info);
    float m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    m_nData = 0;

    Cs137Data = new gammaData("Cs137", 700, 1100, 100);
    Cs137Data->LoadData(); 
    source_name[m_nData] = "Cs137"; 
    gammaData_array[m_nData] = Cs137Data;
    m_nData++;
    gSystem->GetProcInfo(&info);
    m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    Mn54Data  = new gammaData("Mn54", 900, 1300, 100);
    Mn54Data->LoadData(); 
    source_name[m_nData] = "Mn54"; 
    gammaData_array[m_nData] = Mn54Data;
    m_nData++;
    gSystem->GetProcInfo(&info);
    m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    K40Data  = new gammaData("K40", 900, 1300, 100);
    K40Data->LoadData(); 
    source_name[m_nData] = "K40"; 
    gammaData_array[m_nData] = K40Data;
    m_nData++;
    gSystem->GetProcInfo(&info);
    m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    nHData  = new gammaData("nH", 900, 1300, 100);
    nHData->LoadData(); 
    source_name[m_nData] = "nH"; 
    gammaData_array[m_nData] = nHData;
    m_nData++;
    gSystem->GetProcInfo(&info);
    m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    Co60Data  = new gammaData("Co60", 900, 1300, 100);
    Co60Data->LoadData(); 
    source_name[m_nData] = "Co60"; 
    gammaData_array[m_nData] = Co60Data;
    m_nData++;
    gSystem->GetProcInfo(&info);
    m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    Tl208Data  = new gammaData("Tl208", 900, 1300, 100);
    Tl208Data->LoadData(); 
    source_name[m_nData] = "Tl208"; 
    gammaData_array[m_nData] = Tl208Data;
    m_nData++;
    gSystem->GetProcInfo(&info);
    m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    nC12Data  = new gammaData("nC12", 900, 1300, 100);
    nC12Data->LoadData(); 
    source_name[m_nData] = "nC12"; 
    gammaData_array[m_nData] = nC12Data;
    m_nData++;
    gSystem->GetProcInfo(&info);
    m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    O16Data  = new gammaData("O16", 900, 1300, 100);
    O16Data->LoadData(); 
    source_name[m_nData] = "O16"; 
    gammaData_array[m_nData] = O16Data;
    m_nData++;
    gSystem->GetProcInfo(&info);
    m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    nFe56Data  = new gammaData("nFe56", 900, 1300, 100);
    nFe56Data->LoadData(); 
    source_name[m_nData] = "nFe56"; 
    gammaData_array[m_nData] = nFe56Data;
    m_nData++;
    gSystem->GetProcInfo(&info);
    m_momery = info.fMemResident/1024;
    cout << "Memory = " << m_momery << " MB" << endl;

    cout << endl;
    cout << " ++++++++++++++++++++++++++++++++++++++++++++ >> End of Initialization << ++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << endl;
}

double gammaFitter::GetChi2(const double *par)
{
    double chi2 = 0;
    
    electronQuench::setkA(par[0]);
    electronQuench::setBirk1(par[1]);
    electronCerenkov::setkC(par[2]);
    electronCerenkov::setEnergyScale(par[3]);

    for(int iSource=0; iSource<m_nData; iSource++) {
        chi2 += gammaData_array[iSource]->GetChi2();
    }

    return chi2;
}

int gammaFitter::Minimization()
{
    //ROOT::Minuit2::Minuit2Minimizer_new minimum (ROOT::Minuit2::kMigrad);
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
    step[1] = 1e-5;
    step[2] = 0.001;
    step[3] = 0.1;

    // start point :
    double variable[4];
    variable[0] = 0.96;
    variable[1] = 6.5e-3;
    variable[2] = 1.0;
    variable[3] = 3350/2.22;

    minimum->SetFunction(f);

    minimum->SetVariable(0, "kA", variable[0], step[0]);
    minimum->SetVariable(1, "kB", variable[1], step[1]);
    minimum->SetVariable(2, "kC", variable[2], step[2]);
    minimum->SetVariable(3, "energyScale", variable[3], step[3]);

    minimum->SetVariableLimits(0, 0.5, 1.5);
    minimum->SetVariableLimits(1, 5.0e-3, 7.5e-3);
    minimum->SetVariableLimits(2, 0.5, 1.5);
    minimum->SetVariableLimits(3, 1400, 1600);

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

void gammaFitter::Plot()
{
    if (not m_DoFit) Minimization();

    TGraphErrors* gNonlData = new TGraphErrors();
    TGraphErrors* gNonlCalc = new TGraphErrors();
    gNonlData->SetName("gNonlData");
    gNonlCalc->SetName("gNonlCalc");

    int index = 0;
    for(int iData=0; iData<m_nData; iData++) {
        std::string source = source_name[iData];
        gammaData* tmpGammaData = gammaData_array[iData];
        //tmpGammaData->calcGammaNPE();
        double tmp_E       = tmpGammaData->GetEtrue();
        double tmp_pred    = tmpGammaData->GetNonlPred();
        double tmp_data    = tmpGammaData->GetNonlData();
        double tmp_dataErr = tmpGammaData->GetNonlDataErr(); 
        //cout << tmp_E << " " << tmp_data << " " << tmp_pred << endl
        gNonlData->SetPoint(index, tmp_E, tmp_data);
        gNonlData->SetPointError(index, 0, tmp_dataErr);
        gNonlCalc->SetPoint(index, tmp_E, tmp_pred);

        index++;
    }

    TGraphErrors* gNonlNom  = new TGraphErrors();
    gNonlNom->SetName("gNonlNom");
    double nom_kA = 0.96;
    double nom_kB = 6.5e-3;
    double nom_kC = 1;
    double nom_es = 3350./2.22;
    electronQuench::setkA(nom_kA);
    electronQuench::setBirk1(nom_kB);
    electronCerenkov::setkC(nom_kC);
    electronCerenkov::setEnergyScale(nom_es);
    for (int i=0; i<m_nData; i++) {
        gammaData_array[i]->calcGammaResponse();
        gNonlNom->SetPoint(i, gammaData_array[i]->GetEtrue(), gammaData_array[i]->GetNonlPred());
    }

    gNonlData->SetMarkerStyle(20);
    gNonlData->SetMarkerColor(kBlue+1);
    gNonlData->SetLineColor(kBlue+1);
    gNonlData->SetLineWidth(2);
    gNonlData->SetMarkerSize(1.2);
    gNonlCalc->SetMarkerStyle(21);
    gNonlCalc->SetMarkerColor(kRed+1);
    gNonlCalc->SetMarkerSize(1.2);
    gNonlCalc->SetLineColor(kRed+1);
    gNonlCalc->SetLineWidth(2);
    gNonlNom->SetMarkerStyle(21);
    gNonlNom->SetMarkerColor(kViolet+1);
    gNonlNom->SetMarkerSize(1.2);
    gNonlNom->SetLineColor(kViolet+1);
    gNonlNom->SetLineWidth(2);

    TCanvas* c1 = new TCanvas("Nonlinearity", "Nonlinearity");
    c1->cd(); c1->SetGrid();
    gNonlData->SetTitle("Nonlinearity Fitting; Etrue/MeV; Nonlinearity");
    //gNonlData->GetYaxis()->SetRangeUser(0.01,0.045);
    gNonlData->Draw("APL");
    gNonlCalc->Draw("P SAME");
    gNonlNom->Draw("L SAME");
    TLegend* led = new TLegend();
    led->SetFillColor(kWhite);
    led->AddEntry(gNonlData, "data", "PL");
    led->AddEntry(gNonlCalc, "calc", "PL");
    led->AddEntry(gNonlNom, "nominal", "PL");
    led->Draw("SAME");

    c1->SaveAs("GamNLFit.root");    
    
}
