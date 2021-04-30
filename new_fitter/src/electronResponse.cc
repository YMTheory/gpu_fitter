#include "electronResponse.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "junoParameters.hh"

#include <TFile.h>
#include <TGraph.h>

double electronResponse::m_SimEtrue[m_nSimData];
double electronResponse::m_SimNonl[m_nSimData];
double electronResponse::m_scale = 1503.664;
bool electronResponse::m_loadSimFile = false;
bool electronResponse::m_doFit = false;
bool electronResponse::m_loadResol = false;

TGraph* electronResponse::gSimData;
TF1* electronResponse::fElecNonl;
double electronResponse::m_p0 = 1.025;
double electronResponse::m_p1 = 0.1122;
double electronResponse::m_p2 = 1.394;
double electronResponse::m_p3 = 5.55e-4;
double electronResponse::m_ra = -2.17203e+00;
double electronResponse::m_rb = 1.31498e+03;
double electronResponse::m_rc = 1.60508e+02;

TGraphErrors* electronResponse::gMinElecNonl;
TGraphErrors* electronResponse::gElecResol;

double gElecResolFunc(double* x, double* p) {
    double E = x[0];
    double p0 = p[0];
    double p1 = p[1];
    double p2 = p[2];

    double sigma2 = p0 + p1*E + p2*E*E;
    if (sigma2<0)
        return 0;
    else 
        return TMath::Sqrt(sigma2);
}

TF1* electronResponse::fElecResol = new TF1("fElecNonl", gElecResolFunc, 0, 8, 3);

double electronResponse::getElecNonl(double Etrue)
{
    return electronQuench::ScintillatorNL(Etrue) + electronCerenkov::getCerenkovPE(Etrue);
}


void electronResponse::loadSimElecNonl()
{
    gSimData = new TGraph();
    ifstream in; 
    in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt");
    string line;
    double Etrue, nonl, totpe, totpe_sigma;
    int index = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> Etrue >> totpe >> totpe_sigma;
        m_SimEtrue[index] = Etrue/1000.;
        if (Etrue == 0)
            m_SimNonl[index] = 0;
        else{
            nonl = totpe/Etrue/m_scale*1000.;
            m_SimNonl[index] = nonl;
        }
        gSimData->SetPoint(index, m_SimEtrue[index], m_SimNonl[index]);
        index++;
    }
    in.close();

    m_loadSimFile = true;
    return;
}

void electronResponse::loadMinSimElecNonl()
{
    gMinElecNonl = new TGraphErrors();
    ifstream in;
    in.open("./data/electron/electron_response.txt");
    if (!in) std::cout << " >>> No electron response file!!! <<< " << std::endl;
    string line;
    double Etrue, totpe, totpe_err, sigma, sigma_err;
    int index = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> Etrue >> totpe >> totpe_err >> sigma >> sigma_err;
        gMinElecNonl->SetPoint(index, Etrue, totpe/Etrue/m_scale);
        gMinElecNonl->SetPointError(index, 0, 0.001);
        index++;
    }

    in.close();
}






void electronResponse::FuncConstruct()
{
    fElecNonl = new TF1("fElecNonl", "([0]+[3]*x)/(1+[1]*TMath::Exp(-[2]*x))", 0.01, 8);
    SetParameters();
    cout << "Empirical ElecNonl Function has been constructed..." << endl;
}

void electronResponse::SetParameters()
{
    //fElecNonl->SetParameter(0, m_p0);
    //fElecNonl->SetParameter(1, m_p1);
    //fElecNonl->SetParameter(2, m_p2);
    //fElecNonl->SetParameter(3, m_p3);

    fElecResol->SetParameter(0, m_ra);
    fElecResol->SetParameter(1, m_rb);
    fElecResol->SetParameter(2, m_rc);
}


void electronResponse::Plot()
{
    cout << ">>> Draw Analytic Electron Nonlinearity Curve <<< "<< endl;
    if (not m_loadSimFile) loadSimElecNonl();
    //if (not m_doFit) EmpiricalFit();

    gSimData->SetName("elec");
    gSimData->SetLineColor(kBlue+1);
    gSimData->SetLineWidth(2);

    fElecNonl->SetLineColor(kOrange+1);
    fElecNonl->SetLineWidth(2);

    TFile* out = new TFile("simElecNonl.root", "recreate");
    gSimData->Write();
    fElecNonl->Write();
    out->Close();
}


void electronResponse::EmpiricalFit()
{
    if(not m_loadSimFile) loadSimElecNonl();
    fElecNonl->SetParameters(1.025, 0.1122, 1.394, 5.55e-4);
    gSimData->Fit(fElecNonl, "R");

    m_doFit = true;
}


void electronResponse::FitPlot()
{
    TGraph* gNom = new TGraph();
    TGraph* gFit = new TGraph();
    gNom->SetName("nom");
    gFit->SetName("fit");


    electronQuench::setkA(1);
    electronQuench::setEnergyScale(3300.371/2.223);
    electronCerenkov::setkC(1);
    electronCerenkov::setEnergyScale(3300.371/2.223);
    for(int i=0; i<500; i++) {
        double Etrue = 10./500 * (i+1);
        double nonl_nom = electronResponse::getElecNonl(Etrue);
        gNom->SetPoint(i, Etrue, nonl_nom);
    }

    gNom->SetLineWidth(2);
    gNom->SetLineColor(kBlue+1);


    electronQuench::setkA(1.01959e+00);
    electronQuench::setEnergyScale(3300.371/2.223);
    electronCerenkov::setkC(7.09311e-01);
    electronCerenkov::setEnergyScale(3300.371/2.223);
    for(int i=0; i<500; i++) {
        double Etrue = 10./500 * (i+1);
        double nonl_fit = electronResponse::getElecNonl(Etrue);
        gFit->SetPoint(i, Etrue, nonl_fit);
    }

    gFit->SetLineWidth(2);
    gFit->SetLineColor(kGreen+1);


    TFile* out = new TFile("ElecNonlCompare.root", "recreate");
    gNom->Write();
    gFit->Write();
    out->Close();


}

void electronResponse::loadElecResol()
{
    gElecResol = new TGraphErrors();
    ifstream in; in.open(junoParameters::electronResol_File.c_str());
    string line;
    double tmp_E, tmp_mu, tmp_muerr, tmp_sigma, tmp_sigmaerr, tmp_resol, tmp_resolerr;
    int index = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmp_E >> tmp_mu >> tmp_muerr >> tmp_sigma >> tmp_sigmaerr >> tmp_resol >> tmp_resolerr;
        gElecResol->SetPoint(index, tmp_E, tmp_sigma);
        index++;
    }
    in.close();

    m_loadResol = true;
}



