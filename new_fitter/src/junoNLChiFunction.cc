#include "junoNLChiFunction.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"
#include "junoParameters.hh"
#include "junoSpectrum.hh"
#include "junoB12.hh"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

gammaData* junoNLChiFunction::Cs137Data;
gammaData* junoNLChiFunction::Mn54Data;
gammaData* junoNLChiFunction::nHData;
gammaData* junoNLChiFunction::K40Data;
gammaData* junoNLChiFunction::Co60Data;
gammaData* junoNLChiFunction::Tl208Data;
gammaData* junoNLChiFunction::nC12Data;
gammaData* junoNLChiFunction::O16Data;
gammaData* junoNLChiFunction::nFe56Data;
gammaData* junoNLChiFunction::gamma511Data;
gammaData* junoNLChiFunction::gamma4440Data;

gammaResponse* junoNLChiFunction::Cs137;
gammaResponse* junoNLChiFunction::Mn54;
gammaResponse* junoNLChiFunction::K40;
gammaResponse* junoNLChiFunction::Ge68;
gammaResponse* junoNLChiFunction::Co60;
gammaResponse* junoNLChiFunction::nH;
gammaResponse* junoNLChiFunction::AmBe;
gammaResponse* junoNLChiFunction::nC12;
gammaResponse* junoNLChiFunction::AmC;

junoSpectrum* junoNLChiFunction::junoB12data;

junoB12_simplified* junoNLChiFunction::b12data;

double junoNLChiFunction::m_chi2    = 0.;
double junoNLChiFunction::m_chi2Min = 100000;
double junoNLChiFunction::m_chi2B12 = 0;
double junoNLChiFunction::m_chi2C11 = 0;
double junoNLChiFunction::m_chi2C10 = 0;
double junoNLChiFunction::m_chi2Gam = 0;

int junoNLChiFunction::m_nParameter = 3;
double junoNLChiFunction::m_bestFit[20] = {0.};
double junoNLChiFunction::m_bestFitError[20] = {0.};

bool junoNLChiFunction::m_DoFit = false;

int junoNLChiFunction::m_nGam;
int junoNLChiFunction::m_nData;
std::string junoNLChiFunction::source_name[20];
//gammaData* junoNLChiFunction::gammaData_array[20];
gammaResponse* junoNLChiFunction::gammaData_array[10];


bool junoNLChiFunction::m_doGamFit = junoParameters::fitGammaSources;
bool junoNLChiFunction::m_doB12Fit = junoParameters::fitB12;

string junoNLChiFunction::m_nonlMode;

junoNLChiFunction::junoNLChiFunction() {

    cout << ">>>>>>>>>>>> Fitting Mode ==> " << junoParameters::m_calcOption << endl;

    m_nData = 0;
    m_nGam  = 0;

    if (m_doGamFit) {

        //gamma511Data = new gammaData("511keV", 700, 1100, 100);
        //source_name[m_nData] = "gamma511"; 
        //gammaData_array[m_nData] = gamma511Data;
        //m_nData++;
        //m_nGam++;

        //Cs137Data = new gammaData("Cs137", 700, 1100, 100);
        //source_name[m_nData] = "Cs137"; 
        //gammaData_array[m_nData] = Cs137Data;
        //m_nData++;
        //m_nGam++;

        //Mn54Data  = new gammaData("Mn54", 900, 1300, 100);
        //source_name[m_nData] = "Mn54"; 
        //gammaData_array[m_nData] = Mn54Data;
        //m_nData++;
        //m_nGam++;

        //K40Data  = new gammaData("K40", 900, 1300, 100);
        //source_name[m_nData] = "K40"; 
        //gammaData_array[m_nData] = K40Data;
        //m_nData++;
        //m_nGam++;

        //nHData  = new gammaData("nH", 900, 1300, 100);
        //source_name[m_nData] = "nH"; 
        //gammaData_array[m_nData] = nHData;
        //m_nData++;
        //m_nGam++;

        ////Co60Data  = new gammaData("Co60", 900, 1300, 100);
        ////source_name[m_nData] = "Co60"; 
        ////gammaData_array[m_nData] = Co60Data;
        ////m_nData++;

        //Tl208Data  = new gammaData("Tl208", 900, 1300, 100);
        //source_name[m_nData] = "Tl208"; 
        //gammaData_array[m_nData] = Tl208Data;
        //m_nData++;
        //m_nGam++;

        //nC12Data  = new gammaData("nC12", 900, 1300, 100);
        //source_name[m_nData] = "nC12"; 
        //gammaData_array[m_nData] = nC12Data;
        //m_nData++;
        //m_nGam++;

        //O16Data  = new gammaData("O16", 900, 1300, 100);
        //source_name[m_nData] = "O16"; 
        //gammaData_array[m_nData] = O16Data;
        //m_nData++;
        //m_nGam++;

        //nFe56Data  = new gammaData("nFe56", 900, 1300, 100);
        //source_name[m_nData] = "nFe56"; 
        //gammaData_array[m_nData] = nFe56Data;
        //m_nData++;
        //m_nGam++;

        Cs137 = new gammaResponse("Cs137", 100, 600, 1000);
        source_name[m_nData] = "Cs137";
        gammaData_array[m_nData] = Cs137;
        m_nData++;
        m_nGam++;

        Mn54 = new gammaResponse("Mn54", 100, 900, 1300);
        source_name[m_nData] = "Mn54";
        gammaData_array[m_nData] = Mn54;
        m_nData++;
        m_nGam++;

        Ge68 = new gammaResponse("Ge68", 100, 1100, 1500);
        source_name[m_nData] = "Ge68";
        gammaData_array[m_nData] = Ge68;
        m_nData++;
        m_nGam++;

        K40 = new gammaResponse("K40", 100, 1750, 2250);
        source_name[m_nData] = "K40";
        gammaData_array[m_nData] = K40;
        m_nData++;
        m_nGam++;

        nH = new gammaResponse("nH", 100, 2800, 3500);
        source_name[m_nData] = "nH";
        gammaData_array[m_nData] = nH;
        m_nData++;
        m_nGam++;

        Co60 = new gammaResponse("Co60", 100, 3200, 3700);
        source_name[m_nData] = "Co60";
        gammaData_array[m_nData] = Co60;
        m_nData++;
        m_nGam++;

        AmBe = new gammaResponse("AmBe", 100, 6000, 6800);
        source_name[m_nData] = "AmBe";
        gammaData_array[m_nData] = AmBe;
        m_nData++;
        m_nGam++;

        nC12 = new gammaResponse("nC12", 100, 6700, 7600);
        source_name[m_nData] = "nC12";
        gammaData_array[m_nData] = nC12;
        m_nData++;
        m_nGam++;

        AmC = new gammaResponse("AmC", 100, 8400, 9400);
        source_name[m_nData] = "AmC";
        gammaData_array[m_nData] = AmC;
        m_nData++;
        m_nGam++;


    }

    // Nonlinearity mode
    m_nonlMode = junoParameters::m_nonlMode;
    cout << "Nonlinearity formula form " << m_nonlMode << endl;

    if (m_doB12Fit) {
        //junoB12data = new junoSpectrum(1500, 100, 3, 2,
        //                     0, 15, 1, 14, m_nonlMode, "B12");
        b12data = new junoB12_simplified(100, 4500, 17500);
        //b12data = new junoB12();
    }

    electronResponse::FuncConstruct();
    electronResponse::loadElecResol();
}

junoNLChiFunction::~junoNLChiFunction() {
    if (m_doGamFit) {
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
    if(m_doB12Fit)
        //delete junoB12data;
        delete b12data;
}

void junoNLChiFunction::LoadData()
{
    if (m_doGamFit) {
        for (int i=0; i<m_nData; i++) {
            gammaData_array[i]->LoadData();
        }
    }
    if (m_doB12Fit) {
        //junoB12data->LoadData();
        b12data->Initialize();
    }
}

double junoNLChiFunction::GetChi2( double maxChi2 )
{
    double chi2 = 0;

    if (m_doGamFit ) {
        for(int iSource=0; iSource<m_nData; iSource++) {
            chi2 += gammaData_array[iSource]->GetChi2();
        }
    }
    if(m_doB12Fit)
        //chi2 += junoB12data->GetChi2();
        chi2 += b12data->GetChi2();    

    cout << "current total chi2 = " << chi2 << endl;
    return chi2;
}

void junoNLChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}

void junoNLChiFunction::SetParameters(double *par)
{
    //if (junoParameters::scintillatorParameterization == kSimulation) { 
    //    if (m_nonlMode == "histogram") {
    //        electronQuench::setkA               (par[0]);
    //        electronQuench::setBirk1            (par[1]);
    //        electronCerenkov::setkC             (par[2]);
    //        electronCerenkov::setEnergyScale    (par[3]);
    //        junoParameters::m_nuGamma = par[4];
    //    }

    //    if (m_nonlMode == "analytic") {
    //        electronResponse::setp0(par[0]);
    //        electronResponse::setp1(par[1]);
    //        electronResponse::setp2(par[2]);
    //        electronResponse::setp3(par[3]);
    //        electronResponse::SetParameters();
    //    }
    //}

    if (junoParameters::scintillatorParameterization == kSimulation and m_doGamFit) {
        electronQuench::setkA               (par[0]);
        electronQuench::setBirk1            (par[1]);
        electronCerenkov::setkC             (par[2]);
        electronQuench::setEnergyScale      (par[3]);
        electronCerenkov::setEnergyScale    (par[3]);
        junoParameters::m_nuGamma           = par[4];
        Cs137->SetAmp                       (par[5]);
        Mn54->SetAmp                        (par[6]);
        Ge68->SetAmp(par[7]);
        K40->SetAmp(par[8]);
        nH->SetAmp(par[9]);
        Co60->SetAmp(par[10]);
        AmBe->SetAmp(par[11]);
        nC12->SetAmp(par[12]);
        AmC->SetAmp(par[13]);
        electronResponse::setra(par[13]);
        electronResponse::setrb(par[15]);
        electronResponse::setrc(par[16]);

        electronResponse::SetParameters();
    }
    if (junoParameters::scintillatorParameterization == kSimulation and !m_doGamFit) {
        electronQuench::setkA               (par[0]);
        electronQuench::setBirk1            (par[1]);
        electronCerenkov::setkC             (par[2]);
        electronQuench::setEnergyScale      (par[3]);
        electronCerenkov::setEnergyScale    (par[3]);
        junoParameters::m_nuGamma           = par[4];
        electronResponse::setra(par[4]);
        electronResponse::setrb(par[5]);
        electronResponse::setrc(par[6]);

        electronResponse::SetParameters();
    }
}


double junoNLChiFunction::GetChiSquare(double maxChi2)
{
    junoNLMinuit = new TMinuit();
    junoNLMinuit->SetFCN(ChisqFCN);
    junoNLMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    junoNLMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    //if (junoParameters::scintillatorParameterization == kSimulation) {
    //    if (m_nonlMode == "histogram") {
    //        junoNLMinuit->mnparm(iPar, "kA", 1.00, 0.001, 0., 0., ierrflag);              iPar++;
    //        junoNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.1e-3, 7.5e-3, ierrflag);     iPar++;
    //        junoNLMinuit->mnparm(iPar, "kC", 1.0, 0.001, 0., 0., ierrflag);               iPar++;
    //        junoNLMinuit->mnparm(iPar, "energyScale", 3134.078/2.223, 1, 0, 0, ierrflag); iPar++;
    //        junoNLMinuit->mnparm(iPar, "nuGamma", 0.01, 0.0001, 0., 1, ierrflag);         iPar++;
    //    }

    //    if (m_nonlMode == "analytic") {
    //        junoNLMinuit->mnparm(iPar, "p0", 1.025, 0.001, 0.9, 1.1, ierrflag);    iPar++;
    //        junoNLMinuit->mnparm(iPar, "p1", 0.1122,0.0001,  0, 0.2, ierrflag);    iPar++;
    //        junoNLMinuit->mnparm(iPar, "p2", 1.394, 0.001, 1.1, 4.0, ierrflag);    iPar++;
    //        junoNLMinuit->mnparm(iPar, "p3", 5.55e-4, 1e-5, 1e-5, 1e-2, ierrflag); iPar++;
    //    }
    //}

    if (junoParameters::scintillatorParameterization == kSimulation and m_doGamFit) {
        junoNLMinuit->mnparm(iPar, "kA", 1.00, 0.001, 0.9, 1.1, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-4, 5e-3 , 7.5e-3, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "kC", 1.00, 0.001, 0.0, 1.5, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "energyScale", 3134.078/2.223, 1, 1000, 2700, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "nuGamma", 0.0, 0.0001, 0., 1, ierrflag);         iPar++;
        junoNLMinuit->mnparm(iPar, "Cs137Amp", 250, 1, 200, 300, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "Mn54Amp", 240, 1, 200, 300, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "Ge68Amp", 220, 1, 180, 260, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "K40Amp", 180, 1, 140, 220, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "nHAmp", 200, 1, 160, 240, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "Co60Amp", 160, 1, 120, 200, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "AmBeAmp", 180, 1, 140, 220, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "nC12Amp", 180, 1, 140, 220, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "AmCAmp", 180, 1, 140, 220, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "ra", 0, 0.01, -20, 20, ierrflag ); iPar++;
        junoNLMinuit->mnparm(iPar, "rb", 1315, 1, 1250, 1450, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "rc", 160, 1, 100, 220, ierrflag); iPar++;
        
    }

    if (junoParameters::scintillatorParameterization == kSimulation and !m_doGamFit) {
        junoNLMinuit->mnparm(iPar, "kA", 1.00, 0.001, 0.8, 2.2, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "kC", 1.00, 0.001, -3, 3.2, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "energyScale", 3134.078/2.223, 1, 1000, 2700, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "nuGamma", 0.0, 0.0001, 0., 1, ierrflag);         iPar++;
        junoNLMinuit->mnparm(iPar, "ra", 0, 0.01, -20, 20, ierrflag ); iPar++;
        junoNLMinuit->mnparm(iPar, "rb", 1315, 1, 1250, 1450, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "rc", 160, 1, 100, 920, ierrflag); iPar++;
            
        //junoNLMinuit->FixParameter(4);
        //junoNLMinuit->FixParameter(5);
        //junoNLMinuit->FixParameter(6);
        //
    }


    //junoNLMinuit->FixParameter(0);
    //junoNLMinuit->FixParameter(1);
    //junoNLMinuit->FixParameter(2);
    junoNLMinuit->FixParameter(3);

    // Minimization strategy
    junoNLMinuit->SetErrorDef(1);
    arglist[0]=2;
    junoNLMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 5000; //maxCalls
    arglist[1] = 0.01; // tolerance
    junoNLMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    junoNLMinuit->fCstatu.Data();

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    junoNLMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    m_nParameter = junoNLMinuit->GetNumPars();
	for(int i=0; i<m_nParameter; i++)
	{
	    junoNLMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
		//cout<<"curvalue: "<<curvalue<<"	curerror: "<<curerror<<endl;
	}

    m_DoFit = true;

    if (m_doGamFit) {
        Cs137->SaveHist();
        Mn54->SaveHist();
        Ge68->SaveHist();
        K40->SaveHist();
        nH->SaveHist();
        Co60->SaveHist();
        AmBe->SaveHist();
        nC12->SaveHist();
        AmC->SaveHist();
    }
    if (m_doB12Fit)
        //m_nData += junoB12data->getNData();

        cout << " ====================== " << endl;
    cout << "    minChi2: " << min << " with nData = " << m_nData << " and nPar = " << m_nParameter << endl;
    cout << " ====================== " << endl;
    delete junoNLMinuit;
    return min;
}



void junoNLChiFunction::Plot()
{
    //electronResponse::Plot();
    if(m_doGamFit)
    {
        GammaPEPlot();
        GammaPlot();
    }
    if(m_doB12Fit)
        //junoB12data->Plot();
        b12data->Plot();
}


void junoNLChiFunction::GammaPlot()
{
    //cout << " >>> Draw Gamma NL Fitting Results <<< " << endl;

    if (not m_DoFit) {
        cout << "Fitting has not been finished ...";
        return;
    }

    std::cout << " >>>>>>>>>>>> GammaNL Outputs <<<<<<<<<<<< " << std::endl;

    TGraphErrors* gNonlData = new TGraphErrors();
    TGraphErrors* gNonlCalc = new TGraphErrors();
    gNonlData->SetName("gNonlData");
    gNonlCalc->SetName("gNonlCalc");

    double fit_pars[4] = {m_bestFit[0], m_bestFit[1], m_bestFit[2], m_bestFit[3]};
    SetParameters(fit_pars);
    int index = 0;
    for(int iData=0; iData<m_nGam; iData++) {
        std::string source = source_name[iData];
        gammaResponse* tmpGammaResponse = gammaData_array[iData];
        tmpGammaResponse->calcGamResponse();
        double tmp_E       = tmpGammaResponse->GetEtrue();
        double tmp_pred    = tmpGammaResponse->GetNonlCalc();
        double tmp_data    = tmpGammaResponse->GetNonlData();
        double tmp_dataErr = tmpGammaResponse->GetNonlErr(); 
        gNonlData->SetPoint(index, tmp_E, tmp_data);
        gNonlData->SetPointError(index, 0, tmp_dataErr);
        gNonlCalc->SetPoint(index, tmp_E, tmp_pred);

        //tmpGammaData->SaveHist();

        index++;
    }

    gNonlData->SetMarkerStyle(20);
    gNonlData->SetMarkerColor(kBlue+1);
    gNonlData->SetLineColor(kBlue+1);
    gNonlData->SetLineWidth(2);
    gNonlData->SetMarkerSize(1.0);
    gNonlCalc->SetMarkerStyle(21);
    gNonlCalc->SetMarkerColor(kRed+1);
    gNonlCalc->SetMarkerSize(1.0);
    gNonlCalc->SetLineColor(kRed+1);
    gNonlCalc->SetLineWidth(2);

    //cout << "\n";
    //cout << " >>> Nominal Outputs <<< " << endl;
    //double nom_pars[4] = {1.0, 1.0, 3300.371/2.223, 0};
    //SetParameters(nom_pars);
    //TGraphErrors* gNom = new TGraphErrors();
    //index = 0;
    //for(int iData=0; iData<m_nGam; iData++) {
    //    std::string source = source_name[iData];
    //    gammaData* tmpGammaData = gammaData_array[iData];
    //    tmpGammaData->calcGammaResponse();
    //    double tmp_E       = tmpGammaData->GetEtrue();
    //    double tmp_pred    = tmpGammaData->GetNonlPred();
    //    double tmp_pred1   = tmpGammaData->GetNonlPred1();
    //    double tmp_data    = tmpGammaData->GetNonlData();
    //    double tmp_dataErr = tmpGammaData->GetNonlDataErr(); 
    //    double tmp_pedata  = tmpGammaData->GetPEData();
    //    double tmp_pecalc  = tmpGammaData->GetPECalc();
    //    double tmp_cerpecalc = tmpGammaData->GetCerPECalc();
    //    double tmp_sctpecalc = tmpGammaData->GetSctPECalc();
    //    double tmp_sctpedata = tmpGammaData->GetSctPEData();
    //    double tmp_cerpedate = tmpGammaData->GetCerPEData();
    //    //tmpGammaData->calcGammaNPE();
    //    cout << "\n";
    //    cout << source_name[iData] << " " << tmp_E << endl;
    //    cout << tmp_pedata << " " << tmp_pecalc << endl;
    //    cout << tmp_sctpedata << " " << tmp_sctpecalc << endl;
    //    cout << tmp_cerpedate << " " << tmp_cerpecalc << endl;
    //    cout << tmp_data << " " << tmp_pred << endl;
    //    cout << "\n";
    //    gNom->SetPoint(index, tmp_E, tmp_pred1);

    //    index++;
    //}
    //gNom->SetMarkerStyle(21);
    //gNom->SetMarkerColor(kOrange+1);
    //gNom->SetMarkerSize(1.0);
    //gNom->SetLineColor(kOrange+1);
    //gNom->SetLineWidth(2);


    TCanvas* c1 = new TCanvas("Nonlinearity", "Nonlinearity");
    c1->cd(); c1->SetGrid();
    gNonlData->SetTitle("Nonlinearity Fitting; Etrue/MeV; Nonlinearity");
    gNonlData->GetYaxis()->SetRangeUser(0.90, 1.05);
    gNonlData->Draw("APL");
    gNonlCalc->Draw("P SAME");
    //gNom->Draw("LP SAME");
    TLegend* led = new TLegend();
    led->SetFillColor(kWhite);
    led->AddEntry(gNonlData, "data", "PL");
    led->AddEntry(gNonlCalc, "calc", "PL");
    //led->AddEntry(gNom, "nominal", "PL");
    led->Draw("SAME");

    c1->SaveAs("GamNLFit.root");    
    
}



void junoNLChiFunction::GammaPEPlot()
{
    //std::cout << " >>>>>>>>>>>> GammaPE Outputs <<<<<<<<<<<< " << std::endl;

    //TGraphErrors* gPEData = new TGraphErrors();
    //TGraphErrors* gPECalc = new TGraphErrors();
    //TGraphErrors* gPENomi = new TGraphErrors();
    //TGraphErrors* gCerPEData = new TGraphErrors();
    //TGraphErrors* gCerPECalc = new TGraphErrors();
    //TGraphErrors* gCerPENomi = new TGraphErrors();
    //TGraphErrors* gSctPEData = new TGraphErrors();
    //TGraphErrors* gSctPECalc = new TGraphErrors();
    //TGraphErrors* gSctPENomi = new TGraphErrors();
    //TGraphErrors* gPEDiff = new TGraphErrors();
    //

    //gPEData->SetName("PEData");
    //gPECalc->SetName("PECalc");
    //gPENomi->SetName("PENomi");
    //gCerPEData->SetName("CerPEData");
    //gCerPECalc->SetName("CerPECalc");
    //gCerPENomi->SetName("CerPENomi");
    //gSctPEData->SetName("SctPEData");
    //gSctPECalc->SetName("SctPECalc");
    //gSctPENomi->SetName("SctPENomi");
    //gPEDiff->SetName("PEDiff");

    //std::cout << " >>>>>>> GammaPE Fitting Outputs <<<<<<< " << std::endl;
    ////double fit_pars[4] = {m_bestFit[0], m_bestFit[1], m_bestFit[2], m_bestFit[3]};
    //double fit_pars[4] = {1, 1, 3300.371/2.223, 0};
    //SetParameters(fit_pars);
    //int index = 0;
    //for(int iData=0; iData<m_nGam; iData++) {
    //    std::string source = source_name[iData];
    //    gammaData* tmpGammaData = gammaData_array[iData];
    //    tmpGammaData->calcGammaResponse();
    //    double tmp_E       = tmpGammaData->GetEtrue();
    //    double tmp_pred    = tmpGammaData->GetPECalc();
    //    double tmp_data    = tmpGammaData->GetPEData();
    //    double tmp_cerdata = tmpGammaData->GetCerPEData();
    //    double tmp_cerpred = tmpGammaData->GetCerPECalc();
    //    double tmp_sctdata = tmpGammaData->GetSctPEData();
    //    double tmp_sctpred = tmpGammaData->GetSctPECalc();
    //    double tmp_deltape = (tmp_pred - tmp_data) / tmp_pred;

    //    gPEData->SetPoint(index, tmp_E, tmp_data);
    //    gPECalc->SetPoint(index, tmp_E, tmp_pred);
    //    gSctPEData->SetPoint(index, tmp_E, tmp_sctdata);
    //    gSctPECalc->SetPoint(index, tmp_E, tmp_sctpred);
    //    gCerPEData->SetPoint(index, tmp_E, tmp_cerdata);
    //    gCerPECalc->SetPoint(index, tmp_E, tmp_cerpred);
    //    gPEDiff->SetPoint(index, tmp_E, tmp_deltape);

    //    index++;
    //}

    //gPEData->SetMarkerStyle(20);
    //gPEData->SetMarkerColor(kBlue+1);
    //gPEData->SetLineColor(kBlue+1);
    //gPEData->SetLineWidth(2);
    //gPEData->SetMarkerSize(1.0);
    //gPECalc->SetMarkerStyle(21);
    //gPECalc->SetMarkerColor(kRed+1);
    //gPECalc->SetMarkerSize(1.0);
    //gPECalc->SetLineColor(kRed+1);
    //gPECalc->SetLineWidth(2);
    //gSctPEData->SetMarkerStyle(20);
    //gSctPEData->SetMarkerColor(kBlue+1);
    //gSctPEData->SetLineColor(kBlue+1);
    //gSctPEData->SetLineWidth(2);
    //gSctPEData->SetMarkerSize(1.0);
    //gSctPECalc->SetMarkerStyle(21);
    //gSctPECalc->SetMarkerColor(kRed+1);
    //gSctPECalc->SetMarkerSize(1.0);
    //gSctPECalc->SetLineColor(kRed+1);
    //gSctPECalc->SetLineWidth(2);
    //gCerPEData->SetMarkerStyle(20);
    //gCerPEData->SetMarkerColor(kBlue+1);
    //gCerPEData->SetLineColor(kBlue+1);
    //gCerPEData->SetLineWidth(2);
    //gCerPEData->SetMarkerSize(1.0);
    //gCerPECalc->SetMarkerStyle(21);
    //gCerPECalc->SetMarkerColor(kRed+1);
    //gCerPECalc->SetMarkerSize(1.0);
    //gCerPECalc->SetLineColor(kRed+1);
    //gCerPECalc->SetLineWidth(2);
    //gPEDiff->SetMarkerStyle(20);
    //gPEDiff->SetMarkerColor(kBlue+1);
    //gPEDiff->SetLineColor(kBlue+1);
    //gPEDiff->SetLineWidth(2);
    //gPEDiff->SetMarkerSize(1.0);
    ///*
    //std::cout << " >>>>>>> GammaPE Nominal Outputs <<<<<<< " << std::endl;
    //index = 0;
    //double nom_pars[4] = {1, 1, 3300.371/2.223, 0};
    //SetParameters(nom_pars);
    //for(int iData=0; iData<m_nGam; iData++) {
    //    std::string source = source_name[iData];
    //    gammaData* tmpGammaData = gammaData_array[iData];
    //    tmpGammaData->calcGammaResponse();
    //    double tmp_E       = tmpGammaData->GetEtrue();
    //    double tmp_pred    = tmpGammaData->GetPECalc();

    //    cout << tmp_E << " " << tmp_pred << endl;
    //    gPENomi->SetPoint(index, tmp_E, tmp_pred);

    //    index++;
    //}

    //gPENomi->SetMarkerStyle(21);
    //gPENomi->SetMarkerColor(kOrange+1);
    //gPENomi->SetMarkerSize(1.0);
    //gPENomi->SetLineColor(kOrange+1);
    //gPENomi->SetLineWidth(2);
    //*/

    //TFile* out = new TFile("GamPECheck.root", "recreate");
    //gPEData->Write();
    //gPECalc->Write();
    //gCerPEData->Write();
    //gCerPECalc->Write();
    //gSctPEData->Write();
    //gSctPECalc->Write();
    //gPEDiff->Write();
    //out->Close();

}

