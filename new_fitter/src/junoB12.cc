#include "junoB12.hh"
#include "electronResponse.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "gammaResponse.hh"

#include <iostream>

#include <TFile.h>
#include <TTree.h>

using namespace std;

junoB12::junoB12() 
{
}

junoB12::~junoB12()
{

    delete gamma4440;
    delete gamma3215;
}


void junoB12::Initialize()
{
    gamma4440 = new gammaResponse("gamma4440", 100, 6000, 7000);
    gamma3215 = new gammaResponse("gamma3215", 100, 4000, 5000);
    gamma4440->LoadData();
    gamma3215->LoadData();


    m_loadData = false;
    m_loadTheo = false;

    // binning method :
    m_binWidth = (15 - 0.) / 1500.;
    for (int i=0; i<1500; i++ ){
        m_binCenter[i] = m_binWidth/2. + m_binWidth * i;
    }


    m_eTruGam[0][0] = 0;
    m_eTruGam[0][1] = 0;
    m_eTruGam[1][0] = 4.44;
    m_eTruGam[1][1] = 0;
    m_eTruGam[2][0] = 4.44;
    m_eTruGam[2][1] = 3.215;

    LoadDataSpec();
    LoadTheoSpec();

}


void junoB12::LoadDataSpec()
{
    double scale = 3134.078/2.223;
    TH1D* sigH = new TH1D("B12_data", "", 100, 0, 15);

    TFile* ff = new TFile("./data/spectrum/data/B12_data_G4_J19.root", "read");
    //if(!ff) cout << "No such B12 data file " <<  endl;
    TTree* tt = (TTree*)ff->Get("B12");
    int m_totpe;
    tt->SetBranchAddress("totpe", &m_totpe);
    for(int i=0; i<tt->GetEntries(); i++) {
        tt->GetEntry(i);
        double tmp_Evis = m_totpe / scale;
        sigH->Fill(tmp_Evis);
    }
    for (int i=0; i<100; i++) {
		double content = sigH->GetBinContent(i+1);
		double error   = sigH->GetBinError  (i+1);
		m_eData   [i] = content;
		m_eDataErr[i] = error * 1.5;
    }
    
	delete sigH;
    
    delete tt;
    delete ff;
    m_loadData = true;
    cout << " >>> Load Data Spectrum for B12 <<< " << endl;
}


void junoB12::LoadTheoSpec()
{
    //cout <<"Etrue = " << gamma4440->GetEtrue() << endl;

    TFile* ff = new TFile("./data/spectrum/theo/B12_sim.root");
    if (!ff) cout << "No such B12 theo file !" << endl;
    TTree* tt = (TTree*)ff->Get("T");
    int eGammaStr[10];
    double branchRatio;
    double weight;
    tt->SetBranchAddress("BR"       ,&branchRatio);
    tt->SetBranchAddress("photonName", eGammaStr);
    for (int branchIdx=0;branchIdx!=tt->GetEntries();branchIdx++){ 
        tt->GetEntry(branchIdx);

        // beta
        TH1F* electronHist = (TH1F*)ff->Get(Form("hh%d", branchIdx));
        for (int binIdx=0;binIdx!=1500;binIdx++)
        {
            weight = electronHist->Interpolate(m_binCenter[binIdx]);
            m_eTru[branchIdx][binIdx] = branchRatio * weight;
        }
        delete electronHist;
    }

    delete tt;
    delete ff;

    m_loadTheo = true;
    cout << " >>> Load Theo Spectrum for B12 <<< " << endl;
}

int junoB12::ApplyNL()
{
    cout << gamma4440->GetEtrue() << endl;
    return 1.0;
}

void junoB12::ApplyScintillatorNL()
{
    if (not m_loadData) LoadDataSpec();
    if (not m_loadTheo) LoadTheoSpec();

    for (int i=0; i<1500; i++) {
        m_eVis[i] = 0;
    }

    int newBin;
    int newBinLow, newBinHig;
    double bias;

    gamma4440->Prediction();
    gamma3215->Prediction();

    // gamma
    m_eVisGam[0][0] = 0;
    m_eVisGam[0][1] = 0;
    m_eVisGam[1][0] = m_eTruGam[1][0] * gamma4440->GetNonlCalc();
    m_eVisGam[1][1] = 0;
    m_eVisGam[2][0] = m_eTruGam[2][0] * gamma4440->GetNonlCalc();
    m_eVisGam[2][1] = m_eTruGam[2][1] * gamma3215->GetNonlCalc();
    //cout << ">>>>>>>>> Gamma Evis <<<<<<<<< " << endl;
    //cout << m_eVisGam[1][0] << " " << m_eVisGam[2][0] <<" " << m_eVisGam[2][1] << endl;

    // beta
    for (int i=0; i<1500; i++) {
        double eTru = m_binCenter[i];
        //double eVis = eTru * electronResponse::getElecNonl(eTru);
        double eVis = (electronQuench::ScintillatorPE(eTru) + electronCerenkov::getCerPE(eTru))/electronQuench::getEnergyScale(); 
        for (int j=0; j<2; j++) {
            eVis += m_eVisGam[j][0] + m_eVisGam[j][1];

            newBinLow = int((eVis-0)/m_binWidth);
            newBinHig = newBinLow + 1;
            bias      = (eVis - 0 - newBinLow * m_binWidth) / m_binWidth;
            //cout << newBinLow << " " << newBinHig << " " << bias << endl;

            if (newBinLow < 1500) m_eVis[newBinLow] += (1-bias) * m_eTru[j][i];
            if (newBinHig < 1500) m_eVis[newBinHig] += bias * m_eTru[j][i];
        }
    }
}



void junoB12::Normalize()
{
    // Normalize spectrum for data and pred
	int   rebin = 1500/100;
	double binWidthData = m_binWidth * rebin;
	double nTheo = 0;
	double nData = 0;
	for (int i = 0; i < 100; i++)
	{
		m_eTheo[i] = 0;
        for (int j = 0; j < rebin; j++){
			m_eTheo[i] += m_eVis[i*rebin+j];
        } 
            
		if(i*binWidthData>3 && i*binWidthData<12)   // fitting range [3MeV, 12MeV]
		{
			nTheo += m_eTheo[i];
			nData += m_eData[i];
		}
	}
    double scale = 1;
    if( nTheo!=0 ) { scale = nData/nTheo; }
	for (int i = 0; i < 100; i++)
    {
		m_eTheo[i] *= scale;
        
    }
	for (int i = 0; i < 100; i++)
	{
		m_eVis   [i] *= scale;
	}

}


double junoB12::GetChi2()
{
    ApplyScintillatorNL();
    Normalize();

    double chi2 = 0;
    int rebin = 1500 / 100;
    double binWidthData = m_binWidth * rebin;
    int m_nData = 0;
    for(int i=30; i < 90; i++) {
        if(i*binWidthData<3 or binWidthData*i>12) continue;
        if( m_eDataErr[i]!=0 ) {
            chi2 += pow( (m_eData[i] - m_eTheo[i])/m_eDataErr[i], 2); 
            m_nData++;
        }
    }
    cout << "B12 chi2: " << chi2 << " with nData : " << m_nData << endl;
	//if(nDoF>0) chi2 /= double(m_nData - nDoF);
	return chi2;
}


void junoB12::Plot()
{
    TH1D* hData = new TH1D("hData", "", 100, 0, 15);
    TH1D* hTheo = new TH1D("hTheo", "", 100, 0, 15);
    TH1D* hRela = new TH1D("hRela", "", 100, 0, 15);

    hData->SetStats(0);
    hData->SetLineColor(kBlue+1);
    hData->SetLineWidth(2);
    hData->SetMarkerSize(0.8);
    hData->SetMarkerStyle(20);
    hData->SetMarkerColor(kBlue+1);
    hTheo->SetStats(0);
    hTheo->SetLineColor(kRed+1);
    hTheo->SetLineWidth(2);
    hTheo->SetMarkerSize(0.8);
    hTheo->SetMarkerStyle(20);
    hTheo->SetMarkerColor(kRed+1);
    hRela->SetStats(0);
    hRela->SetLineColor(kPink+2);
    hRela->SetLineWidth(2);
    hRela->SetMarkerSize(0.8);
    hRela->SetMarkerStyle(20);

    for(int i=0; i<100; i++) {
        hData->SetBinContent(i+1, m_eData[i]);
        hTheo->SetBinContent(i+1, m_eTheo[i]);
        if(m_eTheo[i]!=0) {
            hRela->SetBinContent(i+1, m_eData[i]/m_eTheo[i]);
            hRela->SetBinError(i+1, m_eDataErr[i]/m_eTheo[i]);
        } else {
            hRela->SetBinContent(i+1, 0);
            hRela->SetBinError(i+1, 0);
        }
    }


    TFile* out = new TFile("spectrum.root", "recreate");
    hData->Write();
    hTheo->Write();
    out->Close();


}

