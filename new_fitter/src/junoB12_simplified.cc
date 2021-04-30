#include "junoB12_simplified.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <iostream>

using namespace std;

junoB12_simplified::junoB12_simplified(int nBinsData, double fitMinPE, double fitMaxPE)
{
    m_nBin = 1000;
    m_nBinData = nBinsData;
    m_fitMinPE = fitMinPE;
    m_fitMaxPE = fitMaxPE;
    m_eMin = 0;
    m_eMax = 15;
    m_peMin = 0;
    m_peMax = 20000;

    m_loadData = false;
    m_loadTheo = false;
}



junoB12_simplified::~junoB12_simplified()
{
    //delete gaus;
}


void junoB12_simplified::Initialize()
{
    m_binWidth = (m_eMax - m_eMin) / m_nBin;
    m_peBinWidth = (m_peMax - m_peMin) / m_nBin;
    for (int i=0; i<m_nBin; i++) {
        m_binCenter[i] = m_binWidth/2. + m_binWidth * i;        // energy region
        m_peBinCenter[i] = m_peBinWidth /2. + m_peBinWidth*i;   // p.e. region
    }

    LoadDataSpec();
    LoadTheoSpec();

    gaus = new TF1("gaus", "1/(TMath::Sqrt(2*TMath::Pi())*[1]) * TMath::Exp(-(x-[0])*(x-[0])/2/[1]/[1])", -100, 30000);
}


void junoB12_simplified::LoadDataSpec()
{  // load in totpe definition

    TH1D* sigH = new TH1D("B12_data", "", m_nBinData, m_peMin, m_peMax);

    TFile* ff = new TFile("./data/spectrum/data/B12_data_G4_J19.root", "read");
    //if(!ff) cout << "No such B12 data file " <<  endl;
    TTree* tt = (TTree*)ff->Get("B12");
    int m_totpe;
    tt->SetBranchAddress("totpe", &m_totpe);
    for(int i=0; i<tt->GetEntries(); i++) {
        tt->GetEntry(i);
        //double tmp_Evis = m_totpe / scale;
        //sigH->Fill(tmp_Evis);
        sigH->Fill(m_totpe);
    }
    for (int i=0; i<m_nBinData; i++) {
		double content = sigH->GetBinContent(i+1);
		double error   = sigH->GetBinError  (i+1);
		m_eData   [i] = content;
		m_eDataErr[i] = error;
    }
    
	delete sigH;
    
    delete tt;
    delete ff;
    m_loadData = true;
    cout << " >>> Load Data Spectrum for B12 <<< " << endl;
    
}

void junoB12_simplified::LoadTheoSpec()
{
    // load in energy region
    TH1D* simH = new TH1D("B12_edep", "", m_nBin, m_eMin, m_eMax);

    TFile* ff = new TFile("./data/spectrum/theo/B12_edep_G4_J19.root");
    if (!ff) cout << "No such B12 theo file !" << endl;
    TTree* tt = (TTree*)ff->Get("B12");
    double m_edep;
    tt->SetBranchAddress("edep", &m_edep);
    for(int i=0; i<tt->GetEntries(); i++) {
        tt->GetEntry(i);
        simH->Fill(m_edep);
    }
    for(int i=0; i<m_nBin; i++) {
        m_eTru[i] = simH->GetBinContent(i);  
    }
    
    delete simH;

    delete tt;
    delete ff;
    m_loadTheo = true;
    cout << " >>> Load Theo Spectrum for B12 <<< " << endl;
}


void junoB12_simplified::ApplyScintillatorNL()
{
    // add LS nonlinearity
    if (not m_loadData) LoadDataSpec();
    if (not m_loadTheo) LoadTheoSpec();

    electronResponse::SetParameters();

    for (int i=0; i<m_nBin; i++) {
        m_eVis[i] = 0;
    }

    int newBin;
    int newBinLow, newBinHig;
    double bias;
    
    for (int i=0; i<m_nBin; i++)  {
        double eTru = m_binCenter[i];
        double tmp_pe = electronQuench::ScintillatorPE(eTru) + electronCerenkov::getCerPE(eTru);
        
        // consider resolution:
        double tmp_sigma = electronResponse::fElecResol->Eval(eTru);
        gaus->SetParameter(0, tmp_pe);
        gaus->SetParameter(1, tmp_sigma);
        int minBin = int((tmp_pe-4.5*tmp_sigma)/m_peBinWidth);
        int maxBin = int((tmp_pe+4.5*tmp_sigma)/m_peBinWidth);
        if (minBin < 0) minBin = 0;
        if(maxBin < m_nBin) maxBin = m_nBin;

        for (int j=minBin; j<maxBin; j++) {
            double tmp_center = m_peBinWidth /2 + m_peBinWidth * j;
            double prob = gaus->Eval(tmp_center);

            m_eVis[j] += prob * m_eTru[i];
        }

        // no resolution applied
        //newBinLow = int((tmp_pe-0)/m_peBinWidth);
        //newBinHig = newBinLow + 1;
        //bias      = (tmp_pe - 0 - newBinLow * m_peBinWidth) / m_peBinWidth;

        //if (newBinLow < m_nBin) m_eVis[newBinLow] += (1-bias) * m_eTru[i];
        //if (newBinHig < m_nBin) m_eVis[newBinHig] += (1-bias) * m_eTru[i];

    }

}



void junoB12_simplified::Normalize()
{
    // Normalize spectrum for data and pred
	int   rebin = m_nBin / m_nBinData;
	double binWidthData = m_peBinWidth * rebin;
	double nTheo = 0;
	double nData = 0;
	for (int i = 0; i < m_nBinData; i++)
	{
		m_eTheo[i] = 0;
        for (int j = 0; j < rebin; j++){
			m_eTheo[i] += m_eVis[i*rebin+j];
        } 
		if(i*binWidthData>m_fitMinPE && i*binWidthData<m_fitMaxPE)   // fitting range [3MeV, 12MeV]
		{
			nTheo += m_eTheo[i];
			nData += m_eData[i];
		}

	}
    double scale = 1;
    if( nTheo!=0 ) { scale = nData/nTheo; }
	for (int i = 0; i < m_nBinData; i++)
    {
		m_eTheo[i] *= scale;
        
    }
	for (int i = 0; i < m_nBinData; i++)
	{
		m_eVis   [i] *= scale;
	}

}


double junoB12_simplified::GetChi2()
{
    ApplyScintillatorNL();
    Normalize();

    double chi2 = 0;
    int rebin = m_nBin / m_nBinData;
	double binWidthData = m_peBinWidth * rebin;
    int m_nData = 0;
    for(int i=0; i < m_nBinData; i++) {
        if(i*binWidthData<m_fitMinPE or binWidthData*i>m_fitMaxPE) continue;
        if( m_eDataErr[i]!=0 ) {
            chi2 += pow( (m_eData[i] - m_eTheo[i])/m_eDataErr[i], 2); 
            m_nData++;
        }
    }
    cout << "simplified B12 chi2: " << chi2 << " with nData : " << m_nData << endl;
	//if(nDoF>0) chi2 /= double(m_nData - nDoF);
	return chi2;
}


void junoB12_simplified::Plot()
{
    TH1D* hData = new TH1D("hData", "", m_nBinData, m_peMin, m_peMax);
    TH1D* hTheo = new TH1D("hTheo", "", m_nBinData, m_peMin, m_peMax);
    TH1D* hRela = new TH1D("hRela", "", m_nBinData, m_peMin, m_peMax);

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
    hRela->Write();
    out->Close();


}






