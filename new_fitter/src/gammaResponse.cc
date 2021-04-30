#include "gammaResponse.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TFile.h>
#include <TRandom.h>
#include <TF1.h>
#include <TTree.h>

gammaResponse::gammaResponse(string name, int nBins, double peMin, double peMax) {
    m_name = name;
    m_nBins = nBins;
    m_peMin = peMin;
    m_peMax = peMax;

    cout << " >>> " << m_name << " " << m_nBins << " bins between [" << m_peMin << ", " << m_peMax <<"]" <<endl;
    hCalc = new TH1D((m_name+"_calc").c_str(), "", m_nBins, m_peMin, m_peMax); 
    hData = new TH1D((m_name+"_data").c_str(), "", m_nBins, m_peMin, m_peMax); 
    func  = new TF1("func", "[0] * TMath::Exp(-(x-[1])*(x-[1])/2/[2]/[2])", m_peMin, m_peMax);

    // init
    m_amp = 180;
    m_loadData  = false;
    m_loadPrm   = false;
    m_doSpecFit = true;
}

gammaResponse::~gammaResponse()
{
    delete hPrmElec;
    delete hPrmPosi;

    delete hCalc;
}


void gammaResponse::LoadData()
{
    cout << " >>> Loading Bared Gamma " << m_name << " Data <<< " << endl;

    ifstream in; in.open(junoParameters::gammaLSNL_File);
    if(!in) cout << "Error: No file "<< junoParameters::gammaLSNL_File << std::endl;
    string line;

    double scale = 3134.078/2.223;
    string tmp_name; double tmp_E, tmp_totPE, tmp_totPESigma, tmp_EvisError, tmp_totPEerr, tmp_totPESigmaerr;
    while(getline(in,line)){
        istringstream ss(line);
        //ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPESigma >> tmp_EvisError ;
        ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPEerr >> tmp_totPESigma >> tmp_totPESigmaerr;
        if(tmp_name == m_name) {
            m_Etrue = tmp_E;
            m_totpeData = tmp_totPE;
            m_nonlData = tmp_totPE/scale/tmp_E;
            //m_nonlDataErr = tmp_totPEerr/scale/tmp_E;
            m_nonlDataErr = 0.001;
            m_Evis = tmp_totPE/scale;
            m_resData = tmp_totPESigma/tmp_totPE;
            m_resDataErr = TMath::Sqrt(tmp_totPESigmaerr*tmp_totPESigmaerr/tmp_totPE/tmp_totPE + tmp_totPEerr*tmp_totPEerr*tmp_totPESigma*tmp_totPESigma/tmp_totPE/tmp_totPE/tmp_totPE/tmp_totPE);
            break;
        }
    }
    in.close();

    if (m_doSpecFit and (m_name != "gamma4440" and m_name!="gamma3215" ) ) {
        TFile* infile = new TFile(("./data/gamma/spectrum/" + m_name + "_totpe.root").c_str(), "read");
        if(!infile) cout << "No such gamma spectrum file!" << endl;
        int m_totpe;
        TTree* tt = (TTree*)infile->Get("evt");
        tt->SetBranchAddress("totalPE", &m_totpe);
        for(int i=0; i<tt->GetEntries(); i++) {
            tt->GetEntry(i);
            hData->Fill(m_totpe);
        }
    }


    LoadPrmBeta();
}

void gammaResponse::LoadPrmBeta()
{
    cout << " >>> Load Primary Electron in Single Event for " << m_name << " <<< " << endl;
    string filename = "./data/gamma/" + m_name + "_J19.root";
    TFile* file = new TFile(filename.c_str(), "read");
    if (!file) cout << " No such input file: " << filename << endl;
    hPrmElec = (TH2D*)file->Get((m_name+"_elec").c_str());
    hPrmPosi = (TH2D*)file->Get((m_name+"_posi").c_str());
    
    m_loadPrm = true;
}


void gammaResponse::preCalculation()
{
    //if (not electronResponse::getLoadResol()) electronResponse::loadElecResol();
    electronResponse::SetParameters();

    for (int index=0; index<m_nSamples; index++) {
        double tmp_pe = 0;
        double tmp_sigma = 0;
        // electron
        for (int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmElec->GetBinContent(index+1, iSec+1);    
            if (tmp_E == 0) break;
            tmp_pe += electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E);
            tmp_sigma += TMath::Power(electronResponse::fElecResol->Eval(tmp_E), 2);
            //tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
        }

        // positron
        for(int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmPosi->GetBinContent(index+1, iSec+1);
            if (tmp_E == 0) break;
            tmp_pe += electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E) + 2*660.8;
            tmp_sigma += TMath::Power(electronResponse::fElecResol->Eval(tmp_E), 2);
            //tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
            tmp_sigma += 2*TMath::Power(27.07, 2);

        }

        tmp_sigma = TMath::Sqrt(tmp_sigma);

        m_pemean[index] = tmp_pe;
        m_pesigma[index] = tmp_sigma;

    }

}


void gammaResponse::Prediction()
{
    preCalculation();
    double tot_mean = 0;
    for (int i=0; i<m_nSamples; i++)  {
        tot_mean += m_pemean[i];
    }
    tot_mean /= m_nSamples;

    double tot_sigma = 0;
    for (int i=0; i<m_nSamples; i++) {
        tot_sigma += (m_pemean[i]-tot_mean)*(m_pemean[i]-tot_mean) + m_pesigma[i]*m_pesigma[i];
    }
    tot_sigma = TMath::Sqrt(tot_sigma/m_nSamples);

    m_totpeCalc = tot_mean;
    m_totpeSigmaCalc = tot_sigma;

    m_nonlCalc = tot_mean / electronQuench::getEnergyScale() / m_Etrue;
    m_resCalc  = tot_sigma / tot_mean;

}



double gammaResponse::SampleGamEnergy(int index)
{
    if (index >= m_nEvents or index < 0) { cout << "Incorrect Index !" << endl; return 0; } 
    else {
        double tmp_pe = 0;
        double tmp_sigma = 0;
        // electron
        for (int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmElec->GetBinContent(index+1, iSec+1);    
            if (tmp_E == 0) break;
            tmp_pe += electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E);
            tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
        }

        // positron
        for(int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmPosi->GetBinContent(index+1, iSec+1);
            if (tmp_E == 0) break;
            tmp_pe += electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E) + 2*660.8;
            tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
            tmp_sigma += 2*TMath::Power(27.07, 2);
            
        }

        tmp_sigma = TMath::Sqrt(tmp_sigma);
    
        double sample_pe = gRandom->Gaus(tmp_pe, tmp_sigma);
        //if (sample_pe < 8000) { // for small p.e. events check...
        //    cout << index << " " << tmp_pe << " " << tmp_sigma <<   endl;
        //}
        return sample_pe;
    }
}




void gammaResponse::calcGamResponse()
{
    if (not electronResponse::getLoadResol()) electronResponse::loadElecResol();
    hCalc->Reset();

    for (int iSample = 0; iSample<m_nSamples; iSample++) {
        int index = gRandom->Integer(5000);
        double sample_pe = SampleGamEnergy(index);
        hCalc->Fill(sample_pe); 
    }
    //hCalc->Fit("gaus", "Q0");
    //double pe_mean  = hCalc->GetFunction("gaus")->GetParameter(1);
    //double pe_sigma = hCalc->GetFunction("gaus")->GetParameter(2);

    //double pe_mean  = hCalc->GetMean();
    //double pe_sigma = hCalc->GetStdDev();

    hCalc->Fit("gaus", "Q0");
    double pe_mean =  hCalc->GetFunction("gaus")->GetParameter(1);
    double pe_sigma = hCalc->GetFunction("gaus")->GetParameter(2);

    m_nonlCalc = pe_mean / electronQuench::getEnergyScale() / m_Etrue;
    m_resCalc = pe_sigma / pe_mean;

}


double gammaResponse::GetChi2()
{
    double chi2 = 0;

    //calcGamResponse();
    Prediction();
    
    func->SetParameters(m_amp, m_totpeCalc, m_totpeSigmaCalc);

    if (not m_doSpecFit) 
        chi2 += TMath::Power((m_nonlData - m_nonlCalc)/m_nonlDataErr, 2);

    else{
        for (int i=0; i<m_nBins; i++) {
            double m_err = hData->GetBinError(i);
            if (m_err != 0) {
                double m_bins = hData->GetBinCenter(i);
                double m_data = hData->GetBinContent(i);
                //double m_calc = hCalc->GetBinContent(i);
                double m_calc = func->Eval(m_bins);
                //cout << "spectrum fitting " << m_bins << " " << m_data << " " << m_calc << " " << m_err<< endl;

                chi2 += TMath::Power((m_calc - m_data)/m_err, 2);
                //cout << hData->GetBinCenter(i) << " " << hCalc->GetBinCenter(i) << " " << m_calc << " " << m_data << " " << m_err << " " << chi2 << endl;
            }
        }
    }
        
    cout << m_name << " " << chi2 << endl;
    return chi2 ;
}


void gammaResponse::SaveHist()
{
    string fileName = m_name + "hist.root";
    TFile* outfile = new TFile(fileName.c_str(), "recreate");
    hData->Write();
    //hCalc->Write();
    func->Write();
    outfile->Close();
}

