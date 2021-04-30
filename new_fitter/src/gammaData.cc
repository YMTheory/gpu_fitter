#include "gammaData.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include <TFile.h>
#include <TH1.h>
# include <TSystem.h>
#include <TMath.h>
#include <TRandom.h>

using namespace std;

//std::string gammaData::m_calcOption = "prmelec";
std::string gammaData::m_calcOption = junoParameters::m_calcOption;

std::string gammaData::m_nonlMode = "histogram";
//std::string gammaData::m_nonlMode = junoParameters::m_nonlMode;

gammaData::gammaData( std::string name,
                        double minPE,
                        double maxPE,
                        int nbins
                      ) 
{
    m_name = name;
    m_minPE = minPE;
    m_maxPE = maxPE;
    m_nbins = nbins;

    m_calcOption = junoParameters::m_calcOption;
    m_nonlMode   = junoParameters::m_nonlMode;

    m_loadData = false;

    //hEmean = new TH1D(m_name.c_str(), "", 10000, 0, 10);
    hPEmean = new TH1D((m_name+"pe").c_str(), "", 5000, 0, 20000);
    //hCerPEmean = new TH1D((m_name+"cerpe").c_str(), "", 500, 0, 2000);
    //hSctPEmean = new TH1D((m_name+"sctpe").c_str(), "", 5000, 0, 20000);

}

gammaData::~gammaData() 
{}

void gammaData::LoadGammaData()
{
    cout << " >>> Loading Naked Gamma " << m_name << " Data <<< " << endl;

    LoadGammaPEComp();

    ifstream in; in.open(junoParameters::gammaLSNL_File);
    string line;

    double scale = 3300.371/2.223;
    string tmp_name; double tmp_E, tmp_totPE, tmp_totPESigma, tmp_EvisError, tmp_totPEerr, tmp_totPESigmaerr;
    while(getline(in,line)){
        istringstream ss(line);
        //ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPESigma >> tmp_EvisError ;
        ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPEerr >> tmp_totPESigma >> tmp_totPESigmaerr;
        if(tmp_name == m_name) {
            m_Etrue = tmp_E;
            m_totpeData = tmp_totPE;
            m_nonlData = tmp_totPE/scale/tmp_E;
            //m_nonlDataErr = tmp_EvisError*tmp_totPE/scale/tmp_E;
            //m_nonlDataErr = tmp_totPEerr/scale/tmp_E;
            m_nonlDataErr = 0.002;
            m_Evis = tmp_totPE/scale;
            m_resData = tmp_totPESigma/tmp_totPE;
            //m_resDataErr = 0.01 * tmp_totPESigma/tmp_totPE;
            m_resDataErr = TMath::Sqrt(tmp_totPESigmaerr*tmp_totPESigmaerr/tmp_totPE/tmp_totPE + tmp_totPEerr*tmp_totPEerr*tmp_totPESigma*tmp_totPESigma/tmp_totPE/tmp_totPE/tmp_totPE/tmp_totPE);
            break;
        }
    }
    in.close();
}


void gammaData::LoadGammaPEComp()
{
    cout << " >>> Load Gamma PE Component File <<< " << endl;
    string line;
    string tmp_name; double tmp_E, tmp_pe, tmp_cerpe, tmp_sctpe;
    ifstream in; in.open(junoParameters::gammaPE_File);
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmp_name >> tmp_E >> tmp_pe >> tmp_cerpe >> tmp_sctpe;
        if(tmp_name == m_name) {
            m_cerPEData = tmp_cerpe;
            m_sctPEData = tmp_sctpe;
        }
   } 
    in.close();
}


void gammaData::LoadPrimElecDist()
{
    cout << " >>> Load Primary Electron Distribution <<< " << endl;

    TFile* file = new TFile(junoParameters::gammaPdf_File.c_str(), "read");
    if (!file) cout << " No such input primary electron file " << endl;
    //string pdfName = m_name;
    string pdfName = "gamma" + m_name;
    TH1D* gGammaPdf = (TH1D*)file->Get(pdfName.c_str());
    //double intvalue = gGammaPdf->Integral() ;
    //gGammaPdf->Scale(1/intvalue);
    if (!gGammaPdf) cout << " No such Pdf : " << pdfName << endl;

    for(int i=0; i<gGammaPdf->GetNbinsX(); i++) {
        m_pdf_eTrue[i] = gGammaPdf->GetBinCenter(i+1);
        m_pdf_prob[i]  = gGammaPdf->GetBinContent(i+1);
        if (m_pdf_prob[i] == 0 ) { 
            m_max_eTrue = i;
            break;
        }
    }
    delete gGammaPdf;

    // Load Positron distribution
    pdfName = "antigamma" + m_name;
    TH1D* gAntiGammaPdf = (TH1D*)file->Get(pdfName.c_str());
    if (!gAntiGammaPdf) cout << "No such positron Pdf : " << pdfName << endl;
    for(int i=0; i<gAntiGammaPdf->GetNbinsX(); i++) {
        m_pdf_eTrue_posi[i] = gAntiGammaPdf->GetBinCenter(i+1);
        m_pdf_prob_posi[i] = gAntiGammaPdf->GetBinContent(i+1);
    }
    delete gAntiGammaPdf;

    file->Close();

}

void gammaData::LoadPrimElecSamples()
{
    cout << " >>> Load Primary Electron in Single Event <<< " << endl;
    string filename = "./data/gamma/" + m_name + "_all.root";
    TFile* file = new TFile(filename.c_str(), "read");
    if (!file) cout << " No such input file: " << filename << endl;
    elec_hist = (TH2D*)file->Get(m_name.c_str());

}




void gammaData::LoadData()
{
    LoadGammaData();

    if (m_calcOption == "prmelec")
        LoadPrimElecDist();

    if (m_calcOption == "twolayer") {
        LoadPrimElecSamples();
    }

    m_loadData = true;
}


void gammaData::calcGammaResponse()
{
    if (!m_loadData) LoadData();
    m_scale = electronQuench::getEnergyScale();

    if (m_calcOption == "prmelec") {
        double numerator = 0;
        double denominator = 0;
        
        double numerator1 = 0;
        double numerator2 = 0;
        double numerator3 = 0;
        double denominator1 = 1e6;

        // electron part
        for(int iBin=1; iBin<m_max_eTrue; iBin++) {
            double E1 = m_pdf_eTrue[iBin-1];
            double E2 = m_pdf_eTrue[iBin];
            double prob1 = m_pdf_prob[iBin-1];
            double prob2 = m_pdf_prob[iBin];

            // based on DYB method...
            double fNL1, fNL2;

            //cout << m_name << " " << iBin << " " << E1 <<" " << prob1 << " " << E2 << " " << prob2 << endl;
            if(m_nonlMode == "histogram") {
                fNL1 = electronQuench::ScintillatorNL(E1) + electronCerenkov::getCerenkovPE(E1);
                fNL2 = electronQuench::ScintillatorNL(E2) + electronCerenkov::getCerenkovPE(E2);
            }

            if (m_nonlMode == "analytic") {
                fNL1 = electronResponse::calcElecNonl(E1);
                fNL2 = electronResponse::calcElecNonl(E2);
            }

            numerator += (prob1*E1*fNL1 + prob2*E2*fNL2) * (E2-E1) / 2.;
            denominator += (prob1*E1 + prob2*E2) * (E2-E1) / 2.;

            double NPE1, NPE2, NPE3;
            if (m_nonlMode == "histogram" ) {
                NPE1 = electronQuench::ScintillatorPE(E1) + electronCerenkov::getCerPE(E1);
                NPE2 = electronQuench::ScintillatorPE(E1) ;
                NPE3 = electronCerenkov::getCerPE(E1);
            }

            numerator1 += NPE1 * prob1;
            numerator2 += NPE2 * prob1;
            numerator3 += NPE3 * prob1;
            //denominator1 += (prob1 + prob2) /2;
        }

        for (int ibin=1; ibin<m_nMaxPdf; ibin++) {
            double E1 = m_pdf_eTrue_posi[ibin-1];
            double E2 = m_pdf_eTrue_posi[ibin];
            double prob1 = m_pdf_prob_posi[ibin-1];
            double prob2 = m_pdf_prob_posi[ibin];

            double fNL1, fNL2;
            fNL1 = electronQuench::ScintillatorNL(E1) + electronCerenkov::getCerenkovPE(E1);
            fNL2 = electronQuench::ScintillatorNL(E2) + electronCerenkov::getCerenkovPE(E2);
        
            numerator += (prob1*E1*fNL1 + prob2*E2*fNL2) * (E2-E1) / 2.;
            denominator += (prob1*E1 + prob2*E2) * (E2-E1) / 2.;

            double NPE1, NPE2, NPE3;
            if (m_nonlMode == "histogram" ) {
                NPE1 = electronQuench::ScintillatorPE(E1) + electronCerenkov::getCerPE(E1) + 695.53*2; 
                NPE2 = electronQuench::ScintillatorPE(E1) ;
                NPE3 = electronCerenkov::getCerPE(E1);
            }

            numerator1 += NPE1 * prob1;
            numerator2 += NPE2 * prob1;
            numerator3 += NPE3 * prob1;
        }

        if (denominator == 0) cout << "Errors Happend While Using GammaPdf Calculation..." << endl;

        m_nonlCalc = numerator / denominator;
        //m_totpeCalc = m_nonlCalc * m_Etrue * m_scale;
        m_totpeCalc = numerator1 / denominator1;
        m_sctPE = numerator2 /denominator1;
        m_cerPE = numerator3 / denominator1;
        m_nonlCalc1 = numerator1 /denominator1 / m_scale / m_Etrue;

        //cout << m_name << " " << m_totpeData << " " << m_totpeCalc << " " << m_nonlData << " " << m_nonlCalc << " " << m_nonlCalc1 << endl;


    }  else if (m_calcOption == "twolayer") {
        if (not electronResponse::getLoadResol()) electronResponse::loadElecResol();

        //hEmean->Reset();
        hPEmean->Reset();
        //hCerPEmean->Reset();
        //hSctPEmean->Reset();

        for (int iSample=0; iSample<m_nSamples; iSample++) {
            // apply Nonlinearity curve
            double tmp_mean = 0;
            double tmp_pe= 0;
            double tmp_cerpe = 0;
            double tmp_sctpe = 0;
            double tmp_sigma = 0;
            for (int iSec=0; iSec<100; iSec++) {
                double tmp_E = elec_hist->GetBinContent(iSample+1, iSec+1);
                if (tmp_E == 0) break;
                double tmp_sigpe;
                double tmp_sigcerpe;
                double tmp_sigsctpe;
                double resol;
                if(m_nonlMode == "histogram") {
                    tmp_pe += electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E) ;
                    //tmp_sigsctpe = tmp_E * (electronQuench::ScintillatorNL(tmp_E)) * m_scale;
                    //tmp_sigcerpe = tmp_E * (electronCerenkov::getCerenkovPE(tmp_E)) * m_scale;

                    tmp_sigma += electronResponse::gElecResol->Eval(tmp_E) * electronResponse::gElecResol->Eval(tmp_E);

                    //tmp_sigpe = (electronQuench::ScintillatorPE(tmp_E) +electronCerenkov::getCerPE(tmp_E)); // alreadly include nonl
                    //tmp_sigcerpe = electronCerenkov::getCerPE(tmp_E);
                    //tmp_sigsctpe = tmp_sigpe - tmp_sigcerpe;
                    //resol = electronResponse::gElecResol->Eval(tmp_E);
                    //tmp_Edep += TMath::Gaus(0, 1) * resol*tmp_Edep;
                    //tmp_sigpe += TMath::Gaus(0, 1) * resol * tmp_sigpe;
                    //tmp_sigcerpe += TMath::Gaus(0, 1) * resol * tmp_sigcerpe;
                    //tmp_sigsctpe += TMath::Gaus(0, 1) * resol * tmp_sigsctpe;
                }
                //if(m_nonlMode == "analytic") {
                //    tmp_Edep = tmp_E * electronResponse::calcElecNonl(tmp_E);
                //}

                //tmp_mean += tmp_Edep;
                //tmp_cerpe += tmp_sigcerpe;
                //tmp_sctpe += tmp_sigsctpe;
                //tmp_pe += tmp_sigpe;

                tmp_sigma = TMath::Sqrt(tmp_sigma);
            }
            double sample_pe = gRandom->Gaus(tmp_pe, tmp_sigma);
            //double sample_pe = tmp_pe;

            //hEmean->Fill(tmp_mean);
            hPEmean->Fill(sample_pe);
            //hCerPEmean->Fill(tmp_cerpe);
            //hSctPEmean->Fill(tmp_sctpe);
        }

        // calculate pe distribution
        //double mean_Edep = 0;
        //double mean_pe = 0;
        //for (int i=0; i<m_nSamples; i++) {
        //    mean_Edep += m_mean[i];
        //    //mean_pe += m_meanpe[i];
        //}
        //mean_Edep /= m_nSamples;
        //mean_pe /= m_nSamples;
        //m_totpeCalc = mean_pe;
        m_totpeCalc = hPEmean->GetMean();
        m_nonlCalc = hPEmean->GetMean() / m_scale / m_Etrue;
        m_nonlCalc1 = hPEmean->GetMean() / m_scale / m_Etrue;
        //m_sctPE = hSctPEmean->GetMean();
        //m_cerPE = hCerPEmean->GetMean();
        //m_nonlCalc = hEmean->GetMean() / m_Etrue;
        //m_nonlCalc1 = m_totpeCalc / m_scale / m_Etrue ;
        //cout << m_name << " " << m_totpeData << " " << m_totpeCalc << " \n" 
        //     << m_nonlData << " " << m_nonlCalc << " " << m_nonlCalc1 << endl;
    }

    // pull term for gamma :
    m_nuGamma = junoParameters::m_nuGamma;
    m_nonlCalc *= (1+m_nuGamma);
 
}


double gammaData::GetChi2()
{
    double chi2 = 0;

    // calculate totpe sigma
    calcGammaResponse();

    chi2 += (m_nonlCalc1 - m_nonlData) * (m_nonlCalc1 - m_nonlData) / m_nonlDataErr / m_nonlDataErr;

    return chi2;
}



void gammaData::SaveHist()
{
    TFile* out = new TFile((m_name+"hist.root").c_str(), "recreate");    
    hEmean->Write();
    hPEmean->Write();
    hCerPEmean->Write();
    hSctPEmean->Write();
    out->Close();
}
