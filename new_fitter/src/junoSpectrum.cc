#include "junoSpectrum.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraphErrors.h>

using namespace std;

junoSpectrum::junoSpectrum( int nMaxBins,
                            int nMaxBinsData,
                            int nMaxBr,
                            int nMaxGam,
                            double eMin,
                            double eMax,
                            double fitMin,
                            double fitMax,
                            string nonlMode,
                            string name)
{

    std::cout << " eMin         = " << eMin << std::endl;
    std::cout << " eMax         = " << eMax << std::endl;
    std::cout << " fitMin         = " << fitMin << std::endl;
    std::cout << " fitMax         = " << fitMax << std::endl;
	std::cout << " nMaxBins     = " << nMaxBins << std::endl;
	std::cout << " nMaxBinsData = " << nMaxBinsData << std::endl;
	std::cout << " nMaxBr       = " << nMaxBr << std::endl;
	std::cout << " nMaxGam      = " << nMaxGam << std::endl;
    std::cout << " NonlMode     = " << nonlMode << std::endl;
	
    m_eMin      = eMin;
    m_eMax      = eMax;
    m_fitMin    = fitMin;
    m_fitMax    = fitMax;
	m_nBins     = nMaxBins;
    m_nBinsData = nMaxBinsData;
	m_nBranch   = nMaxBr;
	m_nGam      = nMaxGam;
    m_nonlMode  = nonlMode;
    m_name      = name;
    m_specTheoMode = junoParameters::specTheoMode;
    m_loadPrmElec  = false;

	m_binCenter = new double[nMaxBins];
    m_eVis      = new double[nMaxBins];
	
	m_eTru      = new double*[nMaxBr];
	
	m_eTruAlp   = new double [nMaxBr];

	for (int branchIdx = 0; branchIdx < nMaxBr; branchIdx++)
	{
		m_eTru   [branchIdx] = new double[nMaxBins];
	}

	m_eTruGam   = new double*[nMaxBr];
    m_eTruGamStr= new int*[nMaxBr];
	for (int branchIdx = 0; branchIdx < nMaxBr; branchIdx++)
	{
		m_eTruGam[branchIdx] = new double[nMaxGam];
        m_eTruGamStr[branchIdx] = new int[nMaxGam];
	}

    m_eData = new double[nMaxBinsData];
    m_eDataErr = new double[nMaxBinsData];

    m_nPdfBins = 1000;

	m_eTheo     = new double[nMaxBinsData];
}

junoSpectrum::~junoSpectrum()
{
	for (int branchIdx = 0; branchIdx < m_nBranch; branchIdx++)
	{
		delete [] m_eTru   [branchIdx];
	}
	for (int branchIdx = 0; branchIdx < m_nBranch; branchIdx++)
	{
		delete [] m_eTruGam[branchIdx];
	}
	delete [] m_binCenter;
    delete [] m_eVis     ;
	delete [] m_eTru     ;
	delete [] m_eTruGam  ;
	delete [] m_eTruAlp  ;
    delete [] m_eData    ;
    delete [] m_eDataErr ;
    delete [] m_eTheo    ;

}

void junoSpectrum::InitTheo()
{
    m_specTheoMode = junoParameters::specTheoMode;
    std::cout << " >>> Loading theoretical " << m_name << " spectrum <<< " << std::endl;
    std::cout << " specCalcMode = " << m_specTheoMode << std::endl;
    //std::cout << "spectrum calc mode " << m_specTheoMode << endl;

    string theofile;
    if (m_specTheoMode == "calc")
        theofile = "./data/spectrum/theo/" + m_name + "_theo.root";
    if (m_specTheoMode == "sim") {
        theofile = "./data/spectrum/theo/" + m_name + "_sim.root"; 
    }

    TheoHistTree(theofile);
}

void junoSpectrum::TheoHistTree(string theofile)
{
        TFile* ff = new TFile(theofile.c_str());
        if(!ff) cout << " >>> No theoretical spectrum file " << theofile << endl;
        TTree* tt = (TTree*)ff->Get("T");
        double eGamma[10];
        int eGammaStr[10];
        int nGamma, branchNumber;
        double branchRatio;
        double weight;
        tt->SetBranchAddress("num"      ,&branchNumber);
        tt->SetBranchAddress("BR"       ,&branchRatio);
        tt->SetBranchAddress("numPhoton",&nGamma);
        tt->SetBranchAddress("photonE"  ,eGamma);
        tt->SetBranchAddress("photonName", eGammaStr);
        m_nBranch = tt->GetEntries();
        cout << " >>> Total Branch Number = " << m_nBranch << endl;
        for (int branchIdx=0;branchIdx!=tt->GetEntries();branchIdx++){ 
            tt->GetEntry(branchIdx);
            cout << " >>> " << branchIdx << " with "<< nGamma << " gamma and branch ratio is " << branchRatio << endl;
            // gamma from this branch :
            //m_nGamma[branchIdx] = nGamma;
            for(int gamIdx=0; gamIdx<nGamma; gamIdx++) {
                cout << "gamma energy = " << eGammaStr[gamIdx] << " keV" << endl;
                m_eTruGam[branchIdx][gamIdx] = eGamma[gamIdx];
                m_eTruGamStr[branchIdx][gamIdx] = eGammaStr[gamIdx];
            }

            // beta
            TH1F* electronHist = (TH1F*)ff->Get(Form("hh%d",branchNumber));
            for (int binIdx=0;binIdx!=m_nBins;binIdx++)
            {
                weight = electronHist->Interpolate(m_binCenter[binIdx]);
                m_eTru[branchIdx][binIdx] = branchRatio * weight;
            }
            delete electronHist;
        }

        delete tt;
        delete ff;

}

void junoSpectrum::InitData()
{
    std::cout << " >>> Loading data " << m_name << "spectrum <<< " << std::endl;
    string fileName = "./data/spectrum/data/" + m_name + "_data.root";
    DataHistTree(fileName);
}

void junoSpectrum::DataHistTree(string fileName)
{
    // Old spectrum data loading ... -> histogram form
    //TFile* ff = new TFile(fileName.c_str());
    //if(!ff) cout << "No such data file " << fileName << endl;
    //TH1D* sigH = (TH1D*)ff->Get(m_name.c_str());
	//for (int i=0;i!=m_nBinsData;i++)
	//{
	//	double content = sigH->GetBinContent(i+1);
	//	double error   = sigH->GetBinError  (i+1);
	//	m_eData   [i] = content;
	//	m_eDataErr[i] = error;
	//}
	//delete sigH;

    // Tree form spectrum data loading...
    double m_scale = 3300.371/2.223;

    TH1D* sigH = new TH1D("B12_data", "", m_nBinsData, m_eMin, m_eMax);

    TFile* ff = new TFile(fileName.c_str());
    if(!ff) cout << "No such data file " << fileName << endl;
    TTree* tt = (TTree*)ff->Get("B12");
    double m_totpe;
    tt->SetBranchAddress("totpe", &m_totpe);
    for(int i=0; i<tt->GetEntries(); i++) {
        tt->GetEntry(i);
        double tmp_Evis = m_totpe / m_scale;
        sigH->Fill(tmp_Evis);
    }
    
    for (int i=0; i<m_nBinsData; i++) {
		double content = sigH->GetBinContent(i+1);
		double error   = sigH->GetBinError  (i+1);
		m_eData   [i] = content;
		m_eDataErr[i] = error ;
    }
    
	delete sigH;
    delete tt;
    delete ff;
}


void junoSpectrum::LoadData()
{
    // Initialization Part
    m_binWidth = (m_eMax - m_eMin) / double(m_nBins);
    m_fitMinBin = int((m_fitMin - m_eMin)/m_binWidth);
    m_fitMinBin = int((m_fitMax - m_eMin)/m_binWidth);

    // theoretical spectrum
    for (int branchIdx=0; branchIdx<m_nBranch; branchIdx++) {
        // beta continuous spectrum
        for (int i=0; i<m_nBins; i++) {
            m_binCenter[i] = m_eMin + m_binWidth*(i+0.5);
            m_eTru[branchIdx][0] = 0;
        }
        // gamma energy
        for (int gamIdx=0; gamIdx<m_nGam; gamIdx++) {
            m_eTruGam[branchIdx][gamIdx] = 0;
        }
        // alpha energy
        m_eTruAlp[branchIdx] = 0;
    }

    // data spectrum
    for (int i=0; i<m_nBinsData; i++) {
        m_eData[i] = 0;
        m_eDataErr[i] = 0;
    }

    // Data Loading ...
    InitTheo();
    InitData();

}


void junoSpectrum::ApplyScintillatorNL()
{
    if (!m_loadPrmElec)
        LoadPrmElecDist();

    for (int i=0; i<m_nBins; i++) {
        m_eVis[i] = 0;
    }

    int newBin;
    int newBinLow, newBinHig;
    double bias;
	double eVisGam[m_nBranch];

    cout << "Current fitting params " << electronQuench::getkA() << " " << electronCerenkov::getkC() <<" " << electronQuench::getEnergyScale()<<endl;

    for(int branchIdx=0; branchIdx<m_nBranch; branchIdx++) {
        // Nonlinearity on gamma
        for (int gamIdx=0; gamIdx<m_nGam; gamIdx++) {
			if(m_eTruGam[branchIdx][gamIdx]==0) break;  // No more gamma in such branch
            eVisGam[branchIdx] += EvisGamma(m_eTruGamStr[branchIdx][gamIdx]) * m_eTruGam[branchIdx][gamIdx] ;
            //cout << "gamma " << m_eTruGam[branchIdx][gamIdx] << " " << EvisGamma(m_eTruGamStr[branchIdx][gamIdx]) << " "<< eVisGam[branchIdx] << endl;
        }
    }

        // Nonlinearity on beta
        for (int i=0; i<m_nBins; i++) {
            double eTru = m_binCenter[i];
            double eVisElec = eTru;
            if (m_nonlMode == "histogram")
                eVisElec *= (electronQuench::ScintillatorNL(eTru) + electronCerenkov::getCerenkovPE(eTru));
            if (m_nonlMode == "analytic")
                eVisElec *= (electronResponse::calcElecNonl(eTru));

            for(int branchIdx=0; branchIdx<m_nBranch; branchIdx++){
                double eVis = eVisElec + eVisGam[branchIdx];
                // if(m_name=="B12" and branchIdx==2) eVis += 0.5; // 3-alpha branch -> No 3-alpha branch in gendecay sim??
                newBinLow     = int((eVis-m_eMin)/m_binWidth);
                newBinHig     = int((eVis-m_eMin)/m_binWidth) + 1;
			    bias          = (eVis - m_eMin - newBinLow * m_binWidth)/m_binWidth;

                if (newBinLow<m_nBins) m_eVis[newBinLow] += (1-bias) * m_eTru[branchIdx][i];
			    if (newBinHig<m_nBins) m_eVis[newBinHig] += bias * m_eTru   [branchIdx][i];
            }
        }
}


void junoSpectrum::LoadPrmElecDist()
{
    cout << " >>> Start Loading Primary electron distribution ..." << endl;

    TFile* file = new TFile(junoParameters::gammaPdf_File.c_str(), "read");
    for (int branchIdx=0; branchIdx<m_nBranch; branchIdx++) {
        
        for (int gamIdx=0; gamIdx<m_nGam; gamIdx++) { 
			if(m_eTruGam[branchIdx][gamIdx]==0) break;  // No more gamma in such branch
            string eTru = to_string(m_eTruGamStr[branchIdx][gamIdx]);

            string pdfName = "gamma"+eTru+"keV";
            std::cout << "Loading PrmElecDist for " << pdfName << std::endl;
            TH1D* gGammaPdf = (TH1D*)file->Get(pdfName.c_str());
            if(!gGammaPdf) cout << "No Such Pdf : " << pdfName << endl;
            
            int tmp_PdfMaxEtrue;
            double* tmp_pdfEtrue = new double[m_nPdfBins];
            double* tmp_pdfProb = new double[m_nPdfBins];

            for(int i=0; i<gGammaPdf->GetNbinsX(); i++)  {
                tmp_pdfEtrue[i] = gGammaPdf->GetBinCenter(i+1);
                tmp_pdfProb[i] = gGammaPdf->GetBinContent(i+1);
                if (tmp_pdfProb[i] == 0) {  tmp_PdfMaxEtrue = i; break;}
            }

            mapPdfMaxEtrue.insert(pair<int, int>(m_eTruGamStr[branchIdx][gamIdx], tmp_PdfMaxEtrue));
            mapPdfEtrue.insert(pair<int, double*>(m_eTruGamStr[branchIdx][gamIdx], tmp_pdfEtrue));
            mapPdfProb.insert(pair<int, double*> (m_eTruGamStr[branchIdx][gamIdx], tmp_pdfProb));

            pdfName = "antigamma"+eTru+"keV";
            std::cout << "Loading PrmPosiDist for " << pdfName << std::endl;
            TH1D* gAntiGammaPdf = (TH1D*)file->Get(pdfName.c_str());
            if(!gGammaPdf) cout << "No Such Pdf : " << pdfName << endl;
            double* tmp_pdfEtruePosi = new double[m_nPdfBins];
            double* tmp_pdfProbPosi = new double[m_nPdfBins];
            for(int i=0; i<gGammaPdf->GetNbinsX(); i++)  {
                tmp_pdfEtruePosi[i] = gAntiGammaPdf->GetBinCenter(i+1);
                tmp_pdfProbPosi[i] = gAntiGammaPdf->GetBinContent(i+1);
            }

            mapPdfEtruePosi.insert(pair<int, double*>(m_eTruGamStr[branchIdx][gamIdx], tmp_pdfEtruePosi));
            mapPdfProbPosi.insert(pair<int, double*> (m_eTruGamStr[branchIdx][gamIdx], tmp_pdfProbPosi));


            delete gGammaPdf;
            delete gAntiGammaPdf;
        }

    }
    
    delete file;

    m_loadPrmElec = true;
}

double junoSpectrum::EvisGamma(int Etrue)
{
    int gamPdfMaxEtrue       = mapPdfMaxEtrue[Etrue];
    double* gamPdfEtrue      = mapPdfEtrue[Etrue];
    double* gamPdfProb       = mapPdfProb[Etrue];
    double* gamPdfEtruePosi  = mapPdfEtruePosi[Etrue];
    double* gamPdfProbPosi   = mapPdfProbPosi[Etrue];

    //cout << Etrue << " " << gamPdfMaxEtrue << endl;
    double numerator = 0.; double denominator = 1e6;
    for(int iBin=0;  iBin<m_nPdfBins; iBin++) {
        double E1 = gamPdfEtrue[iBin];

        double prob1 = gamPdfProb[iBin];

        double NPE = electronQuench::ScintillatorPE(E1) + electronCerenkov::getCerPE(E1);
        //cout << "beta " << E1 << " " << prob1 << " " << NPE << endl;

        numerator   += NPE * prob1;
    } 
    for(int iBin=0;  iBin<m_nPdfBins; iBin++) {
        double E1 = gamPdfEtruePosi[iBin];

        double prob1 = gamPdfProbPosi[iBin];

        double NPE = electronQuench::ScintillatorPE(E1) + electronCerenkov::getCerPE(E1) + 695.53*2; 
        //cout << "beta+ " << E1 << " " << prob1 << " " << NPE << endl;

        numerator   += NPE * prob1;
    } 

    if(denominator ==0) { cout << " >> Error Happens While CalculateGammaNL <<<" << endl; return 0;}
    //return numerator/denominator;   // return totpe prediction value
    //cout << numerator << " " << denominator << " " << electronQuench::getEnergyScale() << " " << Etrue/1000. << " " << numerator/denominator/electronQuench::getEnergyScale()/(Etrue/1000.)<< endl;
    return numerator / denominator / electronQuench::getEnergyScale() / (Etrue/1000.);
}


void junoSpectrum::Normalize()
{
    // Normalize spectrum for data and pred
	int   rebin = m_nBins/m_nBinsData;
	double binWidthData = m_binWidth * rebin;
	double nTheo = 0;
	double nData = 0;
	for (int i = 0; i < m_nBinsData; i++)
	{
		m_eTheo[i] = 0;
        for (int j = 0; j < rebin; j++){
			m_eTheo[i] += m_eVis[i*rebin+j];
        } 
            
		if(i*binWidthData>m_fitMin && i*binWidthData<m_fitMax)
		{
			nTheo += m_eTheo[i];
			nData += m_eData[i];
		}
	}
    double scale = 1;
    if( nTheo!=0 ) { scale = nData/nTheo; }
	for (int i = 0; i < m_nBinsData; i++)
    {
		m_eTheo[i] *= scale;
        
    }
	for (int i = 0; i < m_nBins; i++)
	{
		m_eVis   [i] *= scale;
	}

}


double junoSpectrum::GetChi2()
{
    ApplyScintillatorNL();
    Normalize();
    
    double chi2 = 0;
    int rebin = m_nBins / m_nBinsData;
    double binWidthData = m_binWidth * rebin;
    m_nData = 0;
    for(int i=0; i < m_nBinsData; i++) {
        if(i*binWidthData<m_fitMin or binWidthData*i>m_fitMax-0.1) continue;
        if( m_eDataErr[i]!=0 ) {
            //cout << m_eData[i] << " " << m_eTheo[i] << " " << m_eDataErr[i] << endl;
            chi2 += pow( (m_eData[i] - m_eTheo[i])/m_eDataErr[i], 2); 
            m_nData++;
        }
    }
    cout << m_name << " chi2: " << chi2 << " with nData : " << m_nData << endl;
	//if(nDoF>0) chi2 /= double(m_nData - nDoF);
	return chi2;

}


void junoSpectrum::Plot()
{
    cout << " >>> Plot junoSpectum <<< " << endl;

    TH1D* hData = new TH1D("hData", "", m_nBinsData, m_eMin, m_eMax);
    TH1D* hTheo = new TH1D("hTheo", "", m_nBinsData, m_eMin, m_eMax);
    TH1D* hRela = new TH1D("hRela", "", m_nBinsData, m_eMin, m_eMax);

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
    hRela->SetMarkerColor(kPink+2);


    for(int i=0; i<m_nBinsData; i++) {
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

    //TCanvas * tmpC = new TCanvas();
    //Float_t small = 1e-5;
    //tmpC->Divide(1,2, small, small);
    //gStyle->SetPadBorderMode(0);
    //gStyle->SetFrameBorderMode(0);
    //tmpC->cd(1);
    //gPad->SetBottomMargin(small);
    //hData->Draw("PEX0");
    //hTheo->Draw("PEXO SAME");
    //TLegend* ll = new TLegend();
    //ll->AddEntry(hData, "B12 Data(Sim)", "l");
    //ll->AddEntry(hTheo, "B12 Theo", "l");
    //ll->Draw("SAME");
    //tmpC->cd(2);
    //gPad->SetTopMargin(small);
    //gPad->SetTickx();
    //hRela->Draw("PEX0");

    TFile* out = new TFile("spectrum.root", "recreate");
    hData->Write();
    hTheo->Write();
    hRela->Write();
    out->Close();
}

