#include "electronQuench.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1.h>
#include "junoParameters.hh"

using namespace std;

double electronQuench::m_kA = 1;
double electronQuench::m_birk1     = 6.5e-3; //5.72857e-03;
double electronQuench::m_kBResid  = 0;
double electronQuench::p0 = 1.02561e+00;
double electronQuench::p1 = 1.12245e-01;
double electronQuench::p2 = 1.39421e+00;
double electronQuench::p3 = 5.55117e-04;

double electronQuench::m_quenching_energy[m_nKb][m_nSamples] = {0};
double* electronQuench::m_quenching_energy_low = &m_quenching_energy[0][0];
double* electronQuench::m_quenching_energy_high = &m_quenching_energy[0][0];
double electronQuench::m_quenchingShape1[m_nKb][m_nSamples] = {0};
double* electronQuench::m_quenchingShape1_lowKb = &m_quenchingShape1[0][0];
double* electronQuench::m_quenchingShape1_higKb = &m_quenchingShape1[0][0];

bool electronQuench::m_loadStopPowData = false;
bool electronQuench::m_loadNLData = false;
bool electronQuench::m_loadScintPE = false;

vector<double> electronQuench::m_Etrue;
vector<double> electronQuench::m_StopPow;

double electronQuench::m_edep[1000] = {0.};
double electronQuench::m_nonl[1000] = {0.};

double electronQuench::m_simEtrue[900];
double electronQuench::m_simScintPE[900];
double electronQuench::m_scale = 3134.078/2.223;

TGraph* electronQuench::gNPE_elec = new TGraph();

electronQuench::electronQuench()
{;}

electronQuench::~electronQuench()
{;}

void electronQuench::LoadStopPowData () {
    cout << " >>> Loading Stoppint Power Data <<< " << endl;
    
    ifstream in;
    in.open(junoParameters::stopPow_File.c_str());
    if(!in){
        cout << " >>> Fail to Open Quench File!! <<< " << endl;
        return;
    }

    string line;
    double tmp_E, tmp_StopPow;
    while(getline(in,line))
    {
        istringstream ss(line);
        ss >> tmp_E >> tmp_StopPow;
        m_Etrue.push_back(tmp_E);
        m_StopPow.push_back(tmp_StopPow); 
    }
    in.close();

    m_loadStopPowData = true;
}

double electronQuench::Integral_BirkLaw(double E)
{
    if(!m_loadStopPowData) LoadStopPowData();

    if( m_Etrue.size()==0 || m_StopPow.size() == 0) {
        std::cout << " >>> No Data Loaded in Vector <<< " << endl;
        return -1.0;
    } else if (m_Etrue.size() != m_StopPow.size() ) {
        std::cout << " >>> Quench Integral Vectors have Different Lengths <<< " << std::endl;
        return -1.0;
    } else {
        // numerical integral to get quenched nonlinearity...
        int num = m_Etrue.size();
        double sum = 0.;
        double integ_part1 = 0.; double integ_part2 = 0.;
        for(int i=1; i<num; i++){ //cout << "m_birk1 : " << m_birk1 << endl;
            integ_part1 = 1./(1+m_birk1*m_StopPow[i-1]);
            integ_part2 = 1./(1+m_birk1*m_StopPow[i]);
            if( m_Etrue[i] <= E ){ sum+=(integ_part1+integ_part2)*(m_Etrue[i]-m_Etrue[i-1])/2.; }
            else {break;}
        }
        
        return m_kA*sum/E;
        
    }
    
}

void electronQuench::LoadScintPE() {
    cout << " >>> Load scintillation PE file <<< " << endl;
    ifstream in; in.open(junoParameters::scintFile.c_str());
    if (!in) cout << " >>> No Input File -> " << junoParameters::scintFile << " <<< " << endl;
    string line; double tmpE, tmpPE;
    int index = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmpE >> tmpPE;
        //m_simEtrue[index] = tmpE;
        //m_simScintPE[index] = tmpPE;
        gNPE_elec->SetPoint(index, tmpE, tmpPE);
        index++;
    }
    in.close();

    m_loadScintPE = true;
}


void electronQuench::LoadNLData ()  {
    cout << " >>> Loading Quenching NL Data " << junoParameters::quenchNL_File << " <<< " << endl;
    TFile* quenchingFile = new TFile(junoParameters::quenchNL_File.c_str(), "read"); 
    if(!quenchingFile) { std::cout << " >>> Fail to Open QuenchNL File <<< " << std::endl; }
    for(int kbIdx=51; kbIdx<76; kbIdx++)  {
        //cout << kbIdx << endl;
        stringstream ss; ss << kbIdx;
        TString name1 = "kB"+ss.str();

        //TGraph* quench1G = (TGraph*)quenchingFile->Get(name1);
        TH1D* quench1G = (TH1D*)quenchingFile->Get(name1);
        if(!quench1G) { cout << "No Such A Histogram " << name1 << " in Quench.root File" << endl; return;  }
        //double* quench1 = quench1G->GetY();

        for(int sampleIdx=0; sampleIdx<m_nSamples; sampleIdx++)
        {
			//m_quenchingShape1[kbIdx][sampleIdx] = quench1[sampleIdx];
            m_quenching_energy[kbIdx][sampleIdx] = quench1G->GetBinCenter(sampleIdx+1);
            m_quenchingShape1[kbIdx][sampleIdx] = quench1G->GetBinContent(sampleIdx+1);
        }
        delete quench1G;
    }
	quenchingFile->Close();
    delete quenchingFile;

    m_loadNLData = true;
}

void electronQuench::Update () {
    if(!m_loadNLData) LoadNLData();
    if(m_birk1==0) return;
    else if (m_birk1 > m_birkHigh) m_birk1 = m_birkHigh;
    else if (m_birk1 < m_birkLow ) m_birk1 = m_birkLow;
    int kBIdx = int(m_birk1*1e4);   
	m_kBResid = kBIdx+1 - m_birk1*1e4; 
    m_quenching_energy_low  = &m_quenching_energy[kBIdx] [0];
    m_quenching_energy_high = &m_quenching_energy[kBIdx] [0];
	m_quenchingShape1_lowKb = &m_quenchingShape1[kBIdx]  [0];
	m_quenchingShape1_higKb = &m_quenchingShape1[kBIdx+1][0];
}


double electronQuench::ScintillatorNL    (double eTrue)  {
    
    double nl = ScintillatorShape(eTrue);
    return nl;
}

double electronQuench::ScintillatorPE(double eTrue) {
    if (!m_loadScintPE) LoadScintPE();

    double scintPE = gNPE_elec->Eval(eTrue);
    return m_kA * scintPE;

    double nonl = SimulationNLShape(eTrue);
    return m_kA * nonl * eTrue * m_scale;
}

double electronQuench::ScintillatorShape (double eTrue)  {
    if (junoParameters::scintillatorParameterization == kIntegral) {
        return IntegralNLShape (eTrue);
    } else if (junoParameters::scintillatorParameterization == kSimulation ) { 
        return SimulationNLShape (eTrue);
    } else if (junoParameters::scintillatorParameterization == kEmpirical) {
        return EmpiricalNLShape (eTrue);
    } else if (junoParameters::scintillatorParameterization == kSimulationCalc ) {
        return SimulationNLCalcShape (eTrue);
    }
}

double electronQuench::SimulationNLShape (double eTrue)  {
    if( !m_loadNLData ) LoadNLData();
    Update();
    if (m_birk1 ==0 ) return 1.0;
    int idx = 0;
    if(eTrue<=0.1) { idx = int(eTrue/0.001); }
    else { idx = int((eTrue-0.1)/0.01)+100; }
    
    double quenchNL_low = m_quenchingShape1_lowKb[idx-1] + (m_quenchingShape1_lowKb[idx]-m_quenchingShape1_lowKb[idx-1])*(eTrue-m_quenching_energy_low[idx-1])/ (m_quenching_energy_low[idx]-m_quenching_energy_low[idx-1]);
    double quenchNL_high = m_quenchingShape1_higKb[idx-1] + (m_quenchingShape1_higKb[idx]-m_quenchingShape1_higKb[idx-1])*(eTrue-m_quenching_energy_high[idx-1])/ (m_quenching_energy_high[idx]-m_quenching_energy_high[idx-1]);

    //double quenchNL = m_kA * m_quenchingShape1_lowKb[idx];
    double quenchNL  =  m_kA * ( m_kBResid * quenchNL_low + (1-m_kBResid) * quenchNL_high );
    return quenchNL;
}

double electronQuench::IntegralNLShape   (double eTrue)  {
    if (m_birk1 ==0  ) return 1.0;
    return Integral_BirkLaw(eTrue);
}

double electronQuench::EmpiricalNLShape  (double eTrue) {
    return (p0+p3*eTrue)/(1+p1*TMath::Exp(-p2*eTrue));
}

double electronQuench::SimulationNLCalcShape(double eTrue)
{
    if (!m_loadScintPE) LoadScintPE();

    double deltaE = 0.05; 
    int lowbin  = int(eTrue/0.05);
    int highbin = int(eTrue/0.05) + 1; 
    if ( highbin > 320 ) {
        cout << " >> Energy Beyond Range !!! <<< " << endl;
        return -1;
    }   

    double bias = 1 - (eTrue - m_simEtrue[lowbin]) / deltaE;
    double scintPE = bias * m_simScintPE[lowbin] + (1 - bias) * m_simScintPE[highbin];
    double nonl = scintPE / m_scale / eTrue;

    return m_kA * nonl;
}

void electronQuench::Plot() {
    cout << " >>> Plot Quenching Curve <<< " << endl;
    for(int i=0; i<1000; i++) {
        m_edep[i] = 16./1000.*(i+1);
        m_nonl[i] = electronQuench::ScintillatorNL(m_edep[i]);
    }
    
    TFile* file = new TFile(junoParameters::quenchNL_outFile.c_str(), "recreate");
    TGraph* gQuenchNL = new TGraph(1000, m_edep, m_nonl);
    gQuenchNL->SetLineColor(kBlue+1);
    gQuenchNL->SetMarkerColor(kBlue+1);
    gQuenchNL->SetMarkerSize(0.2);
    gQuenchNL->Write();
    file->Close();

}
