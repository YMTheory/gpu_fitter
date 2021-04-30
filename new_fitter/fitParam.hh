#include <iostream>

#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"

#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"
#include "electronFitter.hh"
#include "gammaData.hh"
#include "gammaFitter.hh"
#include "junoB12.hh"
#include "junoB12_simplified.hh"
#include "gammaResponse.hh"
#include "junoSpectrum.hh"
#include "spectrumFitter.hh"
#include "junoNLChiFunction.hh"

using namespace std;

void SetStyle()
{
	TStyle *style = new TStyle("Modern","Modern Style");
	style->SetTitleFont(43,"xyz");
	style->SetLabelFont(43,"xyz");
	style->SetLegendFont(43);
	style->SetLabelSize(19,"xyz");
	style->SetTitleSize(21,"xyz");
	style->SetLegendBorderSize(0);
	style->SetLegendFillColor(kRed);
	style->SetStatStyle(0);
	style->SetLineColor(kBlue+1);
	style->SetMarkerStyle(20);
	gROOT->SetStyle("Modern"); 
}
