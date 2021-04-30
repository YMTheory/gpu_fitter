#include "fitParam.hh"

int main()
{
    //SetStyle();

    // electron fake nonl fitting :
    //electronFitter* efitter = new electronFitter();
    //efitter->Initialize();
    //efitter->Minimization();
    //efitter->Plot();

    // gamma calibration source
    //gammaFitter*  gfitter = new gammaFitter();
    //gfitter->Initialize();
    //gfitter->Minimization();
    //gfitter->Plot();
    
    //junoSpectrum* B12 = new junoSpectrum(14400, 100,
    //                                     3, 1,
    //                                     0, 15,
    //                                     4, 14,
    //                                     "B12");
    //B12->LoadData();
    //B12->GetChi2();
    //B12->Plot();

    // beta spectum
    //spectrumFitter* fitter = new spectrumFitter();
    //fitter->Initialize();
    //fitter->Minimization();

    junoNLChiFunction* fitter = new junoNLChiFunction();
    fitter->LoadData();
    fitter->GetChiSquare();
    fitter->Plot();
    //

    //electronResponse::FitPlot();

    //junoSpectrum* junoB12data = new junoSpectrum(14400, 100, 3, 2,
    //                         0, 15, 0, 15, "histogram", "B12");
    //junoB12data->LoadData();
    //cout << junoB12data->GetChi2() <<endl;
    //junoB12data->Plot();

    //junoB12_simplified* b12 = new junoB12_simplified(100, 4000, 18000);
    //b12->Initialize();
    ////b12->LoadDataSpec();
    ////b12->LoadTheoSpec();
    //cout << b12->GetChi2() << endl;
    //b12->Plot();

    //gammaData* cs137 = new gammaData("nFe56", 700, 1100, 100);
    //cs137->LoadData();
    //cs137->calcGammaResponse();
    //cout << cs137->GetChi2();
    
    //gammaResponse* gamResArr[9];
    //gammaResponse* Cs137 = new gammaResponse("Cs137", 200, 600, 1000);  gamResArr[0] = Cs137;
    //gammaResponse* Mn54 = new gammaResponse("Mn54", 200, 900, 1300);    gamResArr[1] = Mn54;
    //gammaResponse* Ge68 = new gammaResponse("Ge68", 200, 1100, 1500);   gamResArr[2] = Ge68;
    //gammaResponse* K40 = new gammaResponse("K40", 200, 1750, 2250);     gamResArr[3] = K40;
    //gammaResponse* nH = new gammaResponse("nH", 200, 2800, 3500);       gamResArr[4] = nH;
    //gammaResponse* Co60 = new gammaResponse("Co60", 100, 3200, 3700);   gamResArr[5] = Co60;
    //gammaResponse* AmBe = new gammaResponse("AmBe", 100, 6000, 6800);   gamResArr[6] = AmBe;
    //gammaResponse* nC12 = new gammaResponse("nC12", 100, 6700, 7600);   gamResArr[7] = nC12;
    //gammaResponse* AmC = new gammaResponse("AmC", 100, 8400, 9400);     gamResArr[8] = AmC;
    //for (int i=0; i<1; i++) {
    //    gamResArr[i]->LoadData();
    //    gamResArr[i]->GetChi2();
    ////    cout << gamResArr[i]->GetName() << " " << gamResArr[i]->GetNonlData() << " " << gamResArr[i]->GetNonlCalc() << " " 
    ////         << gamResArr[i]->GetResData()  << " " << gamResArr[i]->GetResCalc() << endl;
    //    gamResArr[i]->SaveHist();
    //}



    return 1.0;
}
