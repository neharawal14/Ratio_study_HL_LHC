#include <iostream>
#include <sstream>
#include "TString.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "RooRealVar.h"
#include "RooCrystalBall.h"
#include "RooFitResult.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TF1.h"
using namespace std;
using namespace RooFit;

class FitFunction{
  public: 
    //void initialise_fit(TString, TString, TString);
    void initialise_fit(TString, TString);
    std::pair<float,float> fitting_BW(TH1D*, TString);
    std::pair<float,float> fitting_DSCB(TH1D* hist, TString, float, float, float);
    std::pair<float, float> fitting_BW_DSCB(TH1D*,TH1D *hist_fit, TString, float );
    TString saving_path;
    TString sample_name;

};
