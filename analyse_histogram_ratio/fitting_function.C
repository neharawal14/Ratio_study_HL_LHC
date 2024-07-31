#include "fitting_function.h"
#include "RooMyPDF_BW.h"
#include "RooFFTConvPdf.h"
//void FitFunction::initialise_fit(TString sample_name_string,TString saving_path_string, TString title_factor ){
void FitFunction::initialise_fit(TString sample_name_string,TString saving_path_string){
  //saving_path = saving_path_string+"scale_"+title_factor+"_";
  saving_path = saving_path_string;
  sample_name = sample_name_string;
}
std::pair<float, float> FitFunction :: fitting_BW(TH1D *hist_fit, TString saving_name){
  // Fit the mass into Breit Wigner
  TH1D *hist_fit_clone = (TH1D*) hist_fit->Clone();
  //RooRealVar mass_var("mass_var","mass_var",70,110);
  double mass_var_low; 
  double mass_var_up;  
 
  std::cout<<"name of the histogram fit "<<hist_fit->GetName()<<std::endl;
  std::cout<<"entries of the histogram fit "<<hist_fit->GetEntries()<<std::endl;

  double mean_low, mean_up, mean_value; 
  double width_low, width_up, width_value; 
  double fit_range_low, fit_range_up;
  if(sample_name=="isHiggs"){
  mass_var_low= 124.9; 
  mass_var_up= 125.1; 
  mean_low = 124.95;
  mean_up = 125.05;
  mean_value = 125; 
  width_low = 0.002;
  width_up = 2;
  width_value = 0.004 ; 
  fit_range_low = 124.996;
  fit_range_up = 125.004;

  }
  else{
  mass_var_low= 70; 
  mass_var_up= 110; 
  mean_low = 88;
  mean_up = 94;
  mean_value = 91; 
  width_low = 1;
  width_up = 3;
  width_value = 2.49 ; 
//  fit_range_low = 86;
//  fit_range_up = 96;

  fit_range_low = 89;
  fit_range_up = 93.5;

  }


  RooRealVar mass_var("mass_var","mass_var", mass_var_low, mass_var_up);
  RooDataHist histo("histo","mass dataset",mass_var,hist_fit_clone);
  RooRealVar mean_mass("mean_mass","Z mass",mean_value, mean_low, mean_up);
  RooRealVar width("width","width of Z mass",width_value, width_low, width_up);
  
  RooMyPDF_BW BW("BW","Breit Wigner fit",mass_var, mean_mass,width);
  RooPlot *xframe=mass_var.frame();
  histo.plotOn(xframe);
  std::cout<<" histo defined "<<std::endl;
  BW.fitTo(histo,Range(fit_range_low, fit_range_up));
  std::cout<<" histo fitted "<<std::endl;
  BW.plotOn(xframe,RooFit::LineColor(kRed+2),Name("BW_sig"));
  BW.paramOn(xframe,RooFit::Layout(0.6,0.9,0.7));
     
   TCanvas *tmp = new TCanvas("tmp","Gen Z mass", 900,600);
   tmp->cd();
   gPad->SetLeftMargin(0.15);
   xframe->getAttText()->SetTextSize(0.025);
   xframe->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
   xframe->SetTitle(hist_fit->GetTitle()); 
//   xframe->GetYaxis()->SetTitle("N/0.5 (GeV)");
   xframe->Draw();
  float chi_square_value = xframe->chiSquare();
  int underflow = hist_fit_clone->GetBinContent(0);
  int overflow = hist_fit_clone->GetBinContent(hist_fit_clone->GetNbinsX() + 1);
  std::cout << "Underflow: " << underflow << std::endl;
  std::cout << "Overflow: " << overflow << std::endl;

  std::pair<float, float> mean_pair;
  mean_pair.first = mean_mass.getVal();
  mean_pair.second = mean_mass.getError();
  TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kBlack);
  leg2->AddEntry("histo","Gen Z", "EP");
  leg2->AddEntry("BW_sig","BW fit","LP");
  leg2->AddEntry("xframe->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
  leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
  leg2->Draw("same");       
  leg2->AddEntry("Underflow",Form("Underflow= %d", underflow),"");
  leg2->AddEntry("Overflow",Form("Overflow= %d", overflow),"");

  gStyle->SetOptStat();
  tmp->SaveAs(saving_path+saving_name+"_BW_fit.pdf");

  TCanvas *c_MC_2 = new TCanvas("c_MC_2", "c_MC_2", 1000, 600);
  c_MC_2->SetFrameFillColor(0);
  //RooPlot* xframe1 = mass_var.frame(RooFit::Title(title_name));
  //histo.plotOn(xframe1);
  c_MC_2->cd();
  c_MC_2->SetTopMargin(1.5);
  xframe->GetYaxis()->SetTitle("");
  xframe->GetXaxis()->SetTitle("");
  xframe->chiSquare();
  xframe->Draw();
  gPad->SetLogy(); 
  gStyle->SetOptStat();
  TString saving_name_2=saving_path+saving_name+"_BW_log.pdf";
  c_MC_2->SaveAs(saving_name_2);
  c_MC_2->Close();

  return mean_pair;
    }
std::pair<float, float> FitFunction :: fitting_DSCB(TH1D *hist, TString saving_name,
    float scale_factor, float fit_low_value, float fit_up_value){
      std::pair<float, float> mean_pair;
      std::pair<float, float> sigma_pair;

     TH1F * hist_clone = (TH1F*) hist->Clone();
   
    double mass_var_low , mass_var_up;
    double cb_mean_low , cb_mean_up , cb_mean;
    double width_mean_low , width_mean_up, width_mean;
    double fit_range_low,   fit_range_up ;
    double alphaL_value, alphaR_value;
    double nL_value, nR_value;
     double nL_min , nL_max ;
     double nR_min , nR_max ;
     double alphaL_min , alphaL_max ;
     double alphaR_min , alphaR_max ;

    double tolerance = 1e-5; // Define a tolerance for floating-point comparison

    if(sample_name=="isHiggs"){
     mass_var_low = 120; 
     mass_var_up = 130; 
     cb_mean = 125;
     width_mean = 2;
     cb_mean_low = 124.2; 
     cb_mean_up = 127.2; 

//     cb_mean_low = 124.2; 
//     cb_mean_up = 125.8; 
     width_mean_low = 0.01; 
     width_mean_up = 5; 
     nL_value = 3;
     nR_value = 3; 
     alphaL_value = 0.8; 
     alphaR_value = 0.8; 
     nL_min = 0.1; nL_max = 20;
     nR_min = 0.1; nR_max = 20;
     alphaL_min = 0.01; alphaL_max = 5;
     alphaR_min = 0.01; alphaR_max = 5;
     //fit_range_low = 123.52; 
     //fit_range_up = 126.5; 
    
     std::cout<<" scale factor :"<<scale_factor<<std::endl; 
     //if(fabs(scale_factor - 0) < tolerance){

     fit_range_low = fit_low_value;
     fit_range_up = fit_up_value;
//     if(scale_factor < 0.005){
//       std::cout<<" scale factor new"<<scale_factor<<std::endl; 
//       fit_range_low = 123.5; 
//       fit_range_up = 126.3; 
//       }
//     else {
//       fit_range_low = 123.4; 
//       fit_range_up = 128.2; 
//
//     }
//     if(fabs(scale_factor - 0.01) < tolerance){
//       std::cout<<" scale factor new"<<scale_factor<<std::endl; 
//       fit_range_low = 123.5; 
//       fit_range_up = 126.3; 
//     }
//     if(fabs(scale_factor - 0.02) < tolerance){
//       std::cout<<" scale factor new"<<scale_factor<<std::endl; 
//       fit_range_low = 123.5; 
//       fit_range_up = 126.4; 
//     }
//
//     if(fabs(scale_factor - 0.03) < tolerance){
//       std::cout<<" scale factor new"<<scale_factor<<std::endl; 
//       fit_range_low = 123.6; 
//       fit_range_up = 126.4; 
//     }
//     if(fabs(scale_factor - 0.04) < tolerance){
//       std::cout<<" scale factor new"<<scale_factor<<std::endl; 
//       fit_range_low = 123.6; 
//       fit_range_up = 126.4; 
//     }
//
//    if( 0.05 <= scale_factor && scale_factor < 0.1){
//    // if(fabs(scale_factor - 0.05) < tolerance){
//       std::cout<<" scale factor new"<<scale_factor<<std::endl; 
//       fit_range_low = 123.6; 
//       fit_range_up = 126.5; 
//     }
//    if( 0.1 <= scale_factor){
//       std::cout<<" scale factor new"<<scale_factor<<std::endl; 
//       fit_range_low = 123.6; 
//       fit_range_up = 126.2; 
//     }

     // fit_range_low = 123.54; 
    // fit_range_up = 126.4; 
    }
    else{
     mass_var_low = 70; 
     mass_var_up = 110; 
     cb_mean = 91;
     cb_mean_low = 88; 
     cb_mean_up = 94;
     width_mean = 2.49; 
     width_mean_low = 0.005; 
     width_mean_up = 8; 

     alphaL_value = 1; 
     alphaR_value = 1; 
     nL_min = 0.1; nL_max = 20;
     nR_min = 0.1; nR_max = 20;
     alphaL_min = 0.01; alphaL_max = 5;
     alphaR_min = 0.01; alphaR_max = 5;

     nL_value = 3;
     nR_value = 3; 
     std::cout<<" scale factor Z :"<<scale_factor<<std::endl; 
    if(fabs(scale_factor - 0) < tolerance){
     std::cout<<" scale factor Z new:"<<scale_factor<<std::endl; 
     fit_range_low = 89; 
     fit_range_up = 93; 
     }
    if(fabs(scale_factor - 0.01) < tolerance){
     std::cout<<" scale factor Z new:"<<scale_factor<<std::endl; 
     fit_range_low = 89; 
     fit_range_up = 93; 
     }
    if(fabs(scale_factor - 0.02) < tolerance){
     std::cout<<" scale factor Z new:"<<scale_factor<<std::endl; 
     fit_range_low = 89; 
     fit_range_up = 93; 
     }

    if(fabs(scale_factor - 0.03) < tolerance){
     std::cout<<" scale factor Z new:"<<scale_factor<<std::endl; 
     fit_range_low = 89; 
     fit_range_up = 93.5; 
     }

    if(fabs(scale_factor - 0.04) < tolerance){
     std::cout<<" scale factor Z new:"<<scale_factor<<std::endl; 
     fit_range_low = 89; 
     fit_range_up = 94; 
     }

    if( 0.05 <=  scale_factor && scale_factor < 0.1){
    //if(fabs(scale_factor - 0.05) < tolerance){
     std::cout<<" scale factor Z new:"<<scale_factor<<std::endl; 
     fit_range_low = 89; 
     fit_range_up = 94; 
     }
    if(fabs(scale_factor - 0.1) < tolerance){
     std::cout<<" scale factor Z new:"<<scale_factor<<std::endl; 
     fit_range_low = 89; 
     fit_range_up = 94; 
     }

    if(scale_factor > 0.1){
     std::cout<<" scale factor Z new:"<<scale_factor<<std::endl; 
     fit_range_low = 89; 
     fit_range_up = 93.5; 

    }
     //fit_range_low = 89; 
     //fit_range_up = 94; 

//     fit_range_low = 89; 
//     fit_range_up = 93; 
    }
     RooRealVar mass_var("mass_var","mass_var",mass_var_low,mass_var_up);
     //RooDataHist histo_gen("histo_gen","mass dataset",mass_var,gen_hist);
     //RooRealVar mean_mass("mean_mass","mean of Z mass",91.04,88,94);
     //RooRealVar width("width","width of Z mass",2.49,1,5);
     //RooMyPDF_BW BW("BW","Breit Wigner fit",mass_var, mean_mass, width);
     //BW.fitTo(histo_gen,Range(88,94));
    
     RooDataHist histo_reco("histo_reco","mass dataset",mass_var,hist_clone);
     //std::cout<<" hist name :::::: "<<gen_hist_clone->GetName()<<std::endl;
     //TString name = gen_hist_clone->GetName();

    RooRealVar cbmean("cbmean", "cbmean" , cb_mean, cb_mean_low, cb_mean_up) ;
    RooRealVar cbsigmaLR("cbsigmaLR", "cbsigmaLR" , width_mean, width_mean_low, width_mean_up) ;
    RooRealVar cbsigmaL("cbsigmaL", "cbsigmaL" , width_mean, width_mean_low, width_mean_up) ;
    RooRealVar cbsigmaR("cbsigmaR", "cbsigmaR" , width_mean, width_mean_low, width_mean_up) ;

//    RooRealVar cbsigmaL("cbsigmaL", "cbsigmaL" , 0.5, 0.01, 10) ;
//    RooRealVar cbsigmaR("cbsigmaR", "cbsigmaR" , 0.5, 0.01,10) ;
    RooRealVar alphaL("alphaL","alphaL", alphaL_value, alphaL_min, alphaL_max);
    RooRealVar alphaR("alphaR","alphaR", alphaR_value, alphaR_min, alphaR_max);
    RooRealVar nL("nL", "nL", nL_value, nL_min, nL_max);
    RooRealVar nR("nR", "nR", nR_value, nR_min, nR_max);
   
    //RooRealVar alphaL("alphaL","alphaL", 5, alphaL_min, alphaL_max);
    //RooRealVar alphaR("alphaR","alphaR", 5, alphaR_min, alphaR_max);
    //RooRealVar nL("nL", "nL", 15, nL_min, nL_max);
    //RooRealVar nR("nR", "nR", 15, nR_min, nR_max);
 //   if(name=="gen_Zmass"){
 //     std::cout<<" ***********************************yes indeed come  here *******************************"<<std::endl;
 //     alphaR_max = 15;
 //     alphaL_max = 7;
 //     nR_min = 1;
 //     nR_max = 5;
 //   } 

//     RooCrystalBall DSCB("DSCB", "DSCB", mass_var, cbmean, cbsigmaLR, alphaL, nL, alphaR, nR);
     RooCrystalBall DSCB("DSCB", "DSCB", mass_var, cbmean, cbsigmaL, cbsigmaR, alphaL, nL, alphaR, nR);
     // making convolution fit 
     //mass_var.setBins(50000, "cache");
     //RooFFTConvPdf BW_DSCB("BW_DSCB", "BW (X) DSCB", mass_var, BW, DSCB);
     // Setting mean and width as constant from BW fit 
     //mean_mass.setConstant(true); 
     //width.setConstant(true); 

     DSCB.fitTo(histo_reco, Range(fit_range_low,fit_range_up));

     std::cout << "***************************************fitted**************************"<< std::endl;
     // Save the fit
     TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
     c_MC->SetFrameFillColor(0);
     c_MC->cd();
     //c_MC->SetTopMargin(1.5);
     RooPlot *frame = mass_var.frame();
     histo_reco.plotOn(frame);
     DSCB.plotOn(frame, LineColor(kBlue), Name("DSCB"));
     DSCB.paramOn(frame, Layout(0.12, 0.7, 0.7));
     
     frame->getAttText()->SetTextSize(0.04);
     frame->SetTitle(hist_clone->GetTitle());
     frame->GetXaxis()->SetTitle("GeV");
     frame->chiSquare();
     frame->Draw();
     double chi_square_value = frame->chiSquare();
     int underflow = hist_clone->GetBinContent(0);
     int overflow = hist_clone->GetBinContent(hist_clone->GetNbinsX() + 1);
     std::cout << "Underflow: " << underflow << std::endl;
     std::cout << "Overflow: " << overflow << std::endl;

     TLegend* leg2 = new TLegend(0.7, 0.7, 0.89, 0.89);
     leg2->AddEntry("DSCB","DSCB fit","l");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",chi_square_value),"");
     leg2->AddEntry("histo_reco.sumEntries()",Form("Events= %.0f",histo_reco.sumEntries()),"");
     leg2->AddEntry("Underflow",Form("Underflow= %d", underflow),"");
     leg2->AddEntry("Overflow",Form("Overflow= %d", overflow),"");
     leg2->Draw("same");

     gStyle->SetOptStat();
     c_MC->SaveAs(saving_path+saving_name+"_DSCB_fit.pdf");
     c_MC->Close();
     
     mean_pair.first   = cbmean.getVal();
     mean_pair.second  = cbmean.getError();
     //sigma_pair.first  = Sigma.getVal();
     //sigma_pair.second = Sigma.getError();

     // Plot log plot also for the fit
     TCanvas *c_MC_2 = new TCanvas("c_MC_2", "c_MC_2", 1000, 600);
     c_MC_2->SetFrameFillColor(0);
     c_MC_2->cd();
     c_MC_2->SetTopMargin(1.5);
     frame->GetYaxis()->SetTitle("");
     frame->GetXaxis()->SetTitle("");
     frame->chiSquare();
     frame->Draw();
     leg2->Draw("same");
     gPad->SetLogy(); 
     gStyle->SetOptStat();
     TString saving_name_2=saving_path+saving_name+"_DSCB_log.pdf";
     c_MC_2->SaveAs(saving_name_2);
     c_MC_2->Close();
    
     return mean_pair;
    }


  std::pair<float, float> FitFunction :: fitting_BW_DSCB(TH1D *hist_gen, TH1D*hist, TString saving_name, float scale_factor){
      std::pair<float, float> mean_pair;
      std::pair<float, float> sigma_pair;

     TH1F * hist_clone = (TH1F*) hist->Clone();
   
    double mass_var_low , mass_var_up;
    double cb_mean_low , cb_mean_up , cb_mean;
    double fit_range_low,   fit_range_up ;

    double BW_mean, BW_mean_low, BW_mean_up;
    double width_mean, width_low, width_up;
    double sigma_mean, sigma_mean_low, sigma_mean_up;
    double BW_range_low, BW_range_up;
    
    double alphaL_value, alphaL_min, alphaL_max;
    double alphaR_value, alphaR_min, alphaR_max;
    double nL_value, nL_min, nL_max;
    double nR_value, nR_min, nR_max;
    if(sample_name =="isHiggs"){
     mass_var_low = 120; 
     mass_var_up = 130; 
//     cb_width_mean = 0.1;
//     cb_width_mean_low = 0.0001; 
//     cb_width_mean_up = 2; 
     nL_value = 3;
     nR_value = 3; 
     alphaL_value = 0.8; 
     alphaR_value = 0.8; 
     nL_min = 0.1; nL_max = 20;
     nR_min = 0.1; nR_max = 20;
     alphaL_min = 0.01; alphaL_max = 5;
     alphaR_min = 0.01; alphaR_max = 5;

     mass_var_low = 120; 
     mass_var_up = 130; 
     cb_mean = 0;
     fit_range_low = 120; 
     fit_range_up = 130; 

     BW_mean = 125; 
     BW_mean_low = 124.9; 
     BW_mean_up = 125.1; 
     width_mean = 0.004; 
     width_low = 0.001;
     width_up = 2;
     BW_range_low = 124.9;
     BW_range_up = 125.1;
     sigma_mean_low = 0.001; 
     sigma_mean_up = 1;
     sigma_mean= 0.004; 

    }
    else{
     mass_var_low = 70; 
     mass_var_up = 110; 
     cb_mean = 0;
     
     mass_var_low = 70; 
     mass_var_up = 110; 

     alphaL_value = 1; 
     alphaR_value = 1; 
     nL_min = 0.1; nL_max = 20;
     nR_min = 0.1; nR_max = 20;
     alphaL_min = 0.01; alphaL_max = 5;
     alphaR_min = 0.01; alphaR_max = 5;

     nL_value = 3;
     nR_value = 3; 
 
     BW_mean = 91.04; 
     BW_mean_low = 88; 
     BW_mean_up = 94; 
     sigma_mean_low = 0.001; 
     sigma_mean_up = 1.5;
     sigma_mean= 0.02; 

     width_mean = 2.49;
     width_low = 1;
     width_up = 5;
     BW_range_low = 89;
     BW_range_up = 93.5;

     fit_range_low = 83; 
     fit_range_up = 98; 
    }
     RooRealVar mass_var("mass_var","mass_var",mass_var_low,mass_var_up);

     RooDataHist histo_gen("histo_gen","mass dataset",mass_var,hist_gen);
     
     RooRealVar mean_mass("mean_mass","mean of Z mass",BW_mean, BW_mean_low, BW_mean_up);
     RooRealVar width("width","width of Z mass",width_mean, width_low, width_up);
     RooMyPDF_BW BW("BW","Breit Wigner fit",mass_var, mean_mass, width);
     BW.fitTo(histo_gen,Range(BW_range_low, BW_range_up));
    
     RooDataHist histo_reco("histo_reco","mass dataset",mass_var,hist_clone);

    RooRealVar cbmean("cbmean", "cbmean" , cb_mean, -1,1) ;
    RooRealVar cbsigmaLR("cbsigmaL", "cbsigmaL" , sigma_mean, sigma_mean_low, sigma_mean_up) ;
//    RooRealVar cbsigmaR("cbsigmaR", "cbsigmaR" , sigma_mean, sigma_mean_low, sigma_mean_up) ;

//    RooRealVar cbsigmaL("cbsigmaL", "cbsigmaL" , 0.5, 0.01, 10) ;
//    RooRealVar cbsigmaR("cbsigmaR", "cbsigmaR" , 0.5, 0.01,10) ;
    //RooRealVar alphaL("alphaL","alphaL", 1, 0.01, 7);
    //RooRealVar alphaR("alphaR","alphaR", 4, 0.01, 15);
    //RooRealVar nL("nL", "nL", 5, 0.1, 100);
    //RooRealVar nR("nR", "nR", 1, 0.1, 100);
    
    RooRealVar alphaL("alphaL","alphaL", alphaL_value, alphaL_min, alphaL_max);
    RooRealVar alphaR("alphaR","alphaR", alphaR_value, alphaR_min, alphaR_max);
    RooRealVar nL("nL", "nL", nL_value, nL_min, nL_max);
    RooRealVar nR("nR", "nR", nR_value, nR_min, nR_max);

    RooCrystalBall DSCB("DSCB", "DSCB", mass_var, cbmean, cbsigmaLR, alphaL, nL, alphaR, nR);
     // making convolution fit 
     mass_var.setBins(50000, "cache");
     RooFFTConvPdf BW_DSCB("BW_DSCB", "BW (X) DSCB", mass_var, BW, DSCB);
     // Setting mean and width as constant from BW fit
     double width_value = width.getVal();  
     width.setConstant(true); 
     cbmean.setConstant(true);
     BW_DSCB.fitTo(histo_reco, Range(fit_range_low,fit_range_up));

     std::cout << "***************************************fitted**************************"<< std::endl;
     // Save the fit
     TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
     c_MC->SetFrameFillColor(0);
     c_MC->cd();
     //c_MC->SetTopMargin(1.5);
     RooPlot *frame = mass_var.frame();
     histo_reco.plotOn(frame);
     BW_DSCB.plotOn(frame, LineColor(kBlue), Name("BW_DSCB"));
     //BW_DSCB.paramOn(frame,Layout(0.12, 0.4, 0.6));
     BW.plotOn(frame, LineColor(kRed+2), LineStyle(kDashed), Name("BW_sig"));
     BW.paramOn(frame, Layout(0.12, 0.4, 0.68));
     DSCB.plotOn(frame, LineColor(kGreen), LineStyle(kDashed), Name("DSCB_sig"));
     DSCB.paramOn(frame, Layout(0.12, 0.4, 0.55));
     
     frame->getAttText()->SetTextSize(0.04);
     frame->GetYaxis()->SetTitle("");
     frame->GetXaxis()->SetTitle("GeV");
     frame->SetTitle(hist->GetTitle());
     frame->chiSquare();
     frame->Draw();
     double chi_square_value = frame->chiSquare("BW_DSCB","h_histo_reco",6);;
     int underflow = hist_clone->GetBinContent(0);
     int overflow = hist_clone->GetBinContent(hist_clone->GetNbinsX() + 1);
     std::cout << "Underflow: " << underflow << std::endl;
     std::cout << "Overflow: " << overflow << std::endl;

     TLegend* leg = new TLegend(0.72, 0.7, 0.9, 0.89);
     leg->AddEntry("width",Form("BW width = %.5f GeV", width_value),"");
     TLegend* leg2 = new TLegend(0.12, 0.7, 0.5, 0.89);
     leg2->AddEntry("BW_DSCB","BW DSCB fit","l");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",chi_square_value),"");
     leg2->AddEntry("histo_reco.sumEntries()",Form("Events= %.0f",histo_reco.sumEntries()),"");
     leg2->AddEntry("Underflow",Form("Underflow= %d", underflow),"");
     leg2->AddEntry("Overflow",Form("Overflow= %d", overflow),"");
     leg2->Draw("same");
     leg->Draw("same");

     gStyle->SetOptStat();
     c_MC->SaveAs(saving_path+saving_name+"_BW_DSCB_fit.pdf");
     c_MC->Close();
     
     mean_pair.first   = mean_mass.getVal();
     mean_pair.second  = mean_mass.getError();
//     sigma_pair.first  = Sigma.getVal();
//     sigma_pair.second = Sigma.getError();

     // Plot log plot also for the fit
     TCanvas *c_MC_2 = new TCanvas("c_MC_2", "c_MC_2", 1000, 600);
     c_MC_2->SetFrameFillColor(0);
     c_MC_2->cd();
     c_MC_2->SetTopMargin(1.5);
     frame->GetYaxis()->SetTitle("");
     frame->GetXaxis()->SetTitle("");
     frame->chiSquare();
     frame->Draw();
     leg2->Draw("same");
     gPad->SetLogy(); 
     gStyle->SetOptStat();
     TString saving_name_2=saving_path+saving_name+"_BW_DSCB_log.pdf";
     c_MC_2->SaveAs(saving_name_2);
     c_MC_2->Close();
     return mean_pair;
    }



///std::pair<float, float> fitting_delta(TH1D *hist_fit, TString saving_name){
///      //TGaxis::SetMaxDigits(20);
///      // delcaring pair to store mean and sigma values after fit to DSCB
///      //std::pair<float, float> mean_pair;
///      //std::pair<float, float> sigma_pair;
///      TH1F * hist_fit_clone = (TH1F*) hist_fit->Clone();
///      // define a initial gaussian fit to find the initial mean and sigma parameters for DSCB fit
///      TF1* f1;
///			double mean_hist, sigma_hist;
///      TString title_hist;
///      mean_hist = hist_fit->GetMean();
///      sigma_hist = hist_fit->GetRMS();
///      title_hist = hist_fit->GetTitle();
///
///      double range_low_x = hist_fit->GetBinLowEdge(0);
///      double range_up_x = hist_fit->GetBinLowEdge(hist_fit->GetNbinsX()+1);
///      double low_x = mean_hist - 3*sigma_hist;
///      double up_x = mean_hist + 3*sigma_hist;
///      
///      hist_fit_clone->GetXaxis()->SetRangeUser(low_x, up_x);
///      f1 = new TF1("f1", "gaus", low_x,up_x);
///      hist_fit_clone->Fit(f1,"R");
///      float mean_gaus, sigma_gaus;
///      mean_gaus = f1->GetParameter(1);
///      sigma_gaus = f1->GetParameter(2);
///
/////      std::cout<<"mean and sigma of gaussian *************"<<mean_gaus<<" sigma "<<sigma_gaus<<std::endl;
///      RooRealVar delta_pT("delta_pT", "deltapT", low_x, up_x, "");
/////      RooRealVar delta_pT("delta_pT", "deltapT", range_low_x, range_up_x, "");
///      RooRealVar Mean("Mean", "Mean",mean_gaus, low_x,up_x);
///      RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.005,2.0);//sigma[decay]);
///
///      RooDataHist histo("histo","dataset with var", delta_pT, hist_fit_clone);
///      RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 1.0, 0.05, 10);//alphaL[decay]);
///      RooRealVar ExpL("n_{L}", "n_{L}", 10, 1, 50);//expL[decay]);
///      RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 1.0, 0.05, 10);//alphaR[decay]);
///      RooRealVar ExpR("n_{R}", "n_{R}", 10, 1, 60);//expR[decay]);
///      RooCrystalBall DSCB("DSCB", "DSCB", delta_pT, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
/// 
///      Int_t color = kRed+2;
///      Double_t size_text = 0.020;
///      // Fit the mass into DSCB
///      TCanvas *c_MC_1 = new TCanvas();
///      c_MC_1->SetFrameFillColor(0);
///      TString title_name = hist_fit_clone->GetTitle();
///      RooPlot *xframe_2=delta_pT.frame(Title(title_name));
///      histo.plotOn(xframe_2);
///      DSCB.fitTo(histo, RooFit::Range(low_x, up_x));
///      DSCB.plotOn(xframe_2, RooFit::LineColor(color),Name("fit"));
///      DSCB.paramOn(xframe_2, RooFit::Layout(0.15, 0.35, 0.90));
///      c_MC_1->cd();
///      c_MC_1->SetTopMargin(1.5);
///      xframe_2->getAttText()->SetTextSize(size_text);
///      xframe_2->getAttText()->SetTextColor(color);
///      xframe_2->GetYaxis()->SetTitle("");
///      xframe_2->GetXaxis()->SetTitle("");
///      xframe_2->chiSquare();
///      double chi_square_value = xframe_2->chiSquare();
///      xframe_2->Draw();
///      int underflow = hist_fit_clone->GetBinContent(0);
///      int overflow = hist_fit_clone->GetBinContent(hist_fit_clone->GetNbinsX() + 1);
///
///      TLegend* leg2_2 = new TLegend(0.68, 0.68, 0.89, 0.89);
///      leg2_2->SetFillColor(kWhite);
///      leg2_2->SetLineColor(kBlack);
///      leg2_2->AddEntry("histo","", "EP");
///      leg2_2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",chi_square_value),"");
///      leg2_2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
///      leg2_2->AddEntry("Underflow", Form("underflow = %d",underflow),"");
///      leg2_2->AddEntry("Overflow", Form("overflow = %d",overflow),"");
///      leg2_2->Draw("same");
///      gStyle->SetOptStat();
///      c_MC_1->SaveAs(saving_path+saving_name+"_DSCB_fit.pdf");// + ".pdf");
///      c_MC_1->Close();
///      std::pair<float, float> sigma_pair;
///     sigma_pair.first  = Sigma.getVal();
///     sigma_pair.second = Sigma.getError();
///     return sigma_pair;
///    }
///
///
///void plot_hist(TH1D *hist,TString title, TString saving_name){
///
///  TCanvas *c = new TCanvas();
///  c->cd();
///  gStyle->SetOptStat(111112211);
///  hist->Draw();
///  hist->SetTitle(title);
///  c->SaveAs(saving_path+saving_name);
///
///  // Plotting smae plot with 5 sigma range
///  TH1D *h_clone = (TH1D*) hist->Clone();
///
///  double mean = h_clone->GetMean();
///  double sigma = h_clone->GetRMS();
///  double x_low = mean-3*sigma; 
///  double x_up = mean+3*sigma; 
///
///  std::cout<<" mean "<<mean<<" sigma "<<sigma<<" x low "<<x_low<<" x up "<<x_up<<std::endl;
///  h_clone->GetXaxis()->SetRangeUser(x_low, x_up);
///  TCanvas *c1 = new TCanvas();
///  c1->cd();
///  gStyle->SetOptStat(111112211);
///  h_clone->Draw();
///  c1->SaveAs(saving_path+"range_"+saving_name);
///
///}
