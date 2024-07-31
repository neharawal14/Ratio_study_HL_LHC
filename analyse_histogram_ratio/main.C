#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "fitting_function.h"
#include "RooMyPDF_BW.h"
#include "TGraphErrors.h"
int main(int argc, char *argv[]){
  TString type_eta = argv[1];
  TString corr = argv[2];
  TString type_sample = argv[3];
  TString output_file_name;
//  TString type_eta = "eta_inclusive"; 
//  TString corr = "noRoch";
//  TString type_sample = "4l";
  TString input_file_name, saving_path;
  TString title;
  if(type_eta=="eta_greater_1"){
   title = "|#eta| >=1";
  }
  else if(type_eta=="eta_lesser_1"){
   title = "|#eta| <=1";
  }
  else if(type_eta=="eta_inclusive"){
   title = "|#eta| <=2.4";
  }


  if(type_sample=="4l"){
    input_file_name = "/eos/home-n/nrawal/RatioStudy/Ratio_taking_Higgs_Z/histogram_ratio_"+corr+"_"+type_eta+".root";
    saving_path = "../results_4l_"+corr+"/"+type_eta+"/";
    output_file_name = "output_graphs_"+corr+"_"+type_eta+"_"+type_sample+".root";
  }
  else if(type_sample=="2l"){
   input_file_name = "/eos/home-n/nrawal/RatioStudy/Ratio_taking_Higgs_Z/histogram_ratio_Z_2l_"+type_eta+"_"+corr+".root";
   saving_path = "../results_2l_"+corr+"/"+type_eta+"/";
  output_file_name = "output_graphs_"+corr+"_"+type_eta+"_"+type_sample+".root";
  }
  TFile *input_file = TFile::Open(input_file_name);

  std::vector<float> fit_low, fit_up;
  if(type_eta =="eta_greater_1" && corr == "noRoch"){
    std::cout<<" here in eta greater 1 "<<std::endl;
//     fit_low = {123.5, 123.2, 123.6, 123.5, 123.6};
//     fit_up = {126.3, 126.4, 126.1 ,126.3, 126.6};
     fit_low = {122.8, 122.8, 123.0, 123.0, 123.2};
     fit_up = {126.8, 126.8, 126.8 ,126.8, 127.6};


  }
  else if(type_eta =="eta_lesser_1" && corr == "noRoch"){
    std::cout<<" here in eta lesser 1 "<<std::endl;
     fit_low = {123.5, 123.5, 123.5, 123.5, 123.6};
     fit_up = {126.3, 126.3, 126.3 ,126.6, 127.4};

  }
  else if(type_eta =="eta_inclusive" && corr == "noRoch"){
    std::cout<<" here in eta inclusive "<<std::endl;
//     fit_low = {122.8, 123.5, 123.0, 123.0, 124.5};
//     fit_up = {126.8, 126.3, 126.5 ,126.5, 127.6};
    fit_low = {123.5, 123.5, 123.5, 123.0, 124.5};
     fit_up = {125.8, 126.3, 126.5 ,126.5, 127.6};

  }

  else if(type_eta =="eta_greater_1" && corr == "Roch"){
    std::cout<<" here in eta greater 1 "<<std::endl;
     fit_low = {123.5, 123.2, 123.6, 123.5, 123.2};
     fit_up = {126.3, 126.4, 126.1 ,126.3, 127.1};
  }
  else if(type_eta =="eta_lesser_1" && corr == "Roch"){
    std::cout<<" here in eta lesser 1 "<<std::endl;
     fit_low = {123.5, 123.5, 123.5, 123.5, 124.0};
     fit_up = {126.3, 126.3, 126.3 ,126.6, 127.1};
  }
  else if(type_eta =="eta_inclusive" && corr == "Roch"){
    std::cout<<" here in eta inclusive "<<std::endl;
     fit_low = {123.5, 123.5, 123.5, 123.6, 124.3};
     fit_up = {126.3, 126.3, 126.3 ,126.7, 127.8};

  }

    // Histograms to read 
   std::vector <TString> hist_list_Higgs = { "massHscaled_0per_Higgs", "massHscaled_0.03per_Higgs", "massHscaled_0.1per_Higgs", "massHscaled_0.3per_Higgs","massHscaled_1per_Higgs" };
  std::vector <TString> hist_list_Z = { "massHscaled_0per_Z", "massHscaled_0.03per_Z", "massHscaled_0.1per_Z", "massHscaled_0.3per_Z", "massHscaled_1per_Z"};
 
  std::vector<TH1D*> hist_massH; 
  std::vector<TH1D*> hist_massZ; 
    for(int i=0; i<hist_list_Higgs.size(); i++){
    TH1D *hH = (TH1D*) input_file->Get(hist_list_Higgs[i]);
    TH1D *hZ = (TH1D*) input_file->Get(hist_list_Z[i]);
    hist_massH.push_back(hH);
    hist_massZ.push_back(hZ);
    }

    TH1D* hist_gen_Higgs = (TH1D*) input_file->Get("GENmass4l_Higgs");
    TH1D* hist_gen_Z = (TH1D*) input_file->Get("GENmass4l_Z");
    TH1D* hist_gen_Higgs_noFSR = (TH1D*) input_file->Get("GENmass4l_Higgs_without_FSR");
    TH1D* hist_gen_Z_noFSR = (TH1D*) input_file->Get("GENmass4l_Z_without_FSR");

    for(int i=0; i<hist_list_Higgs.size(); i++){
    std::cout<<" hist name "<<hist_massH[i]->GetName()<<hist_massH[i]->GetEntries()<<std::endl;
    }


    // FIt the masses
    FitFunction fit_obj;
    FitFunction fit_obj2;
    fit_obj.initialise_fit("isHiggs", saving_path);
    fit_obj2.initialise_fit("isZ", saving_path);

    fit_obj.fitting_BW(hist_gen_Higgs, "gen_mass_Higgs");
    fit_obj2.fitting_BW(hist_gen_Z, "gen_mass_Z");
    fit_obj.fitting_BW(hist_gen_Higgs_noFSR, "gen_mass_Higgs_noFSR");
    fit_obj2.fitting_BW(hist_gen_Z_noFSR, "gen_mass_Z_noFSR");

    std::vector<float> mean_Z_value;
    std::vector<float> mean_error_Z_value;

    std::vector<float> mean_Higgs_value;
    std::vector<float> mean_error_Higgs_value;


    std::vector<float> ratio_vector;
    std::vector<float> ratio_err_vector;
    std::vector<float> mean_Z_vector, mean_Higgs_vector;
    std::vector<float> mean_Z_err_vector, mean_Higgs_err_vector;
    //int scale_string[6] = {0,1,3,5,10,50};
    TString scale_string[12] = {"0", "point_03", "point_1", "point_3","1" };
    std::vector<float> scale_factor = {0, 0.0003, 0.001, 0.003, 0.01};

    for(int i=0; i<hist_list_Higgs.size(); i++){
    std::pair<float, float> mean_pair_Higgs;
    std::pair<float, float> mean_pair_Z;

    TString save_name = TString::Format("scaled_massH_"+scale_string[i]+"per");
    mean_pair_Higgs = fit_obj.fitting_DSCB(hist_massH[i], save_name, scale_factor[i],fit_low[i], fit_up[i]);

    mean_Higgs_value.push_back(mean_pair_Higgs.first);
    mean_error_Higgs_value.push_back(mean_pair_Higgs.second);

    TString save_name_Z = TString::Format("scaled_massZ_"+scale_string[i]+"per");
//    mean_pair_Z = fit_obj2.fitting_BW_DSCB(hist_gen_Z, hist_massZ[i], save_name_Z, scale_factor[i]);
    mean_pair_Z = fit_obj2.fitting_DSCB(hist_massZ[i], save_name_Z, scale_factor[i], fit_low[i], fit_up[i]);

    mean_Z_value.push_back(mean_pair_Z.first);
    mean_error_Z_value.push_back(mean_pair_Z.second);

    double ratio_value = mean_pair_Higgs.first/mean_pair_Z.first;
    float A = mean_pair_Higgs.second/mean_pair_Higgs.first;
    float B = mean_pair_Z.second/mean_pair_Z.first;
    double ratio_error_value = ratio_value * sqrt( pow(A,2) + pow(B,2));

    std::cout<<" ratio value "<<ratio_value<<std::endl;
    std::cout<<" ratio error value "<<ratio_error_value<<std::endl;
    std::cout<<" mH value "<<mean_pair_Higgs.first<<std::endl;
    std::cout<<" mZ value "<<mean_pair_Z.first<<std::endl;
    mean_Higgs_vector.push_back(mean_pair_Higgs.first);
    mean_Higgs_err_vector.push_back(mean_pair_Higgs.second);

    mean_Z_vector.push_back(mean_pair_Z.first);
    mean_Z_err_vector.push_back(mean_pair_Z.second);

    ratio_vector.push_back(ratio_value);
    ratio_err_vector.push_back(ratio_error_value);
    }

    //std::vector<float> scale_factor = {0, 0.01 , 0.03, 0.05, 0.1, 0.5};
    TGraphErrors *gr = new TGraphErrors(5, scale_factor.data(), ratio_vector.data(), 0, ratio_err_vector.data());
    TGraphErrors *gr_Higgs = new TGraphErrors(5, scale_factor.data(), mean_Higgs_vector.data(), 0, mean_Higgs_err_vector.data());
    TGraphErrors *gr_Z = new TGraphErrors(5, scale_factor.data(), mean_Z_vector.data(), 0, mean_Z_err_vector.data());

    TCanvas *c1 = new TCanvas();
    c1->cd();
    gr_Higgs->SetTitle("mH with #mu in "+title+" smeared : "+corr);
    gr_Higgs->Draw("AP*");
    gr_Higgs->SetName("graph_Higgs");
    gr_Higgs->GetYaxis()->SetRangeUser(124.0, 127);
    gr_Higgs->GetXaxis()->SetTitle(" scale "); 
    gr_Higgs->GetYaxis()->SetTitle(" mH "); 
    gr_Higgs->SetLineColor(kRed);
    c1->SaveAs(saving_path+"mean_Higgs_value.pdf");

    TCanvas *c2 = new TCanvas();
    c2->cd();
    gr_Z->SetName("graph_Z");
    gr_Z->SetTitle("mZ with #mu in "+title+ " smeared : "+corr);
    gr_Z->GetYaxis()->SetRangeUser(90.5, 92.5);
    gr_Z->Draw("AP*");
    gr_Z->GetXaxis()->SetTitle(" scale "); 
    gr_Z->GetYaxis()->SetTitle(" mZ "); 
    gr_Z->SetLineColor(kRed);
    c2->SaveAs(saving_path+"mean_Z_value.pdf");

    TCanvas *c = new TCanvas();
    c->cd();
    gr->SetName("graph_ratio");
    gr->SetTitle("mH/mZ with #mu in "+title+" smeared :  "+corr);
    gr->GetYaxis()->SetRangeUser(1.35,1.38);
    gr->Draw("AP*");
    gr->GetXaxis()->SetTitle(" scale "); 
    gr->GetYaxis()->SetTitle(" mH/mZ "); 
    gr->SetLineColor(kRed);
    c->SaveAs(saving_path+"graph_ratio_value.pdf");

//    TFile *output_file =  TFile::Open(output_file_name,"RECREATE");
//    output_file->cd();
//    gr->Write();
//    gr_Higgs->Write();
//    gr_Z->Write();
//    output_file->Close();
//
//    delete output_file;

  return 0; 

}
