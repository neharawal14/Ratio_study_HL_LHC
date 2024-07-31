void superimpose_graphs(){

//  std::vector<TString> type = {"4l", "2l"};
  std::vector<TString> type = {"4l"};
  for(int j=0; j<type.size(); j++){
//  std::vector<TString> type_eta = {"eta_greater_1", "eta_lesser_1", "eta_inclusive"};
  std::vector<TString> type_eta = {"eta_greater_1"};
//  std::vector<TString> eta_title = {"|#eta| >=1", "|#eta| <=1","|#eta| <=2.4"};
  std::vector<TString> eta_title = {"|#eta| >=1"};
  for(int i=0 ; i<type_eta.size(); i++){
   TString file1_name =  "output_graphs_noRoch_"+type_eta[i]+"_"+type[j]+".root";
   TString file2_name =  "output_graphs_Roch_"+type_eta[i]+"_"+type[j]+".root";
   std::cout<<" first file "<<file1_name<<std::endl;
   std::cout<<" second file "<<file2_name<<std::endl;
   TFile *f1 = TFile::Open(file1_name, "READ");
   TFile *f2 = TFile::Open(file2_name, "READ");
   TGraphErrors *gr1 = (TGraphErrors*) f1->Get("graph_ratio");
   TGraphErrors *gr2 = (TGraphErrors*) f2->Get("graph_ratio");
   TGraphErrors *gr_higgs_noRoch = (TGraphErrors*) f1->Get("graph_Higgs");
   TGraphErrors *gr_higgs_Roch = (TGraphErrors*) f2->Get("graph_Higgs");
   TGraphErrors *gr_Z_noRoch = (TGraphErrors*) f1->Get("graph_Z");
   TGraphErrors *gr_Z_Roch = (TGraphErrors*) f2->Get("graph_Z");

   TMultiGraph *mg_Higgs = new TMultiGraph();
   mg_Higgs->Add(gr_higgs_noRoch);
   mg_Higgs->Add(gr_higgs_Roch);
   gr_higgs_noRoch->SetLineColor(kRed);
   gr_higgs_noRoch->SetMarkerColor(kRed);
   gr_higgs_Roch->SetMarkerColor(kBlue);
   gr_higgs_Roch->SetLineColor(kBlue);

   TMultiGraph *mg_Z = new TMultiGraph();
   mg_Z->Add(gr_Z_noRoch);
   mg_Z->Add(gr_Z_Roch);
   gr_Z_noRoch->SetLineColor(kRed);
   gr_Z_noRoch->SetMarkerColor(kRed);
   gr_Z_Roch->SetMarkerColor(kBlue);
   gr_Z_Roch->SetLineColor(kBlue);

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gr1);
   mg->Add(gr2);
   gr1->SetLineColor(kRed);
   gr1->SetMarkerColor(kRed);
   gr2->SetMarkerColor(kBlue);
   gr2->SetLineColor(kBlue);

   TCanvas *c = new TCanvas();
   c->cd();
   c->SetTopMargin(0.1);

   TString title = " ratio(mH/mZ("+ type[j] +")) with #mu in "+eta_title[i]+" smeared ";
   mg->SetTitle(title);
   TLegend *leg = new TLegend(0.7, 0.93, 0.88, 0.98);
   mg->Draw("AP*");
   leg->AddEntry(gr1, "without Roch corr");
   leg->AddEntry(gr2, "with Roch corr");
   leg->Draw("same");
   c->SaveAs("ratio_results/ratio_"+type_eta[i]+"_"+type[j]+".pdf");
  c->Close();

   TCanvas *c1 = new TCanvas();
   c1->cd();
   c1->SetTopMargin(0.1);

   TString title1 = " mH with #mu in "+eta_title[i]+" smeared ";
   mg_Higgs->SetTitle(title1);
   TLegend *leg_Higgs = new TLegend(0.7, 0.93, 0.88, 0.99);
   mg_Higgs->Draw("AP*");
   leg_Higgs->AddEntry(gr_higgs_noRoch, "without Roch corr");
   leg_Higgs->AddEntry(gr_higgs_Roch, "with Roch corr");
   leg_Higgs->Draw("same");
   c1->SaveAs("ratio_results/mass_Higgs_"+type_eta[i]+"_"+type[j]+".pdf");
  c1->Close();

  TCanvas *c2 = new TCanvas();
   c2->cd();
   c2->SetTopMargin(0.1);

   TString title2 = " mZ with #mu in "+eta_title[i]+" smeared ";
   mg_Z->SetTitle(title2);
   TLegend *leg_Z = new TLegend(0.7, 0.93, 0.88, 0.99);
   mg_Z->Draw("AP*");
   leg_Z->AddEntry(gr_Z_noRoch, "without Roch corr");
   leg_Z->AddEntry(gr_Z_Roch, "with Roch corr");
   leg_Higgs->Draw("same");
   c2->SaveAs("ratio_results/mass_Z_"+type_eta[i]+"_"+type[j]+".pdf");
  c2->Close();

  }
}

}
