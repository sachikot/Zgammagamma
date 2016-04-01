#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>

using namespace std;

void Template_SubtructedBkg(){
  gROOT->SetBatch();
  gROOT->Reset();

  TString _iso[7] = {"nomi", "533", "855", "1077", "1299", "151111", "201616"};
  TString fitVal = "Sieie";
  TString path = "/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_WORK_09010/Zgg/Muon/TempHist/TEMP_ROOT_0111/";

  //  for(int i = 0; i < 7; i++){

  TString iso = _iso[6];

  TFile* f1 = new TFile(path + "ROOT_BkgData/muon_Bkg_"+fitVal+"_"+iso+".root", "READ");
  TFile* f2 = new TFile(path + "ROOT_BkgMC/muon_Bkg_"+fitVal+"_"+iso+".root", "READ");
  TFile* file = new TFile("ROOT_Bkg/muon_subtrBkg_"+fitVal+"_"+iso+".root", "recreate");

  gStyle->SetOptStat(0);

  TH1F* h_data_bin1_eb = (TH1F*)f1->Get("h_"+fitVal+"_EB_1");
  TH1F* h_data_bin2_eb = (TH1F*)f1->Get("h_"+fitVal+"_EB_2");
  TH1F* h_data_bin3_eb = (TH1F*)f1->Get("h_"+fitVal+"_EB_3");
  TH1F* h_data_bin1_ee = (TH1F*)f1->Get("h_"+fitVal+"_EE_1");
  TH1F* h_data_bin2_ee = (TH1F*)f1->Get("h_"+fitVal+"_EE_2");
  TH1F* h_data_bin3_ee = (TH1F*)f1->Get("h_"+fitVal+"_EE_3");

  TH1F* h_mc_bin1_eb = (TH1F*)f2->Get("h_"+fitVal+"_EB_1");
  TH1F* h_mc_bin2_eb = (TH1F*)f2->Get("h_"+fitVal+"_EB_2");
  TH1F* h_mc_bin3_eb = (TH1F*)f2->Get("h_"+fitVal+"_EB_3");
  TH1F* h_mc_bin1_ee = (TH1F*)f2->Get("h_"+fitVal+"_EE_1");
  TH1F* h_mc_bin2_ee = (TH1F*)f2->Get("h_"+fitVal+"_EE_2");
  TH1F* h_mc_bin3_ee = (TH1F*)f2->Get("h_"+fitVal+"_EE_3");

  h_data_bin1_eb->Add(h_mc_bin1_eb, -1);
  h_data_bin2_eb->Add(h_mc_bin2_eb, -1);
  h_data_bin3_eb->Add(h_mc_bin3_eb, -1);
  h_data_bin1_ee->Add(h_mc_bin1_ee, -1);
  h_data_bin2_ee->Add(h_mc_bin2_ee, -1);
  h_data_bin3_ee->Add(h_mc_bin3_ee, -1);

  h_data_bin1_eb->Write();
  h_data_bin2_eb->Write();
  h_data_bin3_eb->Write();
  h_data_bin1_ee->Write();
  h_data_bin2_ee->Write();
  h_data_bin3_ee->Write();

  file->Write();
  file->Close();
  //  }
}
