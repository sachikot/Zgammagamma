#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TTree.h>
#include <iostream>

void TemplateHisto_FF(){

  //  gStyle->SetOptStat(0);
  gROOT->SetBatch();
  gROOT->Reset();

  TString sieie_eb[2]  = {"0.01100", "0.02900"};  
  TString sieie_ee[2]  = {"0.03300", "0.08700"};
  TString phoiso_chiso_eb[5] = {"0.700", "5.000", "3.000", "7.000", "9.000"};  
  TString phoiso_chiso_ee[5] = {"1.000", "5.000", "3.000", "7.000", "9.000"};
  TString phoiso_sieie_eb[5] = {"0.700", "9.000", "7.000", "11.000", "16.000"};  
  TString phoiso_sieie_ee[5] = {"1.000", "9.000", "7.000", "11.000", "16.000"};
  TString chiso_eb[5] = {"1.5", "12.0", "10.0", "15.0", "20.0"};
  TString chiso_ee[5] = {"1.2", "12.0", "10.0", "15.0", "20.0"};
  TString sideband[5] = {"", "nomi", "Tight", "Loose", "Loosest"};
  int iso = 1;

  TString path = "/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_WORK_09010/Zgg/file_Josh/file_0106/FakeFakeTemp/";
  TFile* file   = new TFile(path + "job_muon_2012_Jan22rereco.root", "read");
  TTree* tree   = (TTree*)file->Get("ggNtuplizer/EventTree");
  //  TFile* output = new TFile("ROOT_CorrFF/Temp_Corr_FF_"+sideband[iso]+".root", "recreate");
  TFile* output = new TFile("ROOT_CorrFF/Temp_Corr_FF_"+sideband[iso]+".root", "recreate");

  // Histogram to Sieie
  TH1D* h_Sieie_LeadPass_EB_EB_1 = new TH1D("h_Sieie_LeadPass_EB_EB_1", "", 30, 0, 0.03);
  TH1D* h_Sieie_LeadPass_EB_EE_1 = new TH1D("h_Sieie_LeadPass_EB_EE_1", "", 30, 0, 0.09);
  TH1D* h_Sieie_LeadPass_EE_EB_1 = new TH1D("h_Sieie_LeadPass_EE_EB_1", "", 30, 0, 0.03);

  TH1D* h_Sieie_LeadFail_EB_EB_1 = new TH1D("h_Sieie_LeadFail_EB_EB_1", "", 30, 0, 0.03);
  TH1D* h_Sieie_LeadFail_EB_EE_1 = new TH1D("h_Sieie_LeadFail_EB_EE_1", "", 30, 0, 0.09);
  TH1D* h_Sieie_LeadFail_EE_EB_1 = new TH1D("h_Sieie_LeadFail_EE_EB_1", "", 30, 0, 0.03);

  TH1D* h_Sieie_LeadPass_EB_EB_2 = new TH1D("h_Sieie_LeadPass_EB_EB_2", "", 30, 0, 0.03);
  TH1D* h_Sieie_LeadPass_EB_EE_2 = new TH1D("h_Sieie_LeadPass_EB_EE_2", "", 30, 0, 0.09);
  TH1D* h_Sieie_LeadPass_EE_EB_2 = new TH1D("h_Sieie_LeadPass_EE_EB_2", "", 30, 0, 0.03);

  TH1D* h_Sieie_LeadFail_EB_EB_2 = new TH1D("h_Sieie_LeadFail_EB_EB_2", "", 30, 0, 0.03);
  TH1D* h_Sieie_LeadFail_EB_EE_2 = new TH1D("h_Sieie_LeadFail_EB_EE_2", "", 30, 0, 0.09);
  TH1D* h_Sieie_LeadFail_EE_EB_2 = new TH1D("h_Sieie_LeadFail_EE_EB_2", "", 30, 0, 0.03);

  TH1D* h_Sieie_LeadPass_EB_EB_3 = new TH1D("h_Sieie_LeadPass_EB_EB_3", "", 30, 0, 0.03);
  TH1D* h_Sieie_LeadPass_EB_EE_3 = new TH1D("h_Sieie_LeadPass_EB_EE_3", "", 30, 0, 0.09);
  TH1D* h_Sieie_LeadPass_EE_EB_3 = new TH1D("h_Sieie_LeadPass_EE_EB_3", "", 30, 0, 0.03);

  TH1D* h_Sieie_LeadFail_EB_EB_3 = new TH1D("h_Sieie_LeadFail_EB_EB_3", "", 30, 0, 0.03);
  TH1D* h_Sieie_LeadFail_EB_EE_3 = new TH1D("h_Sieie_LeadFail_EB_EE_3", "", 30, 0, 0.09);
  TH1D* h_Sieie_LeadFail_EE_EB_3 = new TH1D("h_Sieie_LeadFail_EE_EB_3", "", 30, 0, 0.03);

  // Histogram to ChIso
  TH1D* h_ChIso_LeadPass_EB_EB_1 = new TH1D("h_ChIso_LeadPass_EB_EB_1", "", 30, 0, 45);
  TH1D* h_ChIso_LeadPass_EB_EE_1 = new TH1D("h_ChIso_LeadPass_EB_EE_1", "", 35, 0, 42);
  TH1D* h_ChIso_LeadPass_EE_EB_1 = new TH1D("h_ChIso_LeadPass_EE_EB_1", "", 30, 0, 45);

  TH1D* h_ChIso_LeadFail_EB_EB_1 = new TH1D("h_ChIso_LeadFail_EB_EB_1", "", 30, 0, 45);
  TH1D* h_ChIso_LeadFail_EB_EE_1 = new TH1D("h_ChIso_LeadFail_EB_EE_1", "", 35, 0, 42);
  TH1D* h_ChIso_LeadFail_EE_EB_1 = new TH1D("h_ChIso_LeadFail_EE_EB_1", "", 30, 0, 45);

  TH1D* h_ChIso_LeadPass_EB_EB_2 = new TH1D("h_ChIso_LeadPass_EB_EB_2", "", 30, 0, 45);
  TH1D* h_ChIso_LeadPass_EB_EE_2 = new TH1D("h_ChIso_LeadPass_EB_EE_2", "", 35, 0, 42);
  TH1D* h_ChIso_LeadPass_EE_EB_2 = new TH1D("h_ChIso_LeadPass_EE_EB_2", "", 30, 0, 45);

  TH1D* h_ChIso_LeadFail_EB_EB_2 = new TH1D("h_ChIso_LeadFail_EB_EB_2", "", 30, 0, 45);
  TH1D* h_ChIso_LeadFail_EB_EE_2 = new TH1D("h_ChIso_LeadFail_EB_EE_2", "", 35, 0, 42);
  TH1D* h_ChIso_LeadFail_EE_EB_2 = new TH1D("h_ChIso_LeadFail_EE_EB_2", "", 30, 0, 45);

  TH1D* h_ChIso_LeadPass_EB_EB_3 = new TH1D("h_ChIso_LeadPass_EB_EB_3", "", 30, 0, 45);
  TH1D* h_ChIso_LeadPass_EB_EE_3 = new TH1D("h_ChIso_LeadPass_EB_EE_3", "", 35, 0, 42);
  TH1D* h_ChIso_LeadPass_EE_EB_3 = new TH1D("h_ChIso_LeadPass_EE_EB_3", "", 30, 0, 45);

  TH1D* h_ChIso_LeadFail_EB_EB_3 = new TH1D("h_ChIso_LeadFail_EB_EB_3", "", 30, 0, 45);
  TH1D* h_ChIso_LeadFail_EB_EE_3 = new TH1D("h_ChIso_LeadFail_EB_EE_3", "", 35, 0, 42);
  TH1D* h_ChIso_LeadFail_EE_EB_3 = new TH1D("h_ChIso_LeadFail_EE_EB_3", "", 30, 0, 45);

  // Sieie Sel
  TString ff_sieie_leadpass_ebeb = "ph_noSIEIEiso1299_n==2 && chIsoCorr_leadph12>"+chiso_eb[0]+"&& chIsoCorr_leadph12<"+chiso_eb[iso]+"&&chIsoCorr_sublph12>"+chiso_eb[0]+"&&chIsoCorr_sublph12<"+chiso_eb[iso]+"&&ph_passNeuIsoCorrMedium[0]&&ph_passNeuIsoCorrMedium[1]&&phoIsoCorr_leadph12>"+phoiso_sieie_eb[0]+"&&phoIsoCorr_leadph12<"+phoiso_sieie_eb[iso]+"&&phoIsoCorr_sublph12>"+phoiso_sieie_eb[0]+"&&phoIsoCorr_sublph12<"+phoiso_sieie_eb[iso]+"&&isEB_leadph12&&isEB_sublph12&&mu_n==1 && el_n==0&& dr_ph1_ph2>0.4&&dr_ph1_leadLep>0.4&&dr_ph2_leadLep>0.4 && sieie_leadph12<0.011";

  TString ff_sieie_leadpass_ebee = "ph_noSIEIEiso1299_n==2 && chIsoCorr_leadph12>"+chiso_eb[0]+"&& chIsoCorr_leadph12<"+chiso_eb[iso]+"&&chIsoCorr_sublph12>"+chiso_ee[0]+"&&chIsoCorr_sublph12<"+chiso_ee[iso]+"&&ph_passNeuIsoCorrMedium[0]&&ph_passNeuIsoCorrMedium[1]&&phoIsoCorr_leadph12>"+phoiso_sieie_eb[0]+"&&phoIsoCorr_leadph12<"+phoiso_sieie_eb[iso]+"&&phoIsoCorr_sublph12>"+phoiso_sieie_ee[0]+"&&phoIsoCorr_sublph12<"+phoiso_sieie_ee[iso]+"&&isEB_leadph12&&isEE_sublph12&&mu_n==1 && el_n==0&& dr_ph1_ph2>0.4&&dr_ph1_leadLep>0.4&&dr_ph2_leadLep>0.4&&sieie_leadph12<0.011";

  TString ff_sieie_leadpass_eeeb = "ph_noSIEIEiso1299_n==2 && chIsoCorr_leadph12>"+chiso_ee[0]+"&& chIsoCorr_leadph12<"+chiso_ee[iso]+"&&chIsoCorr_sublph12>"+chiso_eb[0]+"&&chIsoCorr_sublph12<"+chiso_eb[iso]+"&&ph_passNeuIsoCorrMedium[0]&&ph_passNeuIsoCorrMedium[1]&&phoIsoCorr_leadph12>"+phoiso_sieie_ee[0]+"&&phoIsoCorr_leadph12<"+phoiso_sieie_ee[iso]+"&&phoIsoCorr_sublph12>"+phoiso_sieie_eb[0]+"&&phoIsoCorr_sublph12<"+phoiso_sieie_eb[iso]+"&&isEE_leadph12&&isEB_sublph12&&mu_n==1 && el_n==0&& dr_ph1_ph2>0.4&&dr_ph1_leadLep>0.4&&dr_ph2_leadLep>0.4&&sieie_leadph12<0.033";

  TString ff_sieie_leadfail_ebeb = "ph_noSIEIEiso1299_n==2 && chIsoCorr_leadph12>"+chiso_eb[0]+"&& chIsoCorr_leadph12<"+chiso_eb[iso]+"&&chIsoCorr_sublph12>"+chiso_eb[0]+"&&chIsoCorr_sublph12<"+chiso_eb[iso]+"&&ph_passNeuIsoCorrMedium[0]&&ph_passNeuIsoCorrMedium[1]&&phoIsoCorr_leadph12>"+phoiso_sieie_eb[0]+"&&phoIsoCorr_leadph12<"+phoiso_sieie_eb[iso]+"&&phoIsoCorr_sublph12>"+phoiso_sieie_eb[0]+"&&phoIsoCorr_sublph12<"+phoiso_sieie_eb[iso]+"&&isEB_leadph12&&isEB_sublph12&&mu_n==1 && el_n==0&& dr_ph1_ph2>0.4&&dr_ph1_leadLep>0.4&&dr_ph2_leadLep>0.4 && sieie_leadph12>0.011";

  TString ff_sieie_leadfail_ebee = "ph_noSIEIEiso1299_n==2 && chIsoCorr_leadph12>"+chiso_eb[0]+"&& chIsoCorr_leadph12<"+chiso_eb[iso]+"&&chIsoCorr_sublph12>"+chiso_ee[0]+"&&chIsoCorr_sublph12<"+chiso_ee[iso]+"&&ph_passNeuIsoCorrMedium[0]&&ph_passNeuIsoCorrMedium[1]&&phoIsoCorr_leadph12>"+phoiso_sieie_eb[0]+"&&phoIsoCorr_leadph12<"+phoiso_sieie_eb[iso]+"&&phoIsoCorr_sublph12>"+phoiso_sieie_ee[0]+"&&phoIsoCorr_sublph12<"+phoiso_sieie_ee[iso]+"&&isEB_leadph12&&isEE_sublph12&&mu_n==1 && el_n==0&& dr_ph1_ph2>0.4&&dr_ph1_leadLep>0.4&&dr_ph2_leadLep>0.4&&sieie_leadph12>0.011";

  TString ff_sieie_leadfail_eeeb = "ph_noSIEIEiso1299_n==2 && chIsoCorr_leadph12>"+chiso_ee[0]+"&& chIsoCorr_leadph12<"+chiso_ee[iso]+"&&chIsoCorr_sublph12>"+chiso_eb[0]+"&&chIsoCorr_sublph12<"+chiso_eb[iso]+"&&ph_passNeuIsoCorrMedium[0]&&ph_passNeuIsoCorrMedium[1]&&phoIsoCorr_leadph12>"+phoiso_sieie_ee[0]+"&&phoIsoCorr_leadph12<"+phoiso_sieie_ee[iso]+"&&phoIsoCorr_sublph12>"+phoiso_sieie_eb[0]+"&&phoIsoCorr_sublph12<"+phoiso_sieie_eb[iso]+"&&isEE_leadph12&&isEB_sublph12&&mu_n==1 && el_n==0&& dr_ph1_ph2>0.4&&dr_ph1_leadLep>0.4&&dr_ph2_leadLep>0.4&&sieie_leadph12>0.033";

  ////
  TString ff_chiso_leadpass_ebeb = "mu_n==1 && el_n==0 && dr_ph1_ph2 > 0.4 && dr_ph1_leadLep>0.4 && dr_ph2_leadLep>0.4 && ph_failSIEIEisoNone55_n==2 && sieie_leadph12 >"+sieie_eb[0]+"&& sieie_sublph12 >"+sieie_eb[0]+"&& ph_passNeuIsoCorrMedium[0] && ph_passNeuIsoCorrMedium[1] && phoIsoCorr_leadph12 >"+phoiso_chiso_eb[0]+"&& phoIsoCorr_sublph12 >"+phoiso_chiso_eb[0]+"&& phoIsoCorr_leadph12 <"+phoiso_chiso_eb[iso]+"&& phoIsoCorr_sublph12 <"+phoiso_chiso_eb[iso]+"&& isEB_leadph12 && isEB_sublph12 && chIsoCorr_leadph12<1.5";

  TString ff_chiso_leadpass_ebee = "mu_n==1 && el_n==0 && dr_ph1_ph2 > 0.4 && dr_ph1_leadLep>0.4 && dr_ph2_leadLep>0.4 && ph_failSIEIEisoNone55_n==2 && sieie_leadph12 >"+sieie_eb[0]+"&& sieie_sublph12 >"+sieie_ee[0]+"&& ph_passNeuIsoCorrMedium[0] && ph_passNeuIsoCorrMedium[1] && phoIsoCorr_leadph12 >"+phoiso_chiso_eb[0]+"&& phoIsoCorr_sublph12 >"+phoiso_chiso_ee[0]+"&& phoIsoCorr_leadph12 <"+phoiso_chiso_eb[iso]+"&& phoIsoCorr_sublph12 <"+phoiso_chiso_ee[iso]+"&& isEB_leadph12 && isEE_sublph12 && chIsoCorr_leadph12<1.5";

  TString ff_chiso_leadpass_eeeb = "mu_n==1 && el_n==0 && dr_ph1_ph2 > 0.4 && dr_ph1_leadLep>0.4 && dr_ph2_leadLep>0.4 && ph_failSIEIEisoNone55_n==2 && sieie_leadph12 >"+sieie_ee[0]+"&& sieie_sublph12 >"+sieie_eb[0]+"&& ph_passNeuIsoCorrMedium[0] && ph_passNeuIsoCorrMedium[1] && phoIsoCorr_leadph12 >"+phoiso_chiso_ee[0]+"&& phoIsoCorr_sublph12 >"+phoiso_chiso_eb[0]+"&& phoIsoCorr_leadph12 <"+phoiso_chiso_ee[iso]+"&& phoIsoCorr_sublph12 <"+phoiso_chiso_eb[iso]+"&& isEE_leadph12 && isEB_sublph12 && chIsoCorr_leadph12<1.2";

  TString ff_chiso_leadfail_ebeb = "mu_n==1 && el_n==0 && dr_ph1_ph2 > 0.4 && dr_ph1_leadLep>0.4 && dr_ph2_leadLep>0.4 && ph_failSIEIEisoNone55_n==2 && sieie_leadph12 >"+sieie_eb[0]+"&& sieie_sublph12 >"+sieie_eb[0]+"&& ph_passNeuIsoCorrMedium[0] && ph_passNeuIsoCorrMedium[1] && phoIsoCorr_leadph12 >"+phoiso_chiso_eb[0]+"&& phoIsoCorr_sublph12 >"+phoiso_chiso_eb[0]+"&& phoIsoCorr_leadph12 <"+phoiso_chiso_eb[iso]+"&& phoIsoCorr_sublph12 <"+phoiso_chiso_eb[iso]+"&& isEB_leadph12 && isEB_sublph12 && chIsoCorr_leadph12>1.5";

  TString ff_chiso_leadfail_ebee = "mu_n==1 && el_n==0 && dr_ph1_ph2 > 0.4 && dr_ph1_leadLep>0.4 && dr_ph2_leadLep>0.4 && ph_failSIEIEisoNone55_n==2 && sieie_leadph12 >"+sieie_eb[0]+"&& sieie_sublph12 >"+sieie_ee[0]+"&& ph_passNeuIsoCorrMedium[0] && ph_passNeuIsoCorrMedium[1] && phoIsoCorr_leadph12 >"+phoiso_chiso_eb[0]+"&& phoIsoCorr_sublph12 >"+phoiso_chiso_ee[0]+"&& phoIsoCorr_leadph12 <"+phoiso_chiso_eb[iso]+"&& phoIsoCorr_sublph12 <"+phoiso_chiso_ee[iso]+"&& isEB_leadph12 && isEE_sublph12 && chIsoCorr_leadph12>1.5";

  TString ff_chiso_leadfail_eeeb = "mu_n==1 && el_n==0 && dr_ph1_ph2 > 0.4 && dr_ph1_leadLep>0.4 && dr_ph2_leadLep>0.4 && ph_failSIEIEisoNone55_n==2 && sieie_leadph12 >"+sieie_ee[0]+"&& sieie_sublph12 >"+sieie_eb[0]+"&& ph_passNeuIsoCorrMedium[0] && ph_passNeuIsoCorrMedium[1] && phoIsoCorr_leadph12 >"+phoiso_chiso_ee[0]+"&& phoIsoCorr_sublph12 >"+phoiso_chiso_eb[0]+"&& phoIsoCorr_leadph12 <"+phoiso_chiso_ee[iso]+"&& phoIsoCorr_sublph12 <"+phoiso_chiso_eb[iso]+"&& isEE_leadph12 && isEB_sublph12 && chIsoCorr_leadph12>1.2";

  TString bin1 = "&&pt_leadph12>15.0 && pt_leadph12 < 25.0";
  TString bin2 = "&&pt_leadph12>25.0 && pt_leadph12 < 40.0";
  TString bin3 = "&&pt_leadph12>40.0";

  tree->Draw("sieie_sublph12>>h_Sieie_LeadPass_EB_EB_1", ff_sieie_leadpass_ebeb + bin1);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadPass_EB_EE_1", ff_sieie_leadpass_ebee + bin1);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadPass_EE_EB_1", ff_sieie_leadpass_eeeb + bin1);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadFail_EB_EB_1", ff_sieie_leadfail_ebeb + bin1);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadFail_EB_EE_1", ff_sieie_leadfail_ebee + bin1);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadFail_EE_EB_1", ff_sieie_leadfail_eeeb + bin1);

  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadPass_EB_EB_1", ff_chiso_leadpass_ebeb + bin1);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadPass_EB_EE_1", ff_chiso_leadpass_ebee + bin1);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadPass_EE_EB_1", ff_chiso_leadpass_eeeb + bin1);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadFail_EB_EB_1", ff_chiso_leadfail_ebeb + bin1);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadFail_EB_EE_1", ff_chiso_leadfail_ebee + bin1);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadFail_EE_EB_1", ff_chiso_leadfail_eeeb + bin1);

  tree->Draw("sieie_sublph12>>h_Sieie_LeadPass_EB_EB_2", ff_sieie_leadpass_ebeb + bin2);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadPass_EB_EE_2", ff_sieie_leadpass_ebee + bin2);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadPass_EE_EB_2", ff_sieie_leadpass_eeeb + bin2);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadFail_EB_EB_2", ff_sieie_leadfail_ebeb + bin2);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadFail_EB_EE_2", ff_sieie_leadfail_ebee + bin2);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadFail_EE_EB_2", ff_sieie_leadfail_eeeb + bin2);

  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadPass_EB_EB_2", ff_chiso_leadpass_ebeb + bin2);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadPass_EB_EE_2", ff_chiso_leadpass_ebee + bin2);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadPass_EE_EB_2", ff_chiso_leadpass_eeeb + bin2);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadFail_EB_EB_2", ff_chiso_leadfail_ebeb + bin2);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadFail_EB_EE_2", ff_chiso_leadfail_ebee + bin2);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadFail_EE_EB_2", ff_chiso_leadfail_eeeb + bin2);

  tree->Draw("sieie_sublph12>>h_Sieie_LeadPass_EB_EB_3", ff_sieie_leadpass_ebeb + bin3);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadPass_EB_EE_3", ff_sieie_leadpass_ebee + bin3);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadPass_EE_EB_3", ff_sieie_leadpass_eeeb + bin3);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadFail_EB_EB_3", ff_sieie_leadfail_ebeb + bin3);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadFail_EB_EE_3", ff_sieie_leadfail_ebee + bin3);
  tree->Draw("sieie_sublph12>>h_Sieie_LeadFail_EE_EB_3", ff_sieie_leadfail_eeeb + bin3);

  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadPass_EB_EB_3", ff_chiso_leadpass_ebeb + bin3);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadPass_EB_EE_3", ff_chiso_leadpass_ebee + bin3);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadPass_EE_EB_3", ff_chiso_leadpass_eeeb + bin3);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadFail_EB_EB_3", ff_chiso_leadfail_ebeb + bin3);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadFail_EB_EE_3", ff_chiso_leadfail_ebee + bin3);
  tree->Draw("chIsoCorr_sublph12>>h_ChIso_LeadFail_EE_EB_3", ff_chiso_leadfail_eeeb + bin3);

  h_Sieie_LeadPass_EB_EB_1->Write();
  h_Sieie_LeadPass_EB_EE_1->Write();
  h_Sieie_LeadPass_EE_EB_1->Write();
  h_Sieie_LeadFail_EB_EB_1->Write();
  h_Sieie_LeadFail_EB_EE_1->Write();
  h_Sieie_LeadFail_EE_EB_1->Write();

  h_ChIso_LeadPass_EB_EB_1->Write();
  h_ChIso_LeadPass_EB_EE_1->Write();
  h_ChIso_LeadPass_EE_EB_1->Write();
  h_ChIso_LeadFail_EB_EB_1->Write();
  h_ChIso_LeadFail_EB_EE_1->Write();
  h_ChIso_LeadFail_EE_EB_1->Write();

  h_Sieie_LeadPass_EB_EB_2->Write();
  h_Sieie_LeadPass_EB_EE_2->Write();
  h_Sieie_LeadPass_EE_EB_2->Write();
  h_Sieie_LeadFail_EB_EB_2->Write();
  h_Sieie_LeadFail_EB_EE_2->Write();
  h_Sieie_LeadFail_EE_EB_2->Write();

  h_ChIso_LeadPass_EB_EB_2->Write();
  h_ChIso_LeadPass_EB_EE_2->Write();
  h_ChIso_LeadPass_EE_EB_2->Write();
  h_ChIso_LeadFail_EB_EB_2->Write();
  h_ChIso_LeadFail_EB_EE_2->Write();
  h_ChIso_LeadFail_EE_EB_2->Write();

  h_Sieie_LeadPass_EB_EB_3->Write();
  h_Sieie_LeadPass_EB_EE_3->Write();
  h_Sieie_LeadPass_EE_EB_3->Write();
  h_Sieie_LeadFail_EB_EB_3->Write();
  h_Sieie_LeadFail_EB_EE_3->Write();
  h_Sieie_LeadFail_EE_EB_3->Write();

  h_ChIso_LeadPass_EB_EB_3->Write();
  h_ChIso_LeadPass_EB_EE_3->Write();
  h_ChIso_LeadPass_EE_EB_3->Write();
  h_ChIso_LeadFail_EB_EB_3->Write();
  h_ChIso_LeadFail_EB_EE_3->Write();
  h_ChIso_LeadFail_EE_EB_3->Write();

  output->Close();
}
