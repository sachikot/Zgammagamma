#define TemplateHist_Zgg_cxx
#include "TemplateHist_Zgg.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "iostream"
#include "sstream"
#include <TMath.h>

void TemplateHist_Zgg::Loop(TString name, TString fitVal, int iso, TString sampleDir, TString sample)
{
  TFile* file = TFile::Open(name, "RECREATE");
  const Int_t ptBINS = 1;
  //  Double_t edgesptBINS[ptBINS + 1] = {15,25,40,1000000};
  Double_t edgesptBINS[ptBINS + 1] = {15,1000000};

  std::map<std::string, TH1D*> hM;
  for(int i=0;i<ptBINS;i++){
    std::stringstream ii;
    ii<<i+1;
    if(fitVal == "Sieie"){
      hM[(std::string("h_Sieie_LeadPass_EB_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_Sieie_LeadPass_EB_EB_")+ii.str()).c_str(),"h_Sieie_LeadPass_EB_EB", 30, 0, 0.03);
      hM[(std::string("h_Sieie_LeadPass_EB_EE_")+ii.str()).c_str()]  = new TH1D((std::string("h_Sieie_LeadPass_EB_EE_")+ii.str()).c_str(),"h_Sieie_LeadPass_EB_EE", 30, 0, 0.09);
      hM[(std::string("h_Sieie_LeadPass_EE_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_Sieie_LeadPass_EE_EB_")+ii.str()).c_str(),"h_Sieie_LeadPass_EE_EB", 30, 0, 0.03);

      hM[(std::string("h_Sieie_LeadFail_EB_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_Sieie_LeadFail_EB_EB_")+ii.str()).c_str(),"h_Sieie_LeadFail_EB_EB", 30, 0, 0.03);
      hM[(std::string("h_Sieie_LeadFail_EB_EE_")+ii.str()).c_str()]  = new TH1D((std::string("h_Sieie_LeadFail_EB_EE_")+ii.str()).c_str(),"h_Sieie_LeadFail_EB_EE", 30, 0, 0.09);
      hM[(std::string("h_Sieie_LeadFail_EE_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_Sieie_LeadFail_EE_EB_")+ii.str()).c_str(),"h_Sieie_LeadFail_EE_EB", 30, 0, 0.03);
    }
    else if(fitVal == "ChIso"){
      hM[(std::string("h_ChIso_LeadPass_EB_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_ChIso_LeadPass_EB_EB_")+ii.str()).c_str(),"h_ChIso_LeadPass_EB_EB", 30, 0, 45);
      hM[(std::string("h_ChIso_LeadPass_EB_EE_")+ii.str()).c_str()]  = new TH1D((std::string("h_ChIso_LeadPass_EB_EE_")+ii.str()).c_str(),"h_ChIso_LeadPass_EB_EE", 35, 0, 42);
      hM[(std::string("h_ChIso_LeadPass_EE_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_ChIso_LeadPass_EE_EB_")+ii.str()).c_str(),"h_ChIso_LeadPass_EE_EB", 30, 0, 45);

      hM[(std::string("h_ChIso_LeadFail_EB_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_ChIso_LeadFail_EB_EB_")+ii.str()).c_str(),"h_ChIso_LeadFail_EB_EB", 30, 0, 45);
      hM[(std::string("h_ChIso_LeadFail_EB_EE_")+ii.str()).c_str()]  = new TH1D((std::string("h_ChIso_LeadFail_EB_EE_")+ii.str()).c_str(),"h_ChIso_LeadFail_EB_EE", 35, 0, 42);
      hM[(std::string("h_ChIso_LeadFail_EE_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_ChIso_LeadFail_EE_EB_")+ii.str()).c_str(),"h_ChIso_LeadFail_EE_EB", 30, 0, 45);
    }
    else{
      ;
    }
  }

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    int EBEB=0; int EBEE=0; int EEEB=0;
    if(fabs(pho1_eta)<1.479 && fabs(pho2_eta)<1.479) EBEB=1;
    else if( fabs(pho1_eta)<1.479 && (fabs(pho2_eta)>1.57&&fabs(pho2_eta)<2.5) ) EBEE=1;
    else if( (fabs(pho1_eta)>1.57&&fabs(pho1_eta)<2.5) && fabs(pho2_eta)<1.479 ) EEEB=1;
    else {
      ;
    }

    int Sieie_LeadPass_EBEB = 0; int Sieie_LeadPass_EBEE = 0; int Sieie_LeadPass_EEEB = 0;
    if(EBEB && pho1_Sihih<0.011) Sieie_LeadPass_EBEB = 1;
    else if(EBEE && pho1_Sihih<0.011) Sieie_LeadPass_EBEE = 1;
    else if(EEEB && pho1_Sihih<0.033) Sieie_LeadPass_EEEB = 1;
    else {
      ;
    }
    int Sieie_LeadFail_EBEB = 0; int Sieie_LeadFail_EBEE = 0; int Sieie_LeadFail_EEEB = 0;
    if(EBEB && pho1_Sihih>0.011) Sieie_LeadFail_EBEB = 1;
    else if(EBEE && pho1_Sihih>0.011) Sieie_LeadFail_EBEE = 1;
    else if(EEEB && pho1_Sihih>0.033) Sieie_LeadFail_EEEB = 1;
    else {
      ;
    }

    int ChIso_LeadPass_EBEB = 0; int ChIso_LeadPass_EBEE = 0; int ChIso_LeadPass_EEEB = 0;
    if(EBEB && pho1_ChIso<1.5) ChIso_LeadPass_EBEB = 1;
    else if(EBEE && pho1_ChIso<1.5) ChIso_LeadPass_EBEE = 1;
    else if(EEEB && pho1_ChIso<1.2) ChIso_LeadPass_EEEB = 1;
    else {
      ;
    }
    int ChIso_LeadFail_EBEB = 0; int ChIso_LeadFail_EBEE = 0; int ChIso_LeadFail_EEEB = 0;
    if(EBEB && pho1_ChIso>1.5) ChIso_LeadFail_EBEB = 1;
    else if(EBEE && pho1_ChIso>1.5) ChIso_LeadFail_EBEE = 1;
    else if(EEEB && pho1_ChIso>1.2) ChIso_LeadFail_EEEB = 1;
    else {
      ;
    }

    for(int i = 0; i < ptBINS; ++i){
      if(pho1_pt>edgesptBINS[i] && pho1_pt<edgesptBINS[i+1]){
	std::stringstream ii;
	ii<<i+1;
	if(fitVal == "Sieie"){
	  if(Sieie_LeadPass_EBEB==1) hM[(std::string("h_Sieie_LeadPass_EB_EB_")+ii.str()).c_str()]->Fill(pho2_Sihih);
	  if(Sieie_LeadPass_EBEE==1) hM[(std::string("h_Sieie_LeadPass_EB_EE_")+ii.str()).c_str()]->Fill(pho2_Sihih);
	  if(Sieie_LeadPass_EEEB==1) hM[(std::string("h_Sieie_LeadPass_EE_EB_")+ii.str()).c_str()]->Fill(pho2_Sihih);

	  if(Sieie_LeadFail_EBEB==1) hM[(std::string("h_Sieie_LeadFail_EB_EB_")+ii.str()).c_str()]->Fill(pho2_Sihih);
	  if(Sieie_LeadFail_EBEE==1) hM[(std::string("h_Sieie_LeadFail_EB_EE_")+ii.str()).c_str()]->Fill(pho2_Sihih);
	  if(Sieie_LeadFail_EEEB==1) hM[(std::string("h_Sieie_LeadFail_EE_EB_")+ii.str()).c_str()]->Fill(pho2_Sihih);
	} else if(fitVal == "ChIso"){
	  if(ChIso_LeadPass_EBEB==1) hM[(std::string("h_ChIso_LeadPass_EB_EB_")+ii.str()).c_str()]->Fill(pho2_ChIso);
	  if(ChIso_LeadPass_EBEE==1) hM[(std::string("h_ChIso_LeadPass_EB_EE_")+ii.str()).c_str()]->Fill(pho2_ChIso);
	  if(ChIso_LeadPass_EEEB==1) hM[(std::string("h_ChIso_LeadPass_EE_EB_")+ii.str()).c_str()]->Fill(pho2_ChIso);

	  if(ChIso_LeadFail_EBEB==1) hM[(std::string("h_ChIso_LeadFail_EB_EB_")+ii.str()).c_str()]->Fill(pho2_ChIso);
	  if(ChIso_LeadFail_EBEE==1) hM[(std::string("h_ChIso_LeadFail_EB_EE_")+ii.str()).c_str()]->Fill(pho2_ChIso);
	  if(ChIso_LeadFail_EEEB==1) hM[(std::string("h_ChIso_LeadFail_EE_EB_")+ii.str()).c_str()]->Fill(pho2_ChIso);
	}
	else {
	  ;
	}
      }
    }
  } // loop over for events
  file->Write();
  file->Close();
}

  int main(int argc, char* argv[]){

    TString fitVal = argv[1];
    int iso = atoi(argv[2]);
    TString sampleDir = argv[3];
    TString sample = argv[4];

    for (int i = 0; i < argc; i++ ) {

      TemplateHist_Zgg t("/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_WORK_09010/Zgg/RootForTemp/Zgg"+sampleDir+"/"+sample+"_Zgg_"+fitVal+"_"+_iso[iso]+".root");
      t.Loop("ROOT_Zgg/"+sample+"_Zgg_"+fitVal+"_"+_iso[iso]+".root", fitVal, iso, sampleDir, sample);
    }

    return 0;
  }
