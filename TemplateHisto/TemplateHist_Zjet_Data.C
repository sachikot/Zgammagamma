#define TemplateHist_Zjet_Data_cxx
#include "TemplateHist_Zjet_Data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "iostream"
#include "sstream"
#include <TMath.h>

void TemplateHist_Zjet_Data::Loop(TString name, TString fitVal, int iso)
{
  TFile* file = TFile::Open(name, "RECREATE");
  const Int_t ptBINS = 3;
  Double_t edgesptBINS[ptBINS + 1] = {15,25,40,1000000};

  std::map<std::string, TH1D*> hM;
  for(int i=0;i<ptBINS;i++){
    std::stringstream ii;
    ii<<i+1;
    if(fitVal == "Sieie"){
      hM[(std::string("h_Sieie_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_Sieie_EB_")+ii.str()).c_str(),"h_Sieie_EB",30, 0, 0.03);
      hM[(std::string("h_Sieie_EE_")+ii.str()).c_str()]  = new TH1D((std::string("h_Sieie_EE_")+ii.str()).c_str(),"h_Sieie_EE",30, 0, 0.09);
    }
    else if(fitVal == "ChIso"){
      hM[(std::string("h_ChIso_EB_")+ii.str()).c_str()]  = new TH1D((std::string("h_ChIso_EB_")+ii.str()).c_str(),"h_ChIso_EB",30, 0, 45);
      hM[(std::string("h_ChIso_EE_")+ii.str()).c_str()]  = new TH1D((std::string("h_ChIso_EE_")+ii.str()).c_str(),"h_ChIso_EE",35, 0, 42);
    }
    else {
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

    int EB = 0; int EE = 0;
    if(fabs(pho1_eta)<1.479) EB = 1;
    else if(fabs(pho1_eta) > 1.57 && fabs(pho1_eta) < 2.5) EE = 1;
    else{
      ;
    }

    for(int i = 0; i < ptBINS; ++i){
      if(pho1_pt>edgesptBINS[i] && pho1_pt<edgesptBINS[i+1]){
	std::stringstream ii;
	ii<<i+1;
	if(EB==1){
	  if(fitVal == "Sieie"){
	    hM[(std::string("h_Sieie_EB_")+ii.str()).c_str()]->Fill(pho1_Sihih, totalWeight);
	  }
	  else if(fitVal == "ChIso"){
	    hM[(std::string("h_ChIso_EB_")+ii.str()).c_str()]->Fill(pho1_ChIso, totalWeight);
	  }
	  else {
	    ;
	  }
	} // EB == 1
	if(EE==1){
	  if(fitVal == "Sieie"){
	    hM[(std::string("h_Sieie_EE_")+ii.str()).c_str()]->Fill(pho1_Sihih, totalWeight);
	  }
	  else if(fitVal == "ChIso"){
	    hM[(std::string("h_ChIso_EE_")+ii.str()).c_str()]->Fill(pho1_ChIso, totalWeight);
	  }
	  else {
	    ;
	  }
	} // EE == 1
      }
    }

  } // loop over for events
  file->Write();
  file->Close();
}

  int main(int argc, char* argv[]){

    TString fitVal = argv[1];
    int iso = atoi(argv[2]);

    for (int i = 0; i < argc; i++ ) {

      TemplateHist_Zjet_Data t("/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_WORK_09010/Zgg/RootForTemp/ZjetData/muon_Data_Zjet_"+fitVal+"_"+_iso[iso]+".root");
      t.Loop("ROOT_BkgData/muon_Bkg_"+fitVal+"_"+_iso[iso]+".root", fitVal, iso);
    }

    return 0;
  }
