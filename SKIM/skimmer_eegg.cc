#define skimmer_eegg_cxx
#include "skimmer_eegg.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>
#include "photon2012.cc"
using namespace std;

float skimmer_eegg::geteIDEA(int i, float SCEta){
  float EffectiveArea = -1;

  if (fabs(SCEta) < 1.0 )                              EffectiveArea = 0.130;
  else if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.137;
  else if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.067;
  else if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.089;
  else if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.107;
  else if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.110;
  else if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.138;
  else {
    cout << "Did not get Effective Area for eta " << SCEta << endl;
  }

  return EffectiveArea;
}

bool skimmer_eegg::eID2012(int i){

  float EffectiveArea = geteIDEA(i, (*eleSCEta)[i]);

  if((*elePt)[i] < 10) return false;

  bool pass = true;
  if (fabs((*eleSCEta)[i]) < 1.479) {

    if (!( fabs((*eledEtaAtVtx)[i]) < 0.007 )) pass = false;
    if (!( fabs((*eledPhiAtVtx)[i]) < 0.15 )) pass = false;
    if (!( (*eleSigmaIEtaIEta)[i] < 0.01 )) pass = false;
    if (!( (*eleHoverE)[i] < 0.12 )) pass = false;
    if (!( fabs((*eleD0GV)[i]) < 0.02 )) pass = false;
    if (!( fabs((*eleDzGV)[i]) < 0.2 )) pass = false;
    if (!( fabs(1./(*eleSCEn)[i] - (*eleEoverP)[i]/(*eleSCEn)[i]) < 0.05 )) pass = false;
    if (!( (*eleConvVtxFit)[i] == 0 )) pass = false;
    if (!( (*eleMissHits)[i] <= 1 )) pass = false;
    if (!( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0.,(Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] < 0.15 )) pass = false;

  } else {

    if (!( fabs((*eledEtaAtVtx)[i]) < 0.009 )) pass = false;
    if (!( fabs((*eledPhiAtVtx)[i]) < 0.10 )) pass = false;
    if (!( (*eleSigmaIEtaIEta)[i] < 0.03 )) pass = false;
    if (!( (*eleHoverE)[i] < 0.10 )) pass = false;
    if (!( fabs((*eleD0GV)[i]) < 0.02 )) pass = false;
    if (!( fabs((*eleDzGV)[i]) < 0.2 )) pass = false;
    if (!( fabs(1./(*eleSCEn)[i] - (*eleEoverP)[i]/(*eleSCEn)[i]) < 0.05 )) pass = false;
    if (!( (*eleConvVtxFit)[i] == 0 )) pass = false;
    if (!( (*eleMissHits)[i] <= 1 )) pass = false;

    if ((*elePt)[i] < 20){
      if (!( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0.,(Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] < 0.10 )) pass = false;
    } else {
      if (!( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0.,(Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] < 0.15 )) pass = false;
    }
  }

  return pass;
}

void skimmer_eegg::Loop(TString outputName, float xsScale, int processCode)
{

  TH1F* hEvents = (TH1F*)gDirectory->Get("hEvents");
  float nEvents = hEvents->GetBinContent(1);

  TH1F* hPU     = (TH1F*)gDirectory->Get("hPU");
  TH1F* hPUTrue = (TH1F*)gDirectory->Get("hPUTrue");

  // Get reference histogram
  TFile* fpuData = TFile::Open("Pileup_Observed_69300.root", "READ");
  TFile* fpuData_m5 = TFile::Open("Pileup_observed_69300_m5.root", "READ");
  TFile* fpuData_p5 = TFile::Open("Pileup_observed_69300_p5.root", "READ");
  TH1F* hDataPU_ = (TH1F*)fpuData->Get("pileup");
  TH1F* hDataPU_m5 = (TH1F*)fpuData_m5->Get("pileup");
  TH1F* hDataPU_p5 = (TH1F*)fpuData_p5->Get("pileup");

  // PU Reweighting 
  double npu_probs[101];
  double s = 0;
  vector<double> result(101);
  for(int npu = 0; npu < 101; ++npu)
    npu_probs[npu] = hPU->GetBinContent(npu+1)/hPU->GetEntries();

  for(int npu = 0; npu < 101; ++npu) {
    double npu_estimated = hDataPU_->GetBinContent(hDataPU_->GetXaxis()->FindBin(npu));
    if ( npu_probs[npu] != 0 ) result[npu] = npu_estimated / npu_probs[npu];
    else                       result[npu] = 0;
    s += npu_estimated;
  }     
  for (int npu = 0; npu < 101; ++npu) {
    result[npu] /= s;
  }
  // + 5%
  double s_p5 = 0;
  vector<double> result_p5(101);

  for(int npu = 0; npu < 101; ++npu) {
    double npu_estimated = hDataPU_p5->GetBinContent(hDataPU_p5->GetXaxis()->FindBin(npu));
    if ( npu_probs[npu] != 0 ) result_p5[npu] = npu_estimated / npu_probs[npu];
    else                       result_p5[npu] = 0;
    s_p5 += npu_estimated;
  }     
  for (int npu = 0; npu < 101; ++npu) result_p5[npu] /= s_p5;
  // -5%
  double s_m5 = 0;
  vector<double> result_m5(101);

  for(int npu = 0; npu < 101; ++npu) {
    double npu_estimated = hDataPU_m5->GetBinContent(hDataPU_m5->GetXaxis()->FindBin(npu));
    if ( npu_probs[npu] != 0 ) result_m5[npu] = npu_estimated / npu_probs[npu];
    else                       result_m5[npu] = 0;
    s_m5 += npu_estimated;
  }     
  for (int npu = 0; npu < 101; ++npu) result_m5[npu] /= s_m5;

  TFile* file = TFile::Open(outputName, "RECREATE");
  TTree* MyNewTree = fChain->CloneTree(0);
  MyNewTree->SetBranchStatus("*",0);
  MyNewTree->SetBranchStatus("*LHEWeight*",1);
  MyNewTree->SetBranchStatus("run",1);
  MyNewTree->SetBranchStatus("event",1);
  MyNewTree->SetBranchStatus("lumis",1);
  MyNewTree->SetBranchStatus("nHLT",1);
  MyNewTree->SetBranchStatus("HLT",1);
  MyNewTree->SetBranchStatus("HLTIndex",1);
  MyNewTree->SetBranchStatus("bspotPos",1);
  MyNewTree->SetBranchStatus("nVtx",1);
  MyNewTree->SetBranchStatus("IsVtxGood",1);
  MyNewTree->SetBranchStatus("nPUInfo",1);
  MyNewTree->SetBranchStatus("nPU",1);
  MyNewTree->SetBranchStatus("puBX",1);
  MyNewTree->SetBranchStatus("puTrue",1);
  MyNewTree->SetBranchStatus("nMC",1);
  MyNewTree->SetBranchStatus("mc*",1);
  MyNewTree->SetBranchStatus("pfMET*",1);
  MyNewTree->SetBranchStatus("nEle",1);
  MyNewTree->SetBranchStatus("ele*",1);
  MyNewTree->SetBranchStatus("nMu",1);
  MyNewTree->SetBranchStatus("mu*",1);
  MyNewTree->SetBranchStatus("rho*",1);
  MyNewTree->SetBranchStatus("nPho",1);
  MyNewTree->SetBranchStatus("pho*",1);
  MyNewTree->SetBranchStatus("metFilters*",1);

  MyNewTree->Branch("processCode", &processCode, "processCode/I");
  MyNewTree->Branch("xsScale", &xsScale, "xsScale/F");

  float puWeight, totalWeight;
  float puWeight_p5, totalWeight_p5;
  float puWeight_m5, totalWeight_m5;
  float xsWeight = 19.7*xsScale/hEvents->GetBinContent(1);  
  MyNewTree->Branch("xsWeight", &xsWeight, "xsWeight/F");

  MyNewTree->Branch("puWeight", &puWeight, "puWeight/F");
  MyNewTree->Branch("puWeight_m5", &puWeight_m5, "puWeight_m5/F");
  MyNewTree->Branch("puWeight_p5", &puWeight_p5, "puWeight_p5/F");
  MyNewTree->Branch("totalWeight", &totalWeight, "totalWeight/F");
  MyNewTree->Branch("totalWeight_p5", &totalWeight_p5, "totalWeight_p5/F");
  MyNewTree->Branch("totalWeight_m5", &totalWeight_m5, "totalWeight_m5/F");

  TH1F* h_count = new TH1F("h_count", "h_count", 20, 0, 20);

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  int selectedEvents = 0;
  int nTotalEvents = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    ++nTotalEvents;
    
    if ( jentry%10000 == 0 && jentry != 0 ) 
      cout << "Processed " << nTotalEvents << " events, selected = " << selectedEvents 
	   << "\t" << (float)selectedEvents/nTotalEvents << " selection efficiency" << endl;

    h_count->Fill(0);

    // FOR ELECTRON CHANNEL
    //Trigger
    int trigger = 0;
    if(HLT[HLTIndex[9]] != -1 && HLT[HLTIndex[9]] != 0 ) trigger = 1;
    if(!trigger){
      continue;
    }

    h_count->Fill(1);

    vector <int> ielectrons;
    for (int iele = 0; iele < nEle; ++ iele){

      if ((*elePt)[iele] < 10) continue;

      if ( fabs((*eleSCEta)[iele]) > 2.5 ) continue;

      if( eID2012(iele) ) ielectrons.push_back(iele);
    }
    if (ielectrons.size() > 0 ) h_count->Fill(2);
    if ( ielectrons.size() != 2 ) continue;
    h_count->Fill(3);
    if ( (*elePt)[ielectrons[0]] < 20 || (*elePt)[ielectrons[1]] < 10 ) continue;

    h_count->Fill(4);

    ++selectedEvents;
 
    if ( processCode != 0 ) {
      puWeight = result[(*nPU)[1]];
      puWeight_p5 = result_p5[(*nPU)[1]];
      puWeight_m5 = result_m5[(*nPU)[1]];
      totalWeight = puWeight*xsWeight;
      totalWeight_p5 = puWeight_p5*xsWeight;
      totalWeight_m5 = puWeight_m5*xsWeight;
    } else {
      puWeight = 1.0;
      totalWeight = 1.0;
      puWeight_p5 = 1.0;
      totalWeight_p5 = 1.0;
      puWeight_m5 = 1.0;
      totalWeight_m5 = 1.0;
    }
   
    MyNewTree->Fill();
  }
  h_count->Write();
  MyNewTree->Write();
  hPU->Write();
  hPUTrue->Write();
  hEvents->Write();
  file->Close();
}

int main(int argc, char* argv[]) {

  TString path = argv[1];
  TString listOfFiles = argv[2];

  std::vector<TString> fileNames;
  std::vector<double> nevents;
  std::vector<double> xs;
  std::vector<int> processes;

  TString fileName;
  double xSection;
  long nEvents;
  int processCode;
  // get files to skimmer
  ifstream InputStream(listOfFiles);
  while(!InputStream.eof()) {
    if ( !(InputStream >> fileName) ) break;
    if ( fileName[0] == '#' ) continue;
    if ( !(InputStream >> nEvents) ) break;
    if ( !(InputStream >> xSection) ) break;
    if ( !(InputStream >> processCode) ) break;

    fileNames.push_back(fileName);
    if ( nEvents != 1.0 ) {
      xs.push_back(xSection);
      nevents.push_back(nEvents);
      processes.push_back(processCode);
    } else {
      xs.push_back(1);
      nevents.push_back(1);
      processes.push_back(0);
    }
  }

  for(unsigned int i = 0; i != fileNames.size(); ++i) {
    
    cout << fileNames[i] << '\t' << xs[i] << '\t' << nevents[i] << '\t' << processes[i] << endl;
    
    TString outputName = "";
    for(int ichar = 0; ichar != fileNames[i].Length() - 5; ++ichar) 
      outputName += fileNames[i][ichar];
    outputName += "_skimmed_eegg.root";
    skimmer_eegg t(path + "/" + fileNames[i]);
    t.Loop(outputName, xs[i], processes[i]);   
  }

  return 0;
}

