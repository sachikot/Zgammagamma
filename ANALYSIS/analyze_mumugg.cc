#define analyze_Zgg_cxx
#include "analyze_Zgg.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <TMath.h>
#include <TGraph.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

int nVtx_;
int run_;
long event_;

// electron
float ele1_pt;
float ele1_phi;
float ele1_eta;
float ele2_pt;
float ele2_phi;
float ele2_eta;
float mass_ee;
float mass_eeg;
float mass_eegg;

// photon
float mu1_pt;
float mu1_phi;
float mu1_eta;
float mu2_pt;
float mu2_phi;
float mu2_eta;
float mass_mumu;
float mass_mumug;
float mass_mumugg;

// photon
int   npho;
float pho1_pt;
float pho1_eta;
float pho1_phi;
float pho1_Sihih;
float pho1_ChIso;
float pho1_PhoIso;
float pho1_NeuHadIso;
float pho1_HoverE12;
int   pho1_parentage;

float pho2_pt;
float pho2_eta;
float pho2_phi;
float pho2_Sihih;
float pho2_ChIso;
float pho2_PhoIso;
float pho2_NeuHadIso;
float pho2_HoverE12;
int   pho2_parentage;

double analyze_Zgg::dR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = acos(cos(phi1 - phi2));
  double deta = eta1 - eta2;
  return sqrt(dphi*dphi + deta*deta);
}

bool analyze_Zgg::fidEtaPass(double Eta){
  double fabsEta = TMath::Abs(Eta);
  if( fabsEta > 2.5) return false;
  if( 1.4442 < fabsEta && fabsEta < 1.566) return false;
  return true;
}

int analyze_Zgg::phoRegion(double absEta){
  int region = 0;
  if( absEta >= 1.0  ) region++;
  if( absEta >= 1.479) region++;
  if( absEta >= 2.0  ) region++;
  if( absEta >= 2.2  ) region++;
  if( absEta >= 2.3  ) region++;
  if( absEta >= 2.4  ) region++;
  return region;
}

double analyze_Zgg::phoEffArea03ChHad(double peta){
  double eta = TMath::Abs(peta);
  static const double area[7] = {0.012, 0.010, 0.014, 0.012, 0.016, 0.020, 0.012};
  return area[phoRegion(eta)];
}

double analyze_Zgg::phoEffArea03NeuHad(double peta){
  double eta = TMath::Abs(peta);
  static const double area[7] = {0.030, 0.057, 0.039, 0.015, 0.024, 0.039, 0.072};
  return area[phoRegion(eta)];
}

double analyze_Zgg::phoEffArea03Pho(double peta){
  double eta = TMath::Abs(peta);
  static const double area[7] = {0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};
  return area[phoRegion(eta)];
}

// tight muon ID
bool analyze_Zgg::muID2012(int i){
  bool pass = true;

  if (!( (*muPt)[i] > 5 )) pass = false;
  if (!( fabs((*muEta)[i]) < 2.4 )) pass = false;

  bool is_global_muon = (*muType)[i] & 1<<1;
  bool is_pf_muon = (*muType)[i] & 1<<5;

  if (!( is_global_muon == true )) pass = false;
  if (!( is_pf_muon == true )) pass = false;

  if (!( (*muChi2NDF)[i] < 10 )) pass = false;
  if (!( (*muNumberOfValidMuonHits)[i] > 0 )) pass = false;
  if (!( (*muStations)[i] > 1 )) pass = false;
  if (!( (*muNumberOfValidPixelHits)[i] > 0 )) pass = false;
  if (!( (*muNumberOfValidTrkLayers)[i] > 5 )) pass = false;
  if (!( fabs((*muD0GV)[i]) < 0.2 )) pass = false;
  if (!( fabs((*muDzGV)[i]) < 0.5 )) pass = false;

  if (!( MuIsolation(i) < 0.12 )) pass = false;

  return pass;
}

float analyze_Zgg::MuIsolation(int i){
  float isolation;

  float sum_neu = (*muPFIsoR04_NH)[i] + (*muPFIsoR04_Pho)[i] - 0.5*(*muPFIsoR04_PU)[i];
  if(sum_neu < 0){
    sum_neu = 0.0;
  }
  float corriso = (*muPFIsoR04_CH)[i] + sum_neu;

  isolation = corriso / (*muPt)[i];

  return isolation;
}

void analyze_Zgg::clearVariables(){
  ele1_pt = 0;
  ele1_phi = 0;
  ele1_eta = 0;
  ele2_pt= 0;
  ele2_phi= 0;
  ele2_eta= 0;
  mass_ee = 0;
  mass_eeg = 0;
  mass_eegg = 0;

  mu1_pt = 0;
  mu1_phi = 0;
  mu1_eta = 0;
  mu2_pt= 0;
  mu2_phi= 0;
  mu2_eta= 0;
  mass_mumu = 0;
  mass_mumug = 0;
  mass_mumugg = 0;

  npho = 0;
  pho1_pt = 0;
  pho1_phi = 0;
  pho1_eta = 0;
  pho1_Sihih = 0;
  pho1_ChIso = 0;
  pho1_PhoIso = 0;
  pho1_NeuHadIso = 0;
  pho1_HoverE12 = 0;
  pho1_parentage = 0;

  pho2_pt = 0;
  pho2_phi = 0;
  pho2_eta = 0;
  pho2_Sihih = 0;
  pho2_ChIso = 0;
  pho2_PhoIso = 0;
  pho2_NeuHadIso = 0;
  pho2_HoverE12 = 0;
  pho2_parentage = 0;

  nVtx_ = 0;
  run_ = 0;
  event_ = 0;
}

void analyze_Zgg::Loop(TString name, TString fitVal, int iso)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  TFile* tmp = TFile::Open(name, "RECREATE");

  TTree* miniTree = new TTree("miniTree", "miniTree");

  miniTree->Branch("nVtx_", &nVtx_, "nVtx_/I");
  miniTree->Branch("run_", &run_, "run_/I");
  miniTree->Branch("event_", &event_, "event_/L");
  // electrons
  miniTree->Branch("ele1_pt", &ele1_pt, "ele1_pt/F");
  miniTree->Branch("ele1_phi", &ele1_phi, "ele1_phi/F");
  miniTree->Branch("ele1_eta", &ele1_eta, "ele1_eta/F");
  miniTree->Branch("ele2_pt", &ele2_pt, "ele2_pt/F");
  miniTree->Branch("ele2_phi", &ele2_phi, "ele2_phi/F");
  miniTree->Branch("ele2_eta", &ele2_eta, "ele2_eta/F");
  miniTree->Branch("mass_ee", &mass_ee, "mass_ee/F");
  miniTree->Branch("mass_eeg", &mass_eeg, "mass_eeg/F");
  miniTree->Branch("mass_eegg", &mass_eegg, "mass_eegg/F");
  // muons
  miniTree->Branch("mu1_pt", &mu1_pt, "mu1_pt/F");
  miniTree->Branch("mu1_phi", &mu1_phi, "mu1_phi/F");
  miniTree->Branch("mu1_eta", &mu1_eta, "mu1_eta/F");
  miniTree->Branch("mu2_pt", &mu2_pt, "mu2_pt/F");
  miniTree->Branch("mu2_phi", &mu2_phi, "mu2_phi/F");
  miniTree->Branch("mu2_eta", &mu2_eta, "mu2_eta/F");
  miniTree->Branch("mass_mumu", &mass_mumu, "mass_mumu/F");
  miniTree->Branch("mass_mumug", &mass_mumug, "mass_mumug/F");
  miniTree->Branch("mass_mumugg", &mass_mumugg, "mass_mumugg/F");

  // photons
  miniTree->Branch("npho", &npho, "npho/I");
  miniTree->Branch("pho1_pt", &pho1_pt, "pho1_pt/F");
  miniTree->Branch("pho1_phi", &pho1_phi, "pho1_phi/F");
  miniTree->Branch("pho1_eta", &pho1_eta, "pho1_eta/F");
  miniTree->Branch("pho1_Sihih", &pho1_Sihih, "pho1_Sihih/F");
  miniTree->Branch("pho1_ChIso", &pho1_ChIso, "pho1_ChIso/F");
  miniTree->Branch("pho1_PhoIso", &pho1_PhoIso, "pho1_PhoIso/F");
  miniTree->Branch("pho1_NeuHadIso", &pho1_NeuHadIso, "pho1_NeuHadIso/F");
  miniTree->Branch("pho1_HoverE12", &pho1_HoverE12, "pho1_HoverE12/F");
  miniTree->Branch("pho1_parentage", &pho1_parentage, "pho1_parentage/I");

  miniTree->Branch("pho2_pt", &pho2_pt, "pho2_pt/F");
  miniTree->Branch("pho2_phi", &pho2_phi, "pho2_phi/F");
  miniTree->Branch("pho2_eta", &pho2_eta, "pho2_eta/F");
  miniTree->Branch("pho2_Sihih", &pho2_Sihih, "pho2_Sihih/F");
  miniTree->Branch("pho2_ChIso", &pho2_ChIso, "pho2_ChIso/F");
  miniTree->Branch("pho2_PhoIso", &pho2_PhoIso, "pho2_PhoIso/F");
  miniTree->Branch("pho2_NeuHadIso", &pho2_NeuHadIso, "pho2_NeuHadIso/F");
  miniTree->Branch("pho2_HoverE12", &pho2_HoverE12, "pho2_HoverE12/F");
  miniTree->Branch("pho2_parentage", &pho2_parentage, "pho2_parentage/I");

  miniTree->Branch("totalWeight", &totalWeight, "totalWeight/F");
  miniTree->Branch("totalWeight_m5", &totalWeight_m5, "totalWeight_m5/F");
  miniTree->Branch("totalWeight_p5", &totalWeight_p5, "totalWeight_p5/F");
  miniTree->Branch("processCode", &processCode, "processCode/I");
  miniTree->Branch("xsWeight", &xsWeight, "xsWeight/F");
  miniTree->Branch("puWeight", &puWeight, "puWeight/F");
  miniTree->Branch("puWeight_m5", &puWeight_m5, "puWeight_m5/F");
  miniTree->Branch("puWeight_p5", &puWeight_p5, "puWeight_p5/F");

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);

    clearVariables();

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(jentry % 10000 == 0) cout << "Processed " << jentry
    				 << " events out of " <<nentries<< endl;

    //Trigger
    int trigger = 0;
    if(HLT[HLTIndex[13]]==1 || HLT[HLTIndex[14]]==1) trigger = 1;
    if(trigger!=1){
      continue;
    }

    // LOOP OVER EM OBJECTS
    vector <int> imuons;
    for (int imu = 0; imu < nMu; ++imu){

      // CHECK IF SATISFIES ID
      if( muID2012(imu) ) imuons.push_back(imu);
    }
    if ( imuons.size() != 2 ) continue;
    if( (*muPt)[imuons[0]] < 20. || (*muPt)[imuons[1]] < 10.) continue;
    TLorentzVector mu1, mu2;
    mu1.SetPtEtaPhiM((*muPt)[imuons[0]], (*muEta)[imuons[0]], (*muPhi)[imuons[0]], 0.1057);
    mu2.SetPtEtaPhiM((*muPt)[imuons[1]], (*muEta)[imuons[1]], (*muPhi)[imuons[1]], 0.1057);

    if( (mu1 + mu2).M() < 40 ) continue;

    // FILL MUON VARIABLES
    mu1_pt = (*muPt)[imuons[0]];
    mu1_eta = (*muEta)[imuons[0]];
    mu1_phi = (*muPhi)[imuons[0]];
    mu2_pt = (*muPt)[imuons[1]];
    mu2_eta = (*muEta)[imuons[1]];
    mu2_phi = (*muPhi)[imuons[1]];
    mass_mumu = ((mu1 + mu2).M());

    double ChHadIso_corr = 0; double NeuHadIso_corr = 0; double PhoIso_corr = 0;
    vector <double> chiso, neuiso, phoiso, sieie;
    // LOOP OVER FOR PHOTON
    vector <int> iphotons;
    for (int ipho = 0; ipho < nPho; ++ipho){
      int region = 1;
      if( fabs((*phoSCEta)[ipho])<1.479) region = 0;
      if( fabs((*phoSCEta)[ipho])>2.5) continue;
      if( fabs((*phoSCEta)[ipho])<1.57 && fabs((*phoSCEta)[ipho])>1.44) continue;
      if( (*phoEt)[ipho] < 15. ) continue;

      double Pho03ChHadIso = (*phoPFChIso)[ipho] - rho2012 * phoEffArea03ChHad((*phoSCEta)[ipho]);
      if(Pho03ChHadIso < 0) Pho03ChHadIso = 0;
      double Pho03NeuHadIso = (*phoPFNeuIso)[ipho] - rho2012 * phoEffArea03NeuHad((*phoSCEta)[ipho]);
      if(Pho03NeuHadIso < 0) Pho03NeuHadIso = 0;
      double Pho03PhoIso = (*phoPFPhoIso)[ipho] - rho2012 * phoEffArea03Pho((*phoSCEta)[ipho]);
      if(Pho03PhoIso < 0) Pho03PhoIso = 0;

      ChHadIso_corr = Pho03ChHadIso;
      NeuHadIso_corr = Pho03NeuHadIso - (*phoEt)[ipho] * photonID_RhoCorrR03NeuHadIso_1[region][1];
      PhoIso_corr = Pho03PhoIso - (*phoEt)[ipho] * photonID_RhoCorrR03PhoIso_1[region][1];

      bool phoSel = false;
      if(fitVal=="Sieie"){
	if(
	   ChHadIso_corr < photonID_RhoCorrChHadIso[region][iso] &&
	   NeuHadIso_corr < photonID_RhoCorrNeuHadIso[region][iso] &&
	   PhoIso_corr < photonID_RhoCorrPhoIso[region][iso] &&
	   (*phoEleVeto)[ipho] == 0 &&
	   (*phoHoverE12)[ipho] < photonID_HoverE[region][1] ) phoSel = true;
      }
      else if (fitVal=="ChIso"){
	if(
	   NeuHadIso_corr < photonID_RhoCorrNeuHadIso[region][0] &&
	   PhoIso_corr < photonID_RhoCorrPhoIso[region][0] &&
	   (*phoSigmaIEtaIEta)[ipho] < photonID_SigmaIhIh[region][iso] &&
	   (*phoEleVeto)[ipho] == 0 &&
	   (*phoHoverE12)[ipho] < photonID_HoverE[region][1] ) phoSel = true;
      }
      else if (fitVal=="PhoIso"){
	if(
	   ChHadIso_corr < photonID_RhoCorrChHadIso[region][iso] &&
	   NeuHadIso_corr < photonID_RhoCorrNeuHadIso[region][iso] &&
	   (*phoSigmaIEtaIEta)[ipho] < photonID_SigmaIEtaIEta[region][1] &&
	   (*phoEleVeto)[ipho] == 0 &&
	   (*phoHoverE12)[ipho] < photonID_HoverE[region][1] ) phoSel = true;
      } else {
	"ERROR!! HAVE TO PICK UP VARIAVLE TO FITTING";
      }
      if(!phoSel) continue;

      TLorentzVector pho;
      pho.SetPtEtaPhiE((*phoEt)[ipho], (*phoEta)[ipho], (*phoPhi)[ipho], (*phoE)[ipho]);
      // REJECT PHOTONS THAT ARE TOO CLOSE TO MUON
      if( pho.DeltaR(mu1) < 0.4 || pho.DeltaR(mu2) < 0.4 ) continue;

      sieie.push_back((*phoSigmaIEtaIEta)[ipho]);
      chiso.push_back(ChHadIso_corr);
      neuiso.push_back(NeuHadIso_corr);
      phoiso.push_back(PhoIso_corr);

      iphotons.push_back(ipho);
    } // loop over for photon
    if( iphotons.size() < 2 ) continue;
    npho = iphotons.size();

    TLorentzVector pho1, pho2;
    pho1.SetPtEtaPhiM((*phoEt)[iphotons[0]], (*phoEta)[iphotons[0]], (*phoPhi)[iphotons[0]], 0);
    pho2.SetPtEtaPhiM((*phoEt)[iphotons[1]], (*phoEta)[iphotons[1]], (*phoPhi)[iphotons[1]], 0);

    // FILL PHOTON VARIABLES INTO TREE
    pho1_pt = (*phoEt)[iphotons[0]];
    pho1_eta = (*phoSCEta)[iphotons[0]];
    pho1_phi = (*phoPhi)[iphotons[0]];
    pho1_Sihih = sieie[0];
    pho1_ChIso = chiso[0];
    pho1_PhoIso = phoiso[0];
    pho1_NeuHadIso = neuiso[0];
    pho1_HoverE12 = (*phoHoverE12)[iphotons[0]];
    pho2_pt = (*phoEt)[iphotons[1]];
    pho2_eta = (*phoSCEta)[iphotons[1]];
    pho2_phi = (*phoPhi)[iphotons[1]];
    pho2_Sihih = sieie[1];
    pho2_ChIso = chiso[1];
    pho2_PhoIso = phoiso[1];
    pho2_NeuHadIso = neuiso[1];
    pho2_HoverE12 = (*phoHoverE12)[iphotons[1]];
    mass_mumug = (mu1 + mu2 + pho1).M();
    mass_mumugg = (mu1 + mu2 + pho1 + pho2).M();

    nVtx_ = nVtx;
    run_ = run;
    event_ = event;

    totalWeight = 1.0;
    totalWeight_m5 = 1.0;
    totalWeight_p5 = 1.0;
    puWeight    = 1.0;
    puWeight_m5 = 1.0;
    puWeight_p5 = 1.0;

    miniTree->Fill();
  } // loop over for events
  miniTree->Write();
  tmp->Close();
}

int main(int argc, char* argv[]){

  TString fitVal = argv[1];
  int iso = atoi(argv[2]);

  for (int i = 0; i < argc; i++ ) {

    analyze_Zgg t;
    t.Loop("muon_Zgg_"+fitVal+"_"+Loose[iso]+".root", fitVal, iso);
  }

  return 0;
}
