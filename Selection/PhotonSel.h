const int    photonID_IsConv[2][3]                = { {0, 0, 0} ,             {0, 0, 0}             };
const double photonID_HoverE[2][3]                = { {0.05, 0.05, 0.05} ,    {0.05, 0.05, 0.05}    };
const double photonID_SigmaIEtaIEta[2][3]         = { {0.012, 0.011, 0.011} , {0.034, 0.033, 0.031} };
const double photonID_RhoCorrR03ChHadIso[2][3]    = { {2.6, 1.5, 0.7} ,       {2.3, 1.2, 0.5}       };
const double photonID_RhoCorrR03NeuHadIso_0[2][3] = { {3.5, 1.0, 0.4} ,       {2.9, 1.5, 1.5}       };
const double photonID_RhoCorrR03NeuHadIso_1[2][3] = { {0.04, 0.04, 0.04} ,    {0.04, 0.04, 0.04}    };
const double photonID_RhoCorrR03PhoIso_0[2][3]    = { {1.3, 0.7, 0.5} ,       {999, 1.0, 1.0}       };
const double photonID_RhoCorrR03PhoIso_1[2][3]    = { {0.005, 0.005, 0.005} , {0.005, 0.005, 0.005} };

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxeleESEffSigmaRR = 1;
const Int_t kMaxphoESEffSigmaRR = 1;
const Int_t kMaxnPFPho = 67;
const Int_t kMaxPFPhoEt = 1;
const Int_t kMaxPFPhoEta = 1;
const Int_t kMaxPFPhoPhi = 1;
const Int_t kMaxPFPhoType = 1;
const Int_t kMaxPFPhoIso = 1;

double dR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = acos(cos(phi1 - phi2));
  double deta = eta1 - eta2;
  return sqrt(dphi*dphi + deta*deta);
}

bool fidEtaPass(double Eta){
  double fabsEta = TMath::Abs(Eta);
  if( fabsEta > 2.5) return false;
  if( 1.4442 < fabsEta && fabsEta < 1.566) return false;
  return true;
}

int phoRegion(double absEta){
  int region = 0;
  if( absEta >= 1.0  ) region++;
  if( absEta >= 1.479) region++;
  if( absEta >= 2.0  ) region++;
  if( absEta >= 2.2  ) region++;
  if( absEta >= 2.3  ) region++;
  if( absEta >= 2.4  ) region++;
  return region;
}

double phoEffArea03ChHad(double peta){
  double eta = TMath::Abs(peta);
  static const double area[7] = {0.012, 0.010, 0.014, 0.012, 0.016, 0.020, 0.012};
  return area[phoRegion(eta)];
}

double phoEffArea03NeuHad(double peta){
  double eta = TMath::Abs(peta);
  static const double area[7] = {0.030, 0.057, 0.039, 0.015, 0.024, 0.039, 0.072};
  return area[phoRegion(eta)];
}

double phoEffArea03Pho(double peta){
  double eta = TMath::Abs(peta);
  static const double area[7] = {0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};
  return area[phoRegion(eta)];
}

bool PhotonSel(int ipho){
  int region = 1;
  if(fabs((*phoSCEta)[ipho])<1.479) region = 0;
  double Pho03ChHadIso  = (*phoPFChIso)[ipho] - rho2012 * phoEffArea03ChHad((*phoSCEta)[ipho]);
  double Pho03NeuHadIso = (*phoPFNeuIso)[ipho] - rho2012 * phoEffArea03NeuHad((*phoSCEta)[ipho]);
  double Pho03PhoIso    = (*phoPFPhoIso)[ipho] - rho2012 * phoEffArea03Pho((*phoSCEta)[ipho]);

  if(Pho03ChHadIso < 0)  Pho03ChHadIso = 0;
  if(Pho03NeuHadIso < 0) Pho03NeuHadIso = 0;
  if(Pho03PhoIso < 0)    Pho03PhoIso = 0;

  double ChHadIso_corr  = Pho03ChHadIso;
  double NeuHadIso_corr = Pho03NeuHadIso - (*phoEt)[ipho] * photonID_RhoCorrR03NeuHadIso_1[region][1];
  double PhoIso_corr    = Pho03PhoIso - (*phoEt)[ipho] * photonID_RhoCorrR03PhoIso_1[region][1];

  // photon selection ID
  bool phoSel = 
    ChHadIso_corr < photonID_RhoCorrR03ChHadIso[region][1] &&
    NeuHadIso_corr < photonID_RhoCorrR03NeuHadIso_0[region][1] &&
    PhoIso_corr < photonID_RhoCorrR03PhoIso_0[region][1] &&
    (*phoSigmaIEtaIEta)[ipho] < photonID_SigmaIEtaIEta[region][1] &&
    (*phoEleVeto)[ipho] == 0 &&
    (*phoHoverE12)[ipho] < photonID_HoverE[region][1];

  return phoSel;
}
