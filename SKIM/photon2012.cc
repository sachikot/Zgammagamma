const int    photonID_IsConv[2][3]                = { {0, 0, 0} ,             {0, 0, 0}             };
const double photonID_HoverE[2][3]                = { {0.05, 0.05, 0.05} ,    {0.05, 0.05, 0.05}    };
const double photonID_SigmaIEtaIEta[2][3]         = { {0.012, 0.011, 0.011} , {0.034, 0.033, 0.031} };
const double photonID_RhoCorrR03ChHadIso[2][3]    = { {2.6, 1.5, 0.7} ,       {2.3, 1.2, 0.5}       };
const double photonID_RhoCorrR03NeuHadIso_0[2][3] = { {3.5, 1.0, 0.4} ,       {2.9, 1.5, 1.5}       };
const double photonID_RhoCorrR03NeuHadIso_1[2][3] = { {0.04, 0.04, 0.04} ,    {0.04, 0.04, 0.04}    };
const double photonID_RhoCorrR03PhoIso_0[2][3]    = { {1.3, 0.7, 0.5} ,       {999, 1.0, 1.0}       };
const double photonID_RhoCorrR03PhoIso_1[2][3]    = { {0.005, 0.005, 0.005} , {0.005, 0.005, 0.005} };

double skimmer_zgg::dR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = acos(cos(phi1 - phi2));
  double deta = eta1 - eta2;
  return sqrt(dphi*dphi + deta*deta);
}

bool skimmer_zgg::fidEtaPass(double Eta){
  double fabsEta = TMath::Abs(Eta);
  if( fabsEta > 2.5) return false;
  if( 1.4442 < fabsEta && fabsEta < 1.566) return false;
  return true;
}

int skimmer_zgg::phoRegion(double absEta){
  int region = 0;
  if( absEta >= 1.0  ) region++;
  if( absEta >= 1.479) region++;
  if( absEta >= 2.0  ) region++;
  if( absEta >= 2.2  ) region++;
  if( absEta >= 2.3  ) region++;
  if( absEta >= 2.4  ) region++;
  return region;
}

double skimmer_zgg::phoEffArea03ChHad(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.012, 0.010, 0.014, 0.012, 0.016, 0.020, 0.012};
  return area[phoRegion(eta)];
}

double skimmer_zgg::phoEffArea03NeuHad(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.030, 0.057, 0.039, 0.015, 0.024, 0.039, 0.072};
  return area[phoRegion(eta)];
}

double skimmer_zgg::phoEffArea03Pho(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};
  return area[phoRegion(eta)];
}

bool skimmer_zgg::passPhotonID(int phoInd, int pho_ID_ind = 0) {
  // phoInd - index of the photon in the tree
  // pho_ID_ind: 0 -- loose, 1 -- medium, 2 -- tight
  double eta = (*phoSCEta)[phoInd];
  double et = (*phoEt)[phoInd];
  
  // rho corrected isolations
  double Pho03ChHadIso =     (*phoPFChIso)[phoInd]   - rho2012 * phoEffArea03ChHad(eta);
  double Pho03ChHadSCRIso =  (*phoSCRChIso)[phoInd]  - rho2012 * phoEffArea03ChHad(eta);
  double Pho03NeuHadIso =    (*phoPFNeuIso)[phoInd]  - rho2012 * phoEffArea03NeuHad(eta);
  double Pho03PhoIso =       (*phoPFPhoIso)[phoInd]  - rho2012 * phoEffArea03Pho(eta);
  double Pho03PhoSCRIso =    (*phoSCRPhoIso)[phoInd] - rho2012 * phoEffArea03Pho(eta);
  
  // // manual spike cleaning
  // if (dR((*phoEta)[phoInd], (*phoPhi)[phoInd], -1.76, 1.37) < 0.05) return false;
  // if (dR((*phoEta)[phoInd], (*phoPhi)[phoInd],  2.37, 2.69) < 0.05) return false;
  
  int region = 0; //barrel
  if(TMath::Abs( eta )>1.566) region = 1; //endcap
  bool phoPresel = fidEtaPass( eta ) &&
    et > 10 && //  Et cut
    //    (*phohasPixelSeed)[phoInd] == 0 && // add if needed
    //    (*phoEleVeto)[phoInd] == 0 && // add if needed
    //    (*phoIsConv)[phoInd] == photonID_IsConv[region][pho_ID_ind] && // sachiko:add for photon candidate 
    (*phoHoverE12)[phoInd] < photonID_HoverE[region][pho_ID_ind];// &&
    //Pho03NeuHadIso < (photonID_RhoCorrR03NeuHadIso_0[region][pho_ID_ind] + et * photonID_RhoCorrR03NeuHadIso_1[region][pho_ID_ind]);
  return phoPresel;
}
