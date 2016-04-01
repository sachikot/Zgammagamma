float geteIDEA(int i, float SCEta){
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

// loose electron ID
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification#Barrel_Cuts_eta_supercluster_1_4
bool eID2012(int i){

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
