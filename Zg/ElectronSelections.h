Float_t CorrectPFCombined(Float_t pfCombined, Float_t Eta, Float_t rho2012) {

  if (fabs(Eta) < 1.0)
    return max(pfCombined - 0.208*rho2012, 0.);
  else if (fabs(Eta) > 1.0 && fabs(Eta) < 1.5)
    return max(pfCombined - 0.209*rho2012, 0.);
  else if (fabs(Eta) > 1.5 && fabs(Eta) < 2.0)
    return max(pfCombined - 0.115*rho2012, 0.);
  else if (fabs(Eta) > 2.0 && fabs(Eta) < 2.2)
    return max(pfCombined - 0.143*rho2012, 0.);
  else if (fabs(Eta) > 2.2 && fabs(Eta) < 2.3)
    return max(pfCombined - 0.183*rho2012, 0.);
  else if (fabs(Eta) > 2.3 && fabs(Eta) < 2.4)
    return max(pfCombined - 0.194*rho2012, 0.);
  else if (fabs(Eta) > 2.4)
    return max(pfCombined - 0.261*rho2012, 0.);
  else
    return pfCombined;
}

Float_t CorrectPFCharged(float pfCharged, float SCEta, float rho2012) {

  if (fabs(SCEta) < 1.0)
    return max(pfCharged - 0.012*rho2012, 0.);
  else if (fabs(SCEta) > 1.0 && fabs(SCEta) < 1.4442)
    return max(pfCharged - 0.010*rho2012, 0.);
  else if (fabs(SCEta) > 1.566 && fabs(SCEta) < 2.0)
    return max(pfCharged - 0.014*rho2012, 0.);
  else if (fabs(SCEta) > 2.0 && fabs(SCEta) < 2.2)
    return max(pfCharged - 0.012*rho2012, 0.);
  else if (fabs(SCEta) > 2.2 && fabs(SCEta) < 2.3)
    return max(pfCharged - 0.016*rho2012, 0.);
  else if (fabs(SCEta) > 2.3 && fabs(SCEta) < 2.4)
    return max(pfCharged - 0.020*rho2012, 0.);
  else if (fabs(SCEta) > 2.4 && fabs(SCEta) < 2.5)
    return max(pfCharged - 0.012*rho2012, 0.);
  else
    return pfCharged;
}

Float_t CorrectPFNeutral(float pfNeutral, float SCEta, float rho2012) {

  if (fabs(SCEta) < 1.0)
    return max(pfNeutral - 0.030*rho2012, 0.);
  else if (fabs(SCEta) > 1.0 && fabs(SCEta) < 1.4442)
    return max(pfNeutral - 0.057*rho2012, 0.);
  else if (fabs(SCEta) > 1.566 && fabs(SCEta) < 2.0)
    return max(pfNeutral - 0.039*rho2012, 0.);
  else if (fabs(SCEta) > 2.0 && fabs(SCEta) < 2.2)
    return max(pfNeutral - 0.015*rho2012, 0.);
  else if (fabs(SCEta) > 2.2 && fabs(SCEta) < 2.3)
    return max(pfNeutral - 0.024*rho2012, 0.);
  else if (fabs(SCEta) > 2.3 && fabs(SCEta) < 2.4)
    return max(pfNeutral - 0.039*rho2012, 0.);
  else if (fabs(SCEta) > 2.4 && fabs(SCEta) < 2.5)
    return max(pfNeutral - 0.072*rho2012, 0.);
  else
    return pfNeutral;
}

Float_t CorrectPFPhoton(float pfPhoton, float SCEta, float rho2012) {

  if (fabs(SCEta) < 1.0)
    return max(pfPhoton - 0.148*rho2012, 0.);
  else if (fabs(SCEta) > 1.0 && fabs(SCEta) < 1.4442)
    return max(pfPhoton - 0.130*rho2012, 0.);
  else if (fabs(SCEta) > 1.566 && fabs(SCEta) < 2.0)
    return max(pfPhoton - 0.112*rho2012, 0.);
  else if (fabs(SCEta) > 2.0 && fabs(SCEta) < 2.2)
    return max(pfPhoton - 0.216*rho2012, 0.);
  else if (fabs(SCEta) > 2.2 && fabs(SCEta) < 2.3)
    return max(pfPhoton - 0.262*rho2012, 0.);
  else if (fabs(SCEta) > 2.3 && fabs(SCEta) < 2.4)
    return max(pfPhoton - 0.260*rho2012, 0.);
  else if (fabs(SCEta) > 2.4 && fabs(SCEta) < 2.5)
    return max(pfPhoton - 0.266*rho2012, 0.);
  else
    return pfPhoton;
}

void select_electrons(TreeReader &data, vector<int> &accepted) {

  accepted.clear();
  // Example of electron/positron selection.

  // load relevant branches from TTree/TChain
  Int_t    nEle             = data.GetInt("nEle");
  Float_t* elePt            = data.GetPtrFloat("elePt");
  Float_t* eleSCEta         = data.GetPtrFloat("eleSCEta");
  Float_t* eleSigmaIEtaIEta = data.GetPtrFloat("eleSigmaIEtaIEta");
  Float_t* eleHoverE        = data.GetPtrFloat("eleHoverE");
  Int_t*   eleMissHits      = data.GetPtrInt("eleMissHits");
  Float_t* eleIsoEcalDR03   = data.GetPtrFloat("eleIsoEcalDR03");
  Float_t* eleIsoHcalDR03   = data.GetPtrFloat("eleIsoHcalDR03");
  Float_t* eleIsoTrkDR03    = data.GetPtrFloat("eleIsoTrkDR03");
  Float_t* eleIDMVATrig     = data.GetPtrFloat("eleIDMVATrig");
  Float_t* elePFChIso03     = data.GetPtrFloat("elePFChIso03");
  Float_t* elePFPhoIso03    = data.GetPtrFloat("elePFPhoIso03");
  Float_t* elePFNeuIso03    = data.GetPtrFloat("elePFNeuIso03");
  Float_t  rho2012          = data.GetFloat("rho2012");

  // loop over lepton candidates
  for (int i = 0; i < nEle; ++i) {

    if (elePt[i] < 10.) continue;

    // MVA Preselection
    if (fabs(eleSCEta[i]) < 1.479 && eleSigmaIEtaIEta[i] > 0.014) continue;
    if (fabs(eleSCEta[i]) < 1.479 && eleHoverE[i] > 0.15) continue;
    if (fabs(eleSCEta[i]) < 1.479 && eleMissHits[i] > 0) continue;
    if (fabs(eleSCEta[i]) < 1.479 && eleIsoEcalDR03[i]/elePt[i] > 0.1) continue;
    if (fabs(eleSCEta[i]) < 1.479 && eleIsoHcalDR03[i]/elePt[i] > 0.1) continue;
    if (fabs(eleSCEta[i]) < 1.479 && eleIsoTrkDR03[i]/elePt[i] > 0.1) continue;

    if (fabs(eleSCEta[i]) > 1.479 && eleSigmaIEtaIEta[i] > 0.035) continue;
    if (fabs(eleSCEta[i]) > 1.479 && eleHoverE[i] >= 0.10) continue;
    if (fabs(eleSCEta[i]) > 1.479 && eleMissHits[i] > 0) continue;
    if (fabs(eleSCEta[i]) > 1.479 && eleIsoEcalDR03[i]/elePt[i] > 0.15) continue;
    if (fabs(eleSCEta[i]) > 1.479 && eleIsoHcalDR03[i]/elePt[i] > 0.15) continue;
    if (fabs(eleSCEta[i]) > 1.479 && eleIsoTrkDR03[i]/elePt[i] > 0.15) continue;

    // MVA ID+ISO  selection
    if (elePt[i] < 20 && eleIDMVATrig[i] < -0.9) continue;
    if (elePt[i] > 20 && eleIDMVATrig[i] < -0.5) continue;

    if ((elePFChIso03[i] + CorrectPFCombined(elePFPhoIso03[i]+elePFNeuIso03[i], eleSCEta[i], rho2012))/elePt[i] > 0.3) 
      continue;

    accepted.push_back(i);
  }
}

Float_t geteIDEA(Int_t i, Float_t SCEta) {

  Float_t EffectiveArea = 0;

  if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.013 + 0.122;
  if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.021 + 0.147 ;
  if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.013 + 0.055;
  if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.010 + 0.106;
  if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.024 + 0.138;
  if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.020 + 0.221;
  if (fabs(SCEta) >= 2.4 ) EffectiveArea = 0.019 + 0.211;

  return EffectiveArea;
}

void eID2012(TreeReader &data, vector<int> &accepted, Int_t idWP) {

  accepted.clear();

  Float_t dEtaInCut_B[4]    = {0.007, 0.007, 0.004, 0.004};
  Float_t dPhiInCut_B[4]    = {0.8,   0.15,  0.06,  0.03};
  Float_t sIeIeCut_B[4]     = {0.01,  0.01,  0.01,  0.01};
  Float_t HoECut_B[4]       = {0.15,  0.12,  0.12,  0.12};
  Float_t d0Cut_B[4]        = {0.04,  0.02,  0.02,  0.02};
  Float_t dzCut_B[4]        = {0.2,   0.2,   0.1,   0.1};
  Float_t diffEpCut_B[4]    = {999.,  0.05,  0.05,  0.05};
  Float_t relPFisoCut_B[4]  = {0.15,  0.15,  0.15,  0.10};
  Int_t VtxFitProbCut_B[4]  = {999,   1,     1,     1};
  Int_t missHitsCut_B[4]    = {99,    1,     1,     0};

  Float_t dEtaInCut_EC[4]   = {0.01 , 0.009,  0.007,  0.005};
  Float_t dPhiInCut_EC[4]   = {0.7,   0.10,   0.03,   0.02};
  Float_t sIeIeCut_EC[4]    = {0.03,  0.03,   0.03,   0.03};
  Float_t HoECut_EC[4]      = {999,   0.10,   0.10,   0.10};
  Float_t d0Cut_EC[4]       = {0.04,  0.02,   0.02,   0.02};
  Float_t dzCut_EC[4]       = {0.2,   0.2,    0.1,    0.1};
  Float_t diffEpCut_EC[4]   = {999,   0.05,   0.05,   0.05};
  Float_t relPFisoCut_EC[4] = {0.15,  0.15,   0.15,   0.10};
  Int_t VtxFitProbCut_EC[4] = {999,   1,      1,      1};
  Int_t missHitsCut_EC[4]   = {99,    1,      1,      0};

  Int_t    nEle             = data.GetInt("nEle");
  Float_t* elePt            = data.GetPtrFloat("elePt");
  Float_t* eleSCEta         = data.GetPtrFloat("eleSCEta");
  Float_t* eleSigmaIEtaIEta = data.GetPtrFloat("eleSigmaIEtaIEta");
  Float_t* eleHoverE        = data.GetPtrFloat("eleHoverE");
  Int_t*   eleMissHits      = data.GetPtrInt("eleMissHits");
  Int_t*   eleConvVtxFit    = data.GetPtrInt("eleConvVtxFit");
  Float_t* eledEtaAtVtx     = data.GetPtrFloat("eledEtaAtVtx");
  Float_t* eledPhiAtVtx     = data.GetPtrFloat("eledPhiAtVtx");
  Float_t* eleEcalEn        = data.GetPtrFloat("eleEcalEn");
  Float_t* elePin           = data.GetPtrFloat("elePin");
  Float_t* elePFChIso03     = data.GetPtrFloat("elePFChIso03");
  Float_t* elePFPhoIso03    = data.GetPtrFloat("elePFPhoIso03");
  Float_t* elePFNeuIso03    = data.GetPtrFloat("elePFNeuIso03");
  Float_t  rho2012          = data.GetFloat("rho2012");

  vector<float>* eleD0Vtx   = data.GetPtrVectorFloat("eleD0Vtx", nEle);
  vector<float>* eleDzVtx   = data.GetPtrVectorFloat("eleDzVtx", nEle);

  Double_t elePFIsoSumRel = 0.;
  
  for (int i = 0; i < nEle; ++i) {

    if (elePt[i] < 10.) continue;

    Float_t EffectiveArea = geteIDEA(i, eleSCEta[i]);
    elePFIsoSumRel = (elePFChIso03[i] + TMath::Max((Double_t) 0., (Double_t) elePFPhoIso03[i] + elePFNeuIso03[i] - EffectiveArea*rho2012))/elePt[i];
    
    if (fabs(eleSCEta[i]) < 1.4442) {
      if (fabs(eledEtaAtVtx[i]) > dEtaInCut_B[idWP]) continue;
      if (fabs(eledPhiAtVtx[i]) > dPhiInCut_B[idWP]) continue;
      if (eleSigmaIEtaIEta[i] > sIeIeCut_B[idWP]) continue;
      if (eleHoverE[i] > HoECut_B[idWP]) continue;
      if (fabs(eleD0Vtx[i][0]) > d0Cut_B[idWP]) continue;
      if (fabs(eleDzVtx[i][0]) > dzCut_B[idWP]) continue;
      if (fabs(1./eleEcalEn[i] - 1./elePin[i]) > diffEpCut_B[idWP]) continue;
      if (eleConvVtxFit[i] == VtxFitProbCut_B[idWP]) continue;
      if (eleMissHits[i] > missHitsCut_B[idWP]) continue;
      if (elePFIsoSumRel > relPFisoCut_B[idWP]) continue;
    } else {
      if (fabs(eledEtaAtVtx[i]) > dEtaInCut_EC[idWP]) continue;
      if (fabs(eledPhiAtVtx[i]) > dPhiInCut_EC[idWP]) continue;
      if (eleSigmaIEtaIEta[i] > sIeIeCut_EC[idWP]) continue;
      if (eleHoverE[i] > HoECut_EC[idWP]) continue;
      if (fabs(eleD0Vtx[i][0]) > d0Cut_EC[idWP]) continue;
      if (fabs(eleDzVtx[i][0]) > dzCut_EC[idWP]) continue;
      if (fabs(1./eleEcalEn[i] - 1./elePin[i]) > diffEpCut_EC[idWP]) continue;
      if (eleConvVtxFit[i] == VtxFitProbCut_EC[idWP]) continue;
      if (eleMissHits[i] > missHitsCut_EC[idWP]) continue;
      if (elePFIsoSumRel > relPFisoCut_EC[idWP]) continue;
    }

    accepted.push_back(i);
  }

}


