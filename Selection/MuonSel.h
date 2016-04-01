bool muID2012(int i){
  bool pass = true;

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

float MuIsolation(int i){
  float isolation;

  float sum_neu = (*muPFIsoR04_NH)[i] + (*muPFIsoR04_Pho)[i] - 0.5*(*muPFIsoR04_PU)[i];
  if(sum_neu < 0){
    sum_neu = 0.0;
  }
  float corriso = (*muPFIsoR04_CH)[i] + sum_neu;

  isolation = corriso / (*muPt)[i];

  return isolation;
}
