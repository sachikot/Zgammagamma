#!/bin/csh
#
source /cvmfs/cms.cern.ch/cmsset_default.csh
cd /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/CMSSW_7_0_9/src
cmsenv
cd ${_CONDOR_SCRATCH_DIR}

cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_SKIM/Muon/skimmer_zgg .
cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_SKIM/Muon/Pileup*root .
cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_SKIM/Muon/*.h .
cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_SKIM/Muon/*.cc .
cp "/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_SKIM/Muon/${2}" .

./skimmer_zgg ${1} ${2}

echo "DONE!"




