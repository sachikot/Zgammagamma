#!/bin/csh
#
/bin/rm -f condor_${2}
cat > condor_${2} << +EOF

universe = vanilla
Executable = /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_SKIM/Muon/CrossCheck/makeProduction.csh
output     = /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_SKIM/Muon/CrossCheck/an_${2}.out
error      = /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_SKIM/Muon/CrossCheck/an_${2}.err
log        = /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_SKIM/Muon/CrossCheck/an_${2}.log
Requirements   = (Memory >= 499 && OpSys == "LINUX" && (Arch =="x86_64" || Arch =="INTEL"))
Should_Transfer_Files = YES
When_To_Transfer_Output = ON_EXIT
Arguments = ${1} ${2}
Queue 1

+EOF

condor_submit condor_${2}
