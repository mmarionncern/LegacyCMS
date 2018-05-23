#!/bin/bash


source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.04/x86_64-slc5-gcc43-dbg/root/bin/thisroot.sh



cd /afs/cern.ch/user/m/mmarionn/workspace/private/METStudies
root -b <<EOF
.L RecoilCorrector.cc+
.L ComputeResolution.C+

apf=1;
loadMETTypes("$1","$2");
CompaResoCurves("$3",$4,$5,$6,$7,$8);

EOF
