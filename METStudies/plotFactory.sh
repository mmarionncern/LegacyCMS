#!/bin/bash

root -b <<EOF
.L RecoilCorrector.cc+
.L ComputeResolution.C+

apf=1;
loadMETTypes("$1","$2");
CompaRespCurves(1);

EOF

root -b <<EOF
.L RecoilCorrector.cc+
.L ComputeResolution.C+

apf=1;
loadMETTypes("$1","$2");
CompaResoCurves("para",0,1,0,1,0);

EOF

root -b <<EOF
.L RecoilCorrector.cc+
.L ComputeResolution.C+

apf=1;
loadMETTypes("$1","$2");
CompaResoCurves("para",1,1,0,1,0);

EOF

#root -b <<EOF
#.L RecoilCorrector.cc+
#.L ComputeResolution.C+

#apf=1;
#loadMETTypes("$1","$2");
#CompaResoCurves("para",1,1,1,1,0);

#EOF

root -b <<EOF
.L RecoilCorrector.cc+
.L ComputeResolution.C+

apf=1;
loadMETTypes("$1","$2");
CompaResoCurves("para",0,1,1,1,1);

EOF

#root -b <<EOF
#.L RecoilCorrector.cc+
#.L ComputeResolution.C+

#apf=1;
#loadMETTypes("$1","$2");
#CompaResoCurves("para",0,1,0,1,1);

#EOF

root -b <<EOF
.L RecoilCorrector.cc+
.L ComputeResolution.C+

apf=1;
loadMETTypes("$1","$2");
CompaResoCurves("perp",0,1,0,1,0);

EOF

#root -b <<EOF
#.L RecoilCorrector.cc+
#.L ComputeResolution.C+

#apf=1;
#loadMETTypes("$1","$2");
#CompaResoCurves("perp",1,1,0,1,0);

#EOF

root -b <<EOF
.L RecoilCorrector.cc+
.L ComputeResolution.C+

apf=1;
loadMETTypes("$1","$2");
CompaResoCurves("perp",1,1,1,1,0);

EOF

root -b <<EOF
.L RecoilCorrector.cc+
.L ComputeResolution.C+

apf=1;
loadMETTypes("$1","$2");
CompaResoCurves("perp",0,1,1,1,1);

EOF

#root -b <<EOF
#.L RecoilCorrector.cc+
#.L ComputeResolution.C+

#apf=1;
#loadMETTypes("$1","$2");
#CompaResoCurves("perp",0,1,0,1,1);

#EOF