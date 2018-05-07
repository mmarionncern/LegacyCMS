#!/bin/bash

root -b <<EOF &
.L makeFlipPlots.C+
//splitTree("","MCSamples.root");
splitTree("MVAM","MCSamples.root");
splitTree("MVAVT","MCSamples.root"); 


gROOT->ProcessLine(".q");
EOF

root -b <<EOF &
.L makeFlipPlots.C+
splitTree("","DoubleEG_BCDEFGH.root");  
gROOT->ProcessLine(".q");
EOF

root -b <<EOF &
.L makeFlipPlots.C+
splitTree("MVAM","DoubleEG_BCDEFGH.root"); 
splitTree("MVAVT","DoubleEG_BCDEFGH.root");
gROOT->ProcessLine(".q");
EOF
