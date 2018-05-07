#!/bin/bash

#root -l MCSamplesFlipTrees.root > tmp <<EOF
#.ls
#.q
#EOF


TREES=(
l1Pt_10_l1Eta_0_l2Pt_10_l2Eta_0_SS
l1Pt_10_l1Eta_0_l2Pt_10_l2Eta_0_OS
l1Pt_10_l1Eta_0_l2Pt_10_l2Eta_0.8_SS
l1Pt_10_l1Eta_0_l2Pt_10_l2Eta_0.8_OS
l1Pt_10_l1Eta_0_l2Pt_10_l2Eta_1.479_SS
l1Pt_10_l1Eta_0_l2Pt_10_l2Eta_1.479_OS
l1Pt_10_l1Eta_0_l2Pt_20_l2Eta_0_SS
l1Pt_10_l1Eta_0_l2Pt_20_l2Eta_0_OS
l1Pt_10_l1Eta_0_l2Pt_20_l2Eta_0.8_SS
l1Pt_10_l1Eta_0_l2Pt_20_l2Eta_0.8_OS
l1Pt_10_l1Eta_0_l2Pt_20_l2Eta_1.479_SS
l1Pt_10_l1Eta_0_l2Pt_20_l2Eta_1.479_OS
l1Pt_10_l1Eta_0_l2Pt_50_l2Eta_0_SS
l1Pt_10_l1Eta_0_l2Pt_50_l2Eta_0_OS
l1Pt_10_l1Eta_0_l2Pt_50_l2Eta_0.8_SS
l1Pt_10_l1Eta_0_l2Pt_50_l2Eta_0.8_OS
l1Pt_10_l1Eta_0_l2Pt_50_l2Eta_1.479_SS
l1Pt_10_l1Eta_0_l2Pt_50_l2Eta_1.479_OS
l1Pt_10_l1Eta_0_l2Pt_100_l2Eta_0_SS
l1Pt_10_l1Eta_0_l2Pt_100_l2Eta_0_OS
l1Pt_10_l1Eta_0_l2Pt_100_l2Eta_0.8_SS
l1Pt_10_l1Eta_0_l2Pt_100_l2Eta_0.8_OS
l1Pt_10_l1Eta_0_l2Pt_100_l2Eta_1.479_SS
l1Pt_10_l1Eta_0_l2Pt_100_l2Eta_1.479_OS
l1Pt_10_l1Eta_0_l2Pt_200_l2Eta_0_SS
l1Pt_10_l1Eta_0_l2Pt_200_l2Eta_0_OS
l1Pt_10_l1Eta_0_l2Pt_200_l2Eta_0.8_SS
l1Pt_10_l1Eta_0_l2Pt_200_l2Eta_0.8_OS
l1Pt_10_l1Eta_0_l2Pt_200_l2Eta_1.479_SS
l1Pt_10_l1Eta_0_l2Pt_200_l2Eta_1.479_OS
l1Pt_10_l1Eta_0.8_l2Pt_10_l2Eta_0.8_SS
l1Pt_10_l1Eta_0.8_l2Pt_10_l2Eta_0.8_OS
l1Pt_10_l1Eta_0.8_l2Pt_10_l2Eta_1.479_SS
l1Pt_10_l1Eta_0.8_l2Pt_10_l2Eta_1.479_OS
l1Pt_10_l1Eta_0.8_l2Pt_20_l2Eta_0.8_SS
l1Pt_10_l1Eta_0.8_l2Pt_20_l2Eta_0.8_OS
l1Pt_10_l1Eta_0.8_l2Pt_20_l2Eta_1.479_SS
l1Pt_10_l1Eta_0.8_l2Pt_20_l2Eta_1.479_OS
l1Pt_10_l1Eta_0.8_l2Pt_50_l2Eta_0.8_SS
l1Pt_10_l1Eta_0.8_l2Pt_50_l2Eta_0.8_OS
l1Pt_10_l1Eta_0.8_l2Pt_50_l2Eta_1.479_SS
l1Pt_10_l1Eta_0.8_l2Pt_50_l2Eta_1.479_OS
l1Pt_10_l1Eta_0.8_l2Pt_100_l2Eta_0.8_SS
l1Pt_10_l1Eta_0.8_l2Pt_100_l2Eta_0.8_OS
l1Pt_10_l1Eta_0.8_l2Pt_100_l2Eta_1.479_SS
l1Pt_10_l1Eta_0.8_l2Pt_100_l2Eta_1.479_OS
l1Pt_10_l1Eta_0.8_l2Pt_200_l2Eta_0.8_SS
l1Pt_10_l1Eta_0.8_l2Pt_200_l2Eta_0.8_OS
l1Pt_10_l1Eta_0.8_l2Pt_200_l2Eta_1.479_SS
l1Pt_10_l1Eta_0.8_l2Pt_200_l2Eta_1.479_OS
l1Pt_10_l1Eta_1.479_l2Pt_10_l2Eta_1.479_SS
l1Pt_10_l1Eta_1.479_l2Pt_10_l2Eta_1.479_OS
l1Pt_10_l1Eta_1.479_l2Pt_20_l2Eta_1.479_SS
l1Pt_10_l1Eta_1.479_l2Pt_20_l2Eta_1.479_OS
l1Pt_10_l1Eta_1.479_l2Pt_50_l2Eta_1.479_SS
l1Pt_10_l1Eta_1.479_l2Pt_50_l2Eta_1.479_OS
l1Pt_10_l1Eta_1.479_l2Pt_100_l2Eta_1.479_SS
l1Pt_10_l1Eta_1.479_l2Pt_100_l2Eta_1.479_OS
l1Pt_10_l1Eta_1.479_l2Pt_200_l2Eta_1.479_SS
l1Pt_10_l1Eta_1.479_l2Pt_200_l2Eta_1.479_OS
l1Pt_20_l1Eta_0_l2Pt_20_l2Eta_0_SS
l1Pt_20_l1Eta_0_l2Pt_20_l2Eta_0_OS
l1Pt_20_l1Eta_0_l2Pt_20_l2Eta_0.8_SS
l1Pt_20_l1Eta_0_l2Pt_20_l2Eta_0.8_OS
l1Pt_20_l1Eta_0_l2Pt_20_l2Eta_1.479_SS
l1Pt_20_l1Eta_0_l2Pt_20_l2Eta_1.479_OS
l1Pt_20_l1Eta_0_l2Pt_50_l2Eta_0_SS
l1Pt_20_l1Eta_0_l2Pt_50_l2Eta_0_OS
l1Pt_20_l1Eta_0_l2Pt_50_l2Eta_0.8_SS
l1Pt_20_l1Eta_0_l2Pt_50_l2Eta_0.8_OS
l1Pt_20_l1Eta_0_l2Pt_50_l2Eta_1.479_SS
l1Pt_20_l1Eta_0_l2Pt_50_l2Eta_1.479_OS
l1Pt_20_l1Eta_0_l2Pt_100_l2Eta_0_SS
l1Pt_20_l1Eta_0_l2Pt_100_l2Eta_0_OS
l1Pt_20_l1Eta_0_l2Pt_100_l2Eta_0.8_SS
l1Pt_20_l1Eta_0_l2Pt_100_l2Eta_0.8_OS
l1Pt_20_l1Eta_0_l2Pt_100_l2Eta_1.479_SS
l1Pt_20_l1Eta_0_l2Pt_100_l2Eta_1.479_OS
l1Pt_20_l1Eta_0_l2Pt_200_l2Eta_0_SS
l1Pt_20_l1Eta_0_l2Pt_200_l2Eta_0_OS
l1Pt_20_l1Eta_0_l2Pt_200_l2Eta_0.8_SS
l1Pt_20_l1Eta_0_l2Pt_200_l2Eta_0.8_OS
l1Pt_20_l1Eta_0_l2Pt_200_l2Eta_1.479_SS
l1Pt_20_l1Eta_0_l2Pt_200_l2Eta_1.479_OS
l1Pt_20_l1Eta_0.8_l2Pt_20_l2Eta_0.8_SS
l1Pt_20_l1Eta_0.8_l2Pt_20_l2Eta_0.8_OS
l1Pt_20_l1Eta_0.8_l2Pt_20_l2Eta_1.479_SS
l1Pt_20_l1Eta_0.8_l2Pt_20_l2Eta_1.479_OS
l1Pt_20_l1Eta_0.8_l2Pt_50_l2Eta_0.8_SS
l1Pt_20_l1Eta_0.8_l2Pt_50_l2Eta_0.8_OS
l1Pt_20_l1Eta_0.8_l2Pt_50_l2Eta_1.479_SS
l1Pt_20_l1Eta_0.8_l2Pt_50_l2Eta_1.479_OS
l1Pt_20_l1Eta_0.8_l2Pt_100_l2Eta_0.8_SS
l1Pt_20_l1Eta_0.8_l2Pt_100_l2Eta_0.8_OS
l1Pt_20_l1Eta_0.8_l2Pt_100_l2Eta_1.479_SS
l1Pt_20_l1Eta_0.8_l2Pt_100_l2Eta_1.479_OS
l1Pt_20_l1Eta_0.8_l2Pt_200_l2Eta_0.8_SS
l1Pt_20_l1Eta_0.8_l2Pt_200_l2Eta_0.8_OS
l1Pt_20_l1Eta_0.8_l2Pt_200_l2Eta_1.479_SS
l1Pt_20_l1Eta_0.8_l2Pt_200_l2Eta_1.479_OS
l1Pt_20_l1Eta_1.479_l2Pt_20_l2Eta_1.479_SS
l1Pt_20_l1Eta_1.479_l2Pt_20_l2Eta_1.479_OS
l1Pt_20_l1Eta_1.479_l2Pt_50_l2Eta_1.479_SS
l1Pt_20_l1Eta_1.479_l2Pt_50_l2Eta_1.479_OS
l1Pt_20_l1Eta_1.479_l2Pt_100_l2Eta_1.479_SS
l1Pt_20_l1Eta_1.479_l2Pt_100_l2Eta_1.479_OS
l1Pt_20_l1Eta_1.479_l2Pt_200_l2Eta_1.479_SS
l1Pt_20_l1Eta_1.479_l2Pt_200_l2Eta_1.479_OS
l1Pt_50_l1Eta_0_l2Pt_50_l2Eta_0_SS
l1Pt_50_l1Eta_0_l2Pt_50_l2Eta_0_OS
l1Pt_50_l1Eta_0_l2Pt_50_l2Eta_0.8_SS
l1Pt_50_l1Eta_0_l2Pt_50_l2Eta_0.8_OS
l1Pt_50_l1Eta_0_l2Pt_50_l2Eta_1.479_SS
l1Pt_50_l1Eta_0_l2Pt_50_l2Eta_1.479_OS
l1Pt_50_l1Eta_0_l2Pt_100_l2Eta_0_SS
l1Pt_50_l1Eta_0_l2Pt_100_l2Eta_0_OS
l1Pt_50_l1Eta_0_l2Pt_100_l2Eta_0.8_SS
l1Pt_50_l1Eta_0_l2Pt_100_l2Eta_0.8_OS
l1Pt_50_l1Eta_0_l2Pt_100_l2Eta_1.479_SS
l1Pt_50_l1Eta_0_l2Pt_100_l2Eta_1.479_OS
l1Pt_50_l1Eta_0_l2Pt_200_l2Eta_0_SS
l1Pt_50_l1Eta_0_l2Pt_200_l2Eta_0_OS
l1Pt_50_l1Eta_0_l2Pt_200_l2Eta_0.8_SS
l1Pt_50_l1Eta_0_l2Pt_200_l2Eta_0.8_OS
l1Pt_50_l1Eta_0_l2Pt_200_l2Eta_1.479_SS
l1Pt_50_l1Eta_0_l2Pt_200_l2Eta_1.479_OS
l1Pt_50_l1Eta_0.8_l2Pt_50_l2Eta_0.8_SS
l1Pt_50_l1Eta_0.8_l2Pt_50_l2Eta_0.8_OS
l1Pt_50_l1Eta_0.8_l2Pt_50_l2Eta_1.479_SS
l1Pt_50_l1Eta_0.8_l2Pt_50_l2Eta_1.479_OS
l1Pt_50_l1Eta_0.8_l2Pt_100_l2Eta_0.8_SS
l1Pt_50_l1Eta_0.8_l2Pt_100_l2Eta_0.8_OS
l1Pt_50_l1Eta_0.8_l2Pt_100_l2Eta_1.479_SS
l1Pt_50_l1Eta_0.8_l2Pt_100_l2Eta_1.479_OS
l1Pt_50_l1Eta_0.8_l2Pt_200_l2Eta_0.8_SS
l1Pt_50_l1Eta_0.8_l2Pt_200_l2Eta_0.8_OS
l1Pt_50_l1Eta_0.8_l2Pt_200_l2Eta_1.479_SS
l1Pt_50_l1Eta_0.8_l2Pt_200_l2Eta_1.479_OS
l1Pt_50_l1Eta_1.479_l2Pt_50_l2Eta_1.479_SS
l1Pt_50_l1Eta_1.479_l2Pt_50_l2Eta_1.479_OS
l1Pt_50_l1Eta_1.479_l2Pt_100_l2Eta_1.479_SS
l1Pt_50_l1Eta_1.479_l2Pt_100_l2Eta_1.479_OS
l1Pt_50_l1Eta_1.479_l2Pt_200_l2Eta_1.479_SS
l1Pt_50_l1Eta_1.479_l2Pt_200_l2Eta_1.479_OS
l1Pt_100_l1Eta_0_l2Pt_100_l2Eta_0_SS
l1Pt_100_l1Eta_0_l2Pt_100_l2Eta_0_OS
l1Pt_100_l1Eta_0_l2Pt_100_l2Eta_0.8_SS
l1Pt_100_l1Eta_0_l2Pt_100_l2Eta_0.8_OS
l1Pt_100_l1Eta_0_l2Pt_100_l2Eta_1.479_SS
l1Pt_100_l1Eta_0_l2Pt_100_l2Eta_1.479_OS
l1Pt_100_l1Eta_0_l2Pt_200_l2Eta_0_SS
l1Pt_100_l1Eta_0_l2Pt_200_l2Eta_0_OS
l1Pt_100_l1Eta_0_l2Pt_200_l2Eta_0.8_SS
l1Pt_100_l1Eta_0_l2Pt_200_l2Eta_0.8_OS
l1Pt_100_l1Eta_0_l2Pt_200_l2Eta_1.479_SS
l1Pt_100_l1Eta_0_l2Pt_200_l2Eta_1.479_OS
l1Pt_100_l1Eta_0.8_l2Pt_100_l2Eta_0.8_SS
l1Pt_100_l1Eta_0.8_l2Pt_100_l2Eta_0.8_OS
l1Pt_100_l1Eta_0.8_l2Pt_100_l2Eta_1.479_SS
l1Pt_100_l1Eta_0.8_l2Pt_100_l2Eta_1.479_OS
l1Pt_100_l1Eta_0.8_l2Pt_200_l2Eta_0.8_SS
l1Pt_100_l1Eta_0.8_l2Pt_200_l2Eta_0.8_OS
l1Pt_100_l1Eta_0.8_l2Pt_200_l2Eta_1.479_SS
l1Pt_100_l1Eta_0.8_l2Pt_200_l2Eta_1.479_OS
l1Pt_100_l1Eta_1.479_l2Pt_100_l2Eta_1.479_SS
l1Pt_100_l1Eta_1.479_l2Pt_100_l2Eta_1.479_OS
l1Pt_100_l1Eta_1.479_l2Pt_200_l2Eta_1.479_SS
l1Pt_100_l1Eta_1.479_l2Pt_200_l2Eta_1.479_OS
l1Pt_200_l1Eta_0_l2Pt_200_l2Eta_0_SS
l1Pt_200_l1Eta_0_l2Pt_200_l2Eta_0_OS
l1Pt_200_l1Eta_0_l2Pt_200_l2Eta_0.8_SS
l1Pt_200_l1Eta_0_l2Pt_200_l2Eta_0.8_OS
l1Pt_200_l1Eta_0_l2Pt_200_l2Eta_1.479_SS
l1Pt_200_l1Eta_0_l2Pt_200_l2Eta_1.479_OS
l1Pt_200_l1Eta_0.8_l2Pt_200_l2Eta_0.8_SS
l1Pt_200_l1Eta_0.8_l2Pt_200_l2Eta_0.8_OS
l1Pt_200_l1Eta_0.8_l2Pt_200_l2Eta_1.479_SS
l1Pt_200_l1Eta_0.8_l2Pt_200_l2Eta_1.479_OS
l1Pt_200_l1Eta_1.479_l2Pt_200_l2Eta_1.479_SS
l1Pt_200_l1Eta_1.479_l2Pt_200_l2Eta_1.479_OS
)

for tree in ${TREES[@]}; do
./main -f DataEGBCDFlipTrees.root -m MCSamplesFlipTrees.root -a 1 -s $tree -n dbSingle.db

done
