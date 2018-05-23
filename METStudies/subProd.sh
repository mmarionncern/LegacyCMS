#!/bin/bash


METTypes=( "pat_patType1PhiCorrectedPFMet" "pat_patPFMetMVAPhi" "Calo_caloType1PhiCorrectedMet" )
Comp=( "para" "perp" )


for im in ${METTypes[@]}; do
    echo $im
    for ic in  ${Comp[@]}; do

    bsub -q 1nd source subTask.sh $im $im $ic 0 1 0 1 0
    bsub -q 1nd source subTask.sh $im $im $ic 1 1 0 1 0
    #bsub -q 8nh source subTask.sh $im $im $ic 1 1 1 1 0
    #bsub -q 8nh source subTask.sh $im $im $ic 0 1 1 1 1
    bsub -q 1nd source subTask.sh $im $im $ic 0 1 0 1 1
     
    done
done
