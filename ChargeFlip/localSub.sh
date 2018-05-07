#!/bin/bash

func() {
    root <<EOF
//cout<<"pouet "<<$1<<"  "<<$2<<endl;
.L ZChargeFlip.C+
doFits(0,$1,$2);
.q
EOF
#sleep 3s
}



func 0 48 &
func 48 96 &
func 96 144 &
func 144 192 &
func 192 240 &
func 240 288 &

