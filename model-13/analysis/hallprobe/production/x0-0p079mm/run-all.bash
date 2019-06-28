#!/usr/bin/env bash

func=$1
subfolder=$2
side=$3

folder=/home/imas/repos/si-dipoles-bc/model-13/analysis/hallprobe/production/x0-0p079mm

function f1 {
  subfolder=$1
  side=$2
  cd $folder/BC-02/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-03/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-04/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
}

function f2 {
  subfolder=$1
  side=$2
  cd $folder/BC-05/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-06/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-07/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
}

function f3 {
  subfolder=$1
  side=$2
  cd $folder/BC-08/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-09/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-10/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
}

function f4 {
  subfolder=$1
  side=$2
  cd $folder/BC-11/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-12/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-13/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
}

function f5 {
  subfolder=$1
  side=$2
  cd $folder/BC-14/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-15/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-16/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
}

function f6 {
  subfolder=$1
  side=$2
  cd $folder/BC-17/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-18/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-20/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
  cd $folder/BC-21/$subfolder/$side; fac-fma-analysis.py run; cd ../../../../
}

# function f1 {
#   subfolder=$1
#   side=$2
#   cd $folder/BC-02/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-03/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-04/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
# }

# function f2 {
#   subfolder=$1
#   side=$2
#   cd $folder/BC-05/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-06/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-07/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
# }

# function f3 {
#   subfolder=$1
#   side=$2
#   cd $folder/BC-08/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-09/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-10/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
# }

# function f4 {
#   subfolder=$1
#   side=$2
#   cd $folder/BC-11/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-12/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-13/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
# }

# function f5 {
#   subfolder=$1
#   side=$2
#   cd $folder/BC-14/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-15/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-16/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
# }

# function f6 {
#   subfolder=$1
#   side=$2
#   cd $folder/BC-17/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-18/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-20/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
#   cd $folder/BC-21/$subfolder/$side; fac-fma-analysis.py rawfield; fac-fma-analysis.py trajectory; cd ../../../../
# }

$func $subfolder $side
