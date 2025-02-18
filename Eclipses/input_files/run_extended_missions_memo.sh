#!/bin/bash

export n_trial=20 # number of trials (20 good for reasonable statistics)
export prototype_mode=0 # 0=full sim. 1= 1 tile, 2=10 tiles, 3=~290 tiles (1/10th of sky)

export date_sub="160618_200k_pm$prototype_mode" # substring to put before output files

export ps_pri_file="../../preProcessing/sTIC-selection/shemi_nhemi_4Morb/shemi_nhemi_200k.sav" # 200k over 2yr
export fcam_coord_pri="../cameraPointings/shemi_nhemi_orbits.dat"

#here is where you write all your extended mission names, and their tails (since idk how to bash)
# 1yr runs
e0="shemi_nhemi_nhemi"
e1="shemi_nhemi_npole"
e2="shemi_nhemi_shemiAvoid"
e3="shemi_nhemi_elong"
e4="shemi_nhemi_eshort"
e5="shemi_nhemi_hemis14d"

t0="nhemi"
t1="npole"
t2="shemiAvoid"
t3="elong" # n.b. this is only for the camera pointings call
t4="eshort"
t5="hemis14d"

# declare array variables
declare -a arr=("$e0" "$e1" "$e2" "$e3" "$e4" "$e5")
declare -a tails=("$t0" "$t1" "$t2" "$t3" "$t4" "$t5")

# loop thru above array
ind=0
for i in "${arr[@]}"
do
  echo "Running $i ..."
  export ext_mission_name="$i"

  export fcam_coord_ext="../cameraPointings/"${tails[$ind]}"_orbits.dat"
  export ps_ext_file="../../preProcessing/sTIC-selection/"$i"_4Morb/"$i"_200k.sav"
  idl -e main_ext &> "../output_files/"$date_sub"_"$i"_t"$n_trial".out" &
  let "ind+=1"

done
