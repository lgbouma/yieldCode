#!/bin/bash

export n_trial=1 # number of trials (20 good for reasonable statistics)
export prototype_mode=1 # 0=full sim. 1= 1 tile, 2=10 tiles, 3=~290 tiles (1/10th of sky)

export date_sub="160526_pm$prototype_mode" # substring to put before output files

export ps_pri_file="../../preProcessing/sTIC-selection/shemi_nhemi/shemi_nhemi.sav" # 200k over 2yr
export fcam_coord_pri="../cameraPointings/shemi_nhemi_coord.dat"

#here is where you write all your extended mission names, and their tails (since idk how to bash)
# 1yr runs
e0="shemi_nhemi_nhemi"
e1="shemi_nhemi_npole"
e2="shemi_nhemi_ecliptic"
# 2yr runs
#e3="shemi_nhemi_ecliptic_ecliptic_100perExtYr"
#e4="shemi_nhemi_ecliptic_nhemi_100perExtYr"
#e5="shemi_nhemi_ecliptic_npole_100perExtYr"
#e6="shemi_nhemi_shemi_nhemi"
#e7="shemi_nhemi_npole_spole"
#e8="shemi_nhemi_shemi_shemi_100perExtYr"
#e9="shemi_nhemi_npole_npole_100perExtYr"

t0="nhemi"
t1="npole"
t2="ecliptic"
#t3="ecliptic_ecliptic" # n.b. this is only for the camera pointings call
#t4="ecliptic_nhemi"
#t5="ecliptic_npole"
#t6="shemi_nhemi"
#t7="npole_spole"
#t8="shemi_shemi"
#t9="npole_npole"

# declare array variables
declare -a arr=("$e0" "$e1" "$e2") # "$e3" "$e4" "$e5" "$e6" "$e7" "$e8" "$e9")
declare -a tails=("$t0" "$t1" "$t2") # "$t3" "$t4" "$t5" "$t6" "$t7" "$t8" "$t9")

# loop thru above array
ind=0
for i in "${arr[@]}"
do
  echo "Running $i ..."
  export ext_mission_name="$i"

  export fcam_coord_ext="../cameraPointings/"${tails[$ind]}"_coord.dat"
  export ps_ext_file="../../preProcessing/sTIC-selection/"$i"/"$i".sav"
  idl -e main_ext &> "../output_files/"$date_sub"_"$i"_t"$n_trial".raw" &
  let "ind+=1"

done
