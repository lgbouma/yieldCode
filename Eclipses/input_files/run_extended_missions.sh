#!/bin/bash
# interpret this via bash, always

export n_trial=5 # number of trials (10 good for reasonable statistics)
export prototype_mode=3 # 0=full sim. 1= 1 tile, 2=10 tiles, 3=~290 tiles (1/10th of sky)

export date_sub="160526" # substring to put before output files

export ps_pri_file="../../preProcessing/sTIC-selection/shemi_nhemi/shemi_nhemi.sav" # 200k over 2yr
export fcam_coord_pri="../cameraPointings/shemi_nhemi_coord.dat"

#here is where you write all your extended mission names, and their tails (since idk how to bash)
e0="shemi_nhemi_nhemi"
e1="shemi_nhemi_npole"
t0="nhemi"
t1="npole"

# declare array variables
declare -a arr=("$e0" "$e1")
declare -a tails=("$t0" "$t1")

# loop thru above array
ind=0
for i in "${arr[@]}"
do
  echo "Running $i ..."
  export ext_mission_name="$i"

  export fcam_coord_ext="../cameraPointings/${tails[$ind]}_coord.dat"
  export ps_ext_file="../../preProcessing/sTIC-selection/$i/$i.sav"
  idl -e main_ext &> "../output_files/$i.raw" &
  let "ind+=1"

done
