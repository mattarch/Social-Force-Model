#!/bin/bash

# set variables
options="Benchmark Test"
timesteps=(300)
persons=(8)
os=$(uname)
currentdate=$(date +"%Y-%m-%d")
currenttime=$(date +"%H-%M-%S")
underline="_"
filesuffix=".txt"
timelabel="$underline$currentdate$underline$currenttime$filesuffix"
prog="social_force"
firstarg = $1

# compile program
# if macOS -> use gcc-9, else -> gcc
if [ $os == "Darwin" ]; then
  gcc-9 social_force.c parse_args.c utility.c testing.c social_force_model_basic.c -lm -O3 -ffast-math -g -o $prog
else
  gcc social_force.c parse_args.c utility.c testing.c social_force_model_basic.c -lm -O3 -ffast-math -g -o $prog
fi

# perform benchmark or test
select opt in $options; do
  if [ $opt == "Benchmark" ]; then
    echo "Perform Benchmarks: "
    
    for t in "${timesteps[@]}"
    do
      echo "-Execute Benchmarks for $t time steps:"

      # perform 1 benchmark to create files
			echo "--run benchmark for 4 persons"
      ./$prog --n_timesteps=$t --n_people=4 | while read line
      do
        set -- $line
        filename=$firstarg$timelabel
        echo $line
        echo $line >> ../../benchmark/$filename
      done

      for p in "${persons[@]}"
      do
        echo "--run benchmark for $p persons"
        ./$prog --n_timesteps=$t --n_people=$p | while read line
        do
          set -- $line
          filename=$firstarg$timelabel
          echo $line
          echo $line >> ../../benchmark/$filename
        done
      done
      python3 ../../benchmark/plot-benchmark.py ../../benchmark/$firstarg$timelabel
    done
    exit
  elif [ $opt == "Test" ]; then
    echo "Perform Tests:"
    ./$prog --test
    exit
  else
    echo "Error: chose 1 or 2"
    exit
  fi
done

