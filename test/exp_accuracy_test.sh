#!/bin/bash

# set variables
options="stdc_float stdc_double SIMD_float SIMD_double"
os=$(uname)
filesuffix=".txt"
prog="exp_accuracy_test"

# compile program
# if macOS -> use gcc-9, else -> gcc
if [ $os == "Darwin" ]; then
  gcc-9 *.c -lm -O3 -mavx2 -mfma -ffast-math -g -o $prog
else
  gcc *.c -lm -O3 -mavx2 -mfma -ffast-math -g -o $prog
fi

# perform benchmark or test
select opt in $options; do
  if [ $opt == "stdc_float" ]; then
    echo "Perform stdc_float: "
    ./$prog 0 > "$opt$filesuffix"
  elif [ $opt == "stdc_double" ]; then
    echo "Perform stdc_double: "
    ./$prog 1 > "$opt$filesuffix"
  elif [ $opt == "SIMD_float" ]; then
    echo "Perform SIMD_float: "
    ./$prog 2 > "$opt$filesuffix"
  elif [ $opt == "SIMD_double" ]; then
    echo "Perform SIMD_double: "
    ./$prog 3 > "$opt$filesuffix"
  else
    echo "Error: chose 1, 2, 3 or 4"
    exit
  fi
  # open results
  if [ $os == "Darwin" ] || [ $os == "Linux" ]; then
    python3 plot-errors.py "$opt$filesuffix"
  else
    python plot-errors.py "$opt$filesuffix"
  fi
  exit
done

