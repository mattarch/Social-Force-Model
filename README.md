# Advanced Systems Lab - Spring 2020
This is the GitLab repository of group 48.
Project topic: Social Force Model for Pedestrian Dynamics

The goal was to implement a fast social force model. The final report includes an analysis and results of the work.

This repository contains all C and Python code that was used to generate and analyze the data.


## Compilation
Compile the code with running the bash script 'src/main/run_simulation.sh'. After compilation was succesfull the user has 3 different option:
1. Benchmark
2. Test
3. Visualization

The flags used are:
`gcc -lm -O3 -mavx2 -mfma -ffast-math -g`

## Benchmark
This option runs the benchmark for the implementations discussed in the report for predifined input sizes. Results of the benchmarking process are constantly printed to the command line and also saved in a text file under the folder 'benchmark'. A python script plots the results after the benchmarking is completed for all versions and input sizes.

## Test
The test option allows the user to run testcases and finite-difference tests for the straightforward version. Additional versions will be compared to the values computed by this straightforward implementation. The comparison takes place after eacht time step.

## Visualization
This option outputs a visualization of the straightforward implementation. A text file, saved under the folder 'test', contains all the values computed during the simulation. Later a python script will visualize the results.

A previously run simulation can be visualized by running the visualization script in the folder 'test' with:
`python3 visualization_basic.py <filename>` 




