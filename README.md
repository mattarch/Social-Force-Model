# Advanced Systems Lab - Spring 2020
This is the GitLab repository of group 48.
Project topic: Social Force Model for Pedestrian Dynamics

Goal was to implement a fast social force model. The final report includes an analysis and results of the work.

This repository contains all C and Python code that was used to generate and analyze the data. The optics of the plots were changed specifically for the report. Further, a random element in the code can change the visualization even when all parameters are identical.


## Compilation
Compile the code with running the bash script 'src/main/run_simulation.sh'. After compilation was succesfull the user has 3 different option:
1. Benchmark
2. Test
3. Visualization

The flags used are:
`gcc -lm -O3 -mavx2 -mfma -ffast-math -g`

## Benchmark
This option benchmarks predefined for predifined input sizes. Results of the benchmarking process are constantly printed to the command line and also saved in a text file under the folder 'benchmark'. A python script plot the results after the benchmarking is completed for all versions and input sizes.

## Test
The test option allows the user to run testcases and finite-difference tests for the baseline version. Additional versions will be compared to the values computed by this baseline implementation. The comparison takes place on a one-timestep difference, where the states of every version is adjusted after one timestep to match the basic version.

## Visualization
This option outputs a visualization of the basic implementation. A text file, saved under the folder 'test', contains all the values computed during the simulation. Later a python script will visualize the results.

If the user wants to visualize a previously run simulation, then the output can be visualized in test with:
`python3 visualization_basic.py <filename>`




