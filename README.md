# ASL_project
This is the GitLap repository of group 48.
Project topic: Social Force Model for Pedestrian Dynamics

Goal was to implement a fast social force model. The final report includes an analysis and results of the work.

This repository contains all C and Python code that was used to generate and analyze the data. The optics of the plots were changed specifically for the report. Further, a random element in the code can change the visualization even when all parameters are identical.



%Google Docs Link:
%https://docs.google.com/document/d/1Pd_naY5nFvGuTWZ7cf6Ku0H6xgkG0ytWERB9O8ZGhYY/edit?usp=sharing

%Presentation Link for first meeting with supervisor:
%https://docs.google.com/presentation/d/1ZmCgMJ247X3icRssxIMaQTaob5e6haUbSb5RQ8D-COo/edit?usp=sharing

## Compilation
Compile the code with running the bash script 'run_simulation.sh'. After compilation was succesfull the user has 3 different option:
1. Benchmark
2. Test
3. Visualization and Test
`gcc social_force.c parse_args.c utility.c testing.c social_force_model_basic.c -lm -O3 -ffast-math -g`

## Benchmark
This option benchmarks chosen versions for predifined input sizes. Results of the benchmarking process are constantly prited to the command line and also saved in a text file under the folder 'benchmark'. A python script plot the results after the benchmarking is completed for all versions and input sizes.

##Test
The test option allows the user to run testcases and finite-difference tests for the baseline version. Additional versions will be compared to the values computed by this baseline implementation. The comparison takes place on a one-timestep difference, where the states of every version is adjusted after one timestep to match the basic version.

##Visualization
Next to running the same tests as option "Test", this option also outputs a visualization of the basic implementation. A text file, saved under the folder 'test', contains all the values computed during the simulation. Later a python script will visualize the results.


To debug the application compile the file with the flag `-DDEBUG`.
The output of the debug mode is saved in the folder "test".

Visualize the output in test with:
`python3 visualization_basic.py <filename>`

## Testing mode
To run tests add the flag `--test` to the executable.
To add tests first add parameters in the `test_set.h` file and add a new call inside the `add_tests()` function inside `testing.h`.
To add implementations to check add a new call inside the `add_implementation()` function inside `testing.h`.



