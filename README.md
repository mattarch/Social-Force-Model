# asl_project
Social Force Model for Pedestrian Dynamics

Google Docs Link:
https://docs.google.com/document/d/1Pd_naY5nFvGuTWZ7cf6Ku0H6xgkG0ytWERB9O8ZGhYY/edit?usp=sharing

Presentation Link for first meeting with supervisor:
https://docs.google.com/presentation/d/1ZmCgMJ247X3icRssxIMaQTaob5e6haUbSb5RQ8D-COo/edit?usp=sharing

## compilation
`gcc social_force.c parse_args.c utility.c testing.c social_force_model_basic.c -lm -O3 -ffast-math -g`

## Debug mode
To debug the application compile the file with the flag `-DDEBUG`.
The output of the debug mode is saved in the folder "test".

Visualize the output in test with:
`python3 visualization_basic.py <filename>`

## Testing mode
To run tests add the flag `--test` to the executable.
To add tests first add parameters in the `test_set.h` file and add a new call inside the `add_tests()` function inside `testing.h`.
To add implementations to check add a new call inside the `add_implementation()` function inside `testing.h`.



