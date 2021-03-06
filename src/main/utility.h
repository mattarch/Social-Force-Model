/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef UTILITY_H_ /* Include guard */
#define UTILITY_H_

extern struct arguments arguments;
extern char filename_global[80];


/* function defined in the header file itself */
int contains_substring(char* haystack, char *needle);

void get_filename();
void output_to_file_initial_state(char *filename, float *position, float *actual_speed, float *desired_direction, float *final_destination, int n, int n_timestep);
void output_to_file_persons(char *filename, float *position, float *actual_speed, float *desired_direction, float *final_destination, int n, int n_timestep);
void output_to_file_constants(char *);
void output_to_file_initial_state_double(char *filename, double *position, double *actual_speed, double *desired_direction, double *final_destination, int n, int n_timestep);
void output_to_file_persons_double(char *filename, double *position, double *actual_speed, double *desired_direction, double *final_destination, int n, int n_timestep);
void output_to_file_constants_double(char *);
void free_all(int n, ...);
float exp_fast_float(float x);
double exp_fast_double(double x);

float exp_taylor(float x);

float sampleNormal(float sigma, float mu);
void set_zero(float *p, int size);
void set_zero_double(double *p, int size);

#endif