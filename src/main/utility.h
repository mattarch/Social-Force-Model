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
void get_filename();
void output_to_file_initial_state(char *filename, float *position, float *actual_speed, float *desired_direction, float *final_destination, int n, int n_timestep);
void output_to_file_persons(char *filename, float *position, float *actual_speed, float *desired_direction, float *final_destination, int n, int n_timestep);
void output_to_file_constants(char *);
void free_all(int n, ...);
float exp_fast(float x);
float exp_taylor(float x);

float sampleNormal(float sigma, float mu);
void set_zero(float *p, int size);

#endif