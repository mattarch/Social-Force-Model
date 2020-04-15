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
void output_to_file_initial_state(char *filename, double *position, double *actual_speed, double *desired_direction, double *final_destination, int n, int n_timestep);
void output_to_file_persons(char *filename, double *position, double *actual_speed, double *desired_direction, double *final_destination, int n, int n_timestep);
void output_to_file_constants(char *);
void free_all(int n, ...);

double sampleNormal(double sigma, double mu);

#endif