/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef BASIC_H_ /* Include guard */
#define BASIC_H_

extern char filename_global[80];

void initialize_people(double *position, double *desired_direction, double *final_destination, double *desired_speed, int n);
void initialize_borders(double *borders, int n_borders);
void update_desired_direction(double *position, double *final_destination, double *desired_direction, int n);
void update_acceleration_term(double *desired_direction, double *acceleration_term, double *actual_velocity, double *desired_speed, int n);
void compute_actual_velocity(double *actual_speed, double *desired_direction, double *actual_velocity, int n);
void update_people_repulsion_term(double *position, double *desired_direction, double *actual_speed, double *Repulsion_term, int n);
void update_border_repulsion_term(double *position, double *borders, double *border_repulsion_term, int n, int n_borders);
void compute_social_force(double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term, double *social_force, int n, int n_borders);
void update_position(double *position, double *desired_direction, double *actual_speed, double *social_force, double *actual_velocity, double *desired_speed, int n);
void simulation_basic(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination, double *borders, double *actual_velocity, double *acceleration_term,
                      double *people_repulsion_term, double *border_repulsion_term, double *social_force, double *desired_speed, double *desired_max_speed);
void test_simulation_basic(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination, double *borders, double *actual_velocity, double *acceleration_term,
                           double *people_repulsion_term, double *border_repulsion_term, double *social_force, double *desired_speed, double *desired_max_speed);

#endif