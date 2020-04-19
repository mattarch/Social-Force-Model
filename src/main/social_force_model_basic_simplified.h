/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef SIMPLIFIED_H_ /* Include guard */
#define SIMPLIFIED_H_

extern char fiename_global[80];

void initialize_people_simplified(double *position, double *desired_direction, double *final_destination, double *desired_speed, int n);
void initialize_borders_simplified(double *borders, int n_borders);
void update_desired_direction_simplified(double *position, double *final_destination, double *desired_direction, int n);
void update_acceleration_term_simplified(double *desired_direction, double *acceleration_term, double *actual_velocity, double *desired_speed, int n);
//void compute_actual_velocity_simplified(double *actual_speed, double *desired_direction, double *actual_velocity, int n);
void update_people_repulsion_term_simplified(double *position, double *desired_direction, double *actual_speed, double *Repulsion_term, int n);
void update_border_repulsion_term_simplified(double *position, double *borders, double *border_repulsion_term, int n, int n_borders);
void compute_social_force_simplified(double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term, double *social_force, int n, int n_borders);
void update_position_simplified(double *position, double *desired_direction, double *actual_speed, double *social_force, double *actual_velocity, double *desired_speed, int n);
void simulation_basic_simplified(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination, double *borders, double *actual_velocity, double *acceleration_term,
                      double *people_repulsion_term, double *border_repulsion_term, double *social_force, double *desired_speed);
void test_simulation_basic_simplified(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination, double *borders, double *actual_velocity, double *acceleration_term,
                           double *people_repulsion_term, double *border_repulsion_term, double *social_force, double *desired_speed);

#endif