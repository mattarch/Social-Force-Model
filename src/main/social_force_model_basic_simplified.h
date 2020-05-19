/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef SIMPLIFIED_H_ /* Include guard */
#define SIMPLIFIED_H_

//void initialize_max_speed_float(float *desired_speed, float *desired_max_speed, int n);
void update_desired_direction_simplified(float *position, float *final_destination, float *desired_direction, int n);
void update_acceleration_term_simplified(float *desired_direction, float *acceleration_term, float *actual_velocity, float *desired_speed, int n);
//void compute_actual_velocity_simplified(float *actual_speed, float *desired_direction, float *actual_velocity, int n);
void update_people_repulsion_term_simplified(float *position, float *desired_direction, float *actual_speed, float *Repulsion_term, int n);
void update_border_repulsion_term_simplified(float *position, float *borders, float *border_repulsion_term, int n, int n_borders);
void compute_social_force_simplified(float *acceleration_term, float *people_repulsion_term, float *border_repulsion_term, float *social_force, int n, int n_borders);
void update_position_simplified(float *position, float *desired_direction, float *actual_speed, float *social_force, float *actual_velocity, float *desired_max_speed, int n);
void simulation_basic_simplified(int number_of_people, int n_timesteps, float *position, float *speed, float *desired_direction, float *final_destination, float *borders, float *actual_velocity, float *acceleration_term,
                      float *people_repulsion_term, float *border_repulsion_term, float *social_force, float *desired_speed, float* desired_max_speed);

void test_simulation_basic_simplified(int number_of_people, int n_timesteps, float *position, float *speed, float *desired_direction, float *final_destination, float *borders, float *actual_velocity, float *acceleration_term,
                           float *people_repulsion_term, float *border_repulsion_term, float *social_force, float *desired_speed, float* desired_max_speed);

#endif