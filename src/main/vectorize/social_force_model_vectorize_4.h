/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef VECTORIZE_4_H_ /* Include guard */
#define VECTORIZE_4_H_

#include <immintrin.h>

void update_desired_direction_vectorize_4(float *position, float *final_destination, float *desired_direction, int n);
void update_acceleration_term_vectorize_4(float *desired_direction, float *acceleration_term, float *actual_velocity, float *desired_speed, int n);
//void compute_actual_velocity_vectorize_4(float *actual_speed, float *desired_direction, float *actual_velocity, int n);
void update_people_repulsion_term_vectorize_4(float *position, float *desired_direction, float *actual_speed, float *Repulsion_term, int n);
void update_border_repulsion_term_vectorize_4(float *position, float *borders, float *border_repulsion_term, int n, int n_borders);
void compute_social_force_vectorize_4(float *acceleration_term, float *people_repulsion_term, float *border_repulsion_term, float *social_force, int n, int n_borders);
void update_position_vectorize_4(float *position, float *desired_direction, float *actual_speed, float *social_force, float *actual_velocity, float *desired_max_speed, int n);
void simulation_basic_vectorize_4(int number_of_people, int n_timesteps, float *position, float *speed, float *desired_direction, float *final_destination, float *borders, float *actual_velocity, float *acceleration_term,
                                  float *people_repulsion_term, float *border_repulsion_term, float *social_force, float *desired_speed, float *desired_max_speed);

#endif