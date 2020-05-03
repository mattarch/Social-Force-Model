/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef VECTORIZE_2_5_1_H_ /* Include guard */
#define VECTORIZE_2_5_1_H_

void simulation_basic_vectorize_2_5_1(int number_of_people, int n_timesteps, float *position, float *speed, float *desired_direction, float *final_destination,
                                 float *borders, float *actual_velocity, float *acceleration_term, float *people_repulsion_term, float *border_repulsion_term,
                                 float *social_force, float *desired_speed, float *desired_max_speed);

#endif