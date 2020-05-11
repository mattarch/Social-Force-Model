/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef OPTV_2_5_1_H_ /* Include guard */
#define OPTV_2_5_1_H_

void simulation_basic_optv_2_5_1(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination,
                                 double *borders, double *actual_velocity, double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term,
                                 double *social_force, double *desired_speed, double *desired_max_speed);

#endif