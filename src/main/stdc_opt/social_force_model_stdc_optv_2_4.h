/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef OPTV_3_2_H_ /* Include guard */
#define OPTV_3_2_H_

extern char fiename_global[80];
void simulation_basic_optv_2_4(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination, double *borders, double *actual_velocity, double *position_temp,
                               double *speed_temp, double *desired_direction_temp, double *social_force, double *desired_speed, double *desired_max_speed);

#endif