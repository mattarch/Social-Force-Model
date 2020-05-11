/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#include <stdio.h>
#include <math.h>

#include "social_force_model_stdc_optv_3_2.h"
#include "../social_force.h"
#include "../utility.h"

extern char filename_global[80];
extern struct arguments arguments;

void simulation_basic_optv_3_2(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination, double *borders, double *actual_velocity, double *position_temp,
                               double *speed_temp, double *desired_direction_temp, double *social_force, double *desired_speed, double *desired_max_speed)
{
  // start simulation
  int n = number_of_people;
  const double inv_sigma = 1 / SIGMA; // 1 div -> 1 flop
  double border_0 = borders[0];
  double border_1 = borders[1];
  // simulate steps
  for (int step = 0; step < n_timesteps; step = step + 2)
  {
    for (int i = 0; i < n; i=i+2)
    {
      speed[i] *= TIMESTEP;
      speed[i+1] *= TIMESTEP;
    }
    // iterate over every person for step % 2 == 0
    for (int i = 0; i < n; i=i+2)
    {
      /************************************************/
      // load values
      /************************************************/
      double social_force_x = social_force[IndexX(i)];
      double social_force_y = social_force[IndexY(i,n)];
      double rx_a = position[IndexX(i)];
      double ry_a = position[IndexY(i,n)];
      double ex_a = desired_direction[IndexX(i)];
      double ey_a = desired_direction[IndexY(i,n)];
      double actual_velocity_x = actual_velocity[IndexX(i)];
      double actual_velocity_y = actual_velocity[IndexY(i,n)];
      double max_speed = desired_max_speed[i];
      double desired_speed_value = desired_speed[i];
      double target_x = final_destination[IndexX(i)];
      double target_y = final_destination[IndexY(i,n)];

      /************************************************/

      double social_force_x_1 = social_force[IndexX(i+1)];// 0.0;
      double social_force_y_1 = social_force[IndexY(i+1,n)];// 0.0;
      double rx_a_1 = position[IndexX(i+1)];
      double ry_a_1 = position[IndexY(i+1,n)];
      double ex_a_1 = desired_direction[IndexX(i+1)];
      double ey_a_1 = desired_direction[IndexY(i+1,n)];
      double actual_velocity_x_1 = actual_velocity[IndexX(i+1)];
      double actual_velocity_y_1 = actual_velocity[IndexY(i+1,n)];
      double max_speed_1 = desired_max_speed[(i+1)];
      double desired_speed_value_1 = desired_speed[(i+1)];
      double target_x_1 = final_destination[IndexX(i+1)];
      double target_y_1 = final_destination[IndexY(i+1,n)];

       /************************************************/
       // UPDATE BORDER REPULSION TERM
       /************************************************/
       //unrolled borders since there are only 2 in the scenario
      
        
      double ry_aB_0 = ry_a - border_0; //1 add
      double ry_aB_1 = ry_a - border_1; //1 add

      double r_aB_0_norm = ry_aB_0 > 0 ? ry_aB_0 : -ry_aB_0;
      double r_aB_1_norm = ry_aB_1 > 0 ? ry_aB_1 : -ry_aB_1;

      double shared_expression_0 = exp_fast((-r_aB_0_norm) * INV_R) * UTIMESR / r_aB_0_norm; //1 exp, 2 mult, 1 div
      double shared_expression_1 = exp_fast((-r_aB_1_norm) * INV_R) * UTIMESR / r_aB_1_norm; //1 exp, 2 mult, 1 div

      double repulsion_0 = shared_expression_0 * ry_aB_0; //1 mult
      double repulsion_1 = shared_expression_1 * ry_aB_1; //1 mult

      social_force_y += repulsion_0 + repulsion_1; //1 add
      
      /************************************************/

      double ry_a_1B_0 = ry_a_1 - border_0; //1 add
      double ry_a_1B_1 = ry_a_1 - border_1; //1 add

      double r_a_1B_0_norm = ry_a_1B_0 > 0 ? ry_a_1B_0 : -ry_a_1B_0;
      double r_a_1B_1_norm = ry_a_1B_1 > 0 ? ry_a_1B_1 : -ry_a_1B_1;

      double shared_expression_0_1 = exp_fast((-r_a_1B_0_norm) * INV_R) * UTIMESR / r_a_1B_0_norm; //1 exp, 2 mult, 1 div
      double shared_expression_1_1 = exp_fast((-r_a_1B_1_norm) * INV_R) * UTIMESR / r_a_1B_1_norm; //1 exp, 2 mult, 1 div

      double repulsion_0_1 = shared_expression_0_1 * ry_a_1B_0; //1 mult
      double repulsion_1_1 = shared_expression_1_1 * ry_a_1B_1; //1 mult

      social_force_y_1 += repulsion_0_1 + repulsion_1_1; //1 add

      /************************************************/
      // UPDATE PEOPLE REPULSION TERM
      /************************************************/

      // solve common repulsion of alpha0 and alpha1
      double delta_a = speed[i];
      double delta_a_1 = speed[i+1];

      double rx_a0a1 = rx_a - rx_a_1;
      double ry_a0a1 = ry_a - ry_a_1;
      double r_a0a1_norm = sqrt(rx_a0a1 * rx_a0a1 + ry_a0a1 * ry_a0a1);

      // from now on compute repulsion term directly for both alpha and beta
      //me stands for "minus e"
      double rx_a0a1_mex = rx_a0a1 - delta_a_1 * ex_a_1; //1 add, 1 mult
      double rx_a1a0_mex = -rx_a0a1 - delta_a * ex_a; //1 add, 1 mult

      double ry_a0a1_mey = ry_a0a1 - delta_a_1 * ey_a_1; //1 add, 1 mult
      double ry_a1a0_mey = -ry_a0a1 - delta_a * ey_a; //1 add, 1 mult

      double r_a0a1_me_norm = sqrt(rx_a0a1_mex * rx_a0a1_mex + ry_a0a1_mey * ry_a0a1_mey); //1 add, 2 mult, 1 sqrt
      double r_a1a0_me_norm = sqrt(rx_a1a0_mex * rx_a1a0_mex + ry_a1a0_mey * ry_a1a0_mey); //1 add, 2 mult, 1 sqrt

      double norm_sum_a = r_a0a1_norm + r_a0a1_me_norm;                                //1 add
      double norm_sum_a1 = r_a0a1_norm + r_a1a0_me_norm;                                //1 add

      double repulsion_x = rx_a0a1 / r_a0a1_norm + rx_a0a1_mex / r_a0a1_me_norm; //1 add, 2 div
      double repulsion_x_a1 = -rx_a0a1 / r_a0a1_norm + rx_a1a0_mex / r_a1a0_me_norm; //1 add, 2 div

      double repulsion_y = ry_a0a1 / r_a0a1_norm + ry_a0a1_mey / r_a0a1_me_norm; //1 add, 2 div
      double repulsion_y_a1 = -ry_a0a1 / r_a0a1_norm + ry_a1a0_mey / r_a1a0_me_norm; //1 add, 2 div

      double b_a = sqrt(norm_sum_a * norm_sum_a - delta_a_1 * delta_a_1) * 0.5; //1 add, 3 mult, 1 sqrt
      double b_a1 = sqrt(norm_sum_a1 * norm_sum_a1 - delta_a * delta_a) * 0.5; //1 add, 3 mult, 1 sqrt

      double common_factor_a = exp_fast(-b_a * inv_sigma) * norm_sum_a * DIV_FACTOR / b_a; //2 mult, 2 div, 1 exp
      double common_factor_a1 = exp_fast(-b_a1 * inv_sigma) * norm_sum_a1 * DIV_FACTOR / b_a1; //2 mult, 2 div, 1 exp

      repulsion_x *= common_factor_a; //1 mult
      repulsion_x_a1 *= common_factor_a1; //1 mult

      repulsion_y *= common_factor_a; //1 mult
      repulsion_y_a1 *= common_factor_a1; //1 mult

      double check = ex_a * repulsion_x + ey_a * repulsion_y;                                             //1 add, 2 mult
      double check_a1 = ex_a_1 * repulsion_x_a1 + ey_a_1 * repulsion_y_a1;                                             //1 add, 2 mult

      double threshold = sqrt(repulsion_x * repulsion_x + repulsion_y * repulsion_y) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
      double threshold_a1 = sqrt(repulsion_x_a1 * repulsion_x_a1 + repulsion_y_a1 * repulsion_y_a1) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt

      double w = -check >= threshold ? 1 : INFLUENCE;
      double w_a1 = -check_a1 >= threshold_a1 ? 1 : INFLUENCE;

      social_force_x += repulsion_x * w; // 1 add => 1 flop
      social_force_y += repulsion_y * w; // 1 add => 1 flop
      social_force_x_1 += repulsion_x_a1 * w_a1;
      social_force_y_1 += repulsion_y_a1 * w_a1;


      //iterate over all people
      
      for (int j = i+2; j < n; j++)
      {
        double rx_ab = rx_a - position[IndexX(j)];
        double ry_ab = ry_a - position[IndexY(j,n)];
        double ex_b = desired_direction[IndexX(j)];
        double ey_b = desired_direction[IndexY(j,n)];
        double delta_b = speed[j];
        double r_ab_norm = sqrt(rx_ab * rx_ab + ry_ab * ry_ab);

        double rx_a1b = rx_a_1 - position[IndexX(j)];
        double ry_a1b = ry_a_1 - position[IndexY(j,n)];
        double r_a1b_norm = sqrt(rx_a1b * rx_a1b + ry_a1b * ry_a1b);

        // from now on compute repulsion term directly for both alpha and beta
        //me stands for "minus e"
        double rx_ab_mex = rx_ab - delta_b * ex_b; //1 add, 1 mult
        double rx_ba_mex = -rx_ab - delta_a * ex_a; //1 add, 1 mult
        double rx_a1b_mex = rx_a1b - delta_b * ex_b; //1 add, 1 mult
        double rx_ba1_mex = -rx_a1b - delta_a_1 * ex_a_1; //1 add, 1 mult

        double ry_ab_mey = ry_ab - delta_b * ey_b; //1 add, 1 mult
        double ry_ba_mey = -ry_ab - delta_a * ey_a; //1 add, 1 mult
        double ry_a1b_mey = ry_a1b - delta_b * ey_b; //1 add, 1 mult
        double ry_ba1_mey = -ry_a1b - delta_a_1 * ey_a_1; //1 add, 1 mult

        double r_ab_me_norm = sqrt(rx_ab_mex * rx_ab_mex + ry_ab_mey * ry_ab_mey); //1 add, 2 mult, 1 sqrt
        double r_ba_me_norm = sqrt(rx_ba_mex * rx_ba_mex + ry_ba_mey * ry_ba_mey); //1 add, 2 mult, 1 sqrt
        double r_a1b_me_norm = sqrt(rx_a1b_mex * rx_a1b_mex + ry_a1b_mey * ry_a1b_mey); //1 add, 2 mult, 1 sqrt
        double r_ba1_me_norm = sqrt(rx_ba1_mex * rx_ba1_mex + ry_ba1_mey * ry_ba1_mey); //1 add, 2 mult, 1 sqrt

        double norm_sum_a = r_ab_norm + r_ab_me_norm;                                //1 add
        double norm_sum_b = r_ab_norm + r_ba_me_norm;                                //1 add
        double norm_sum_a1b = r_a1b_norm + r_a1b_me_norm;                                //1 add
        double norm_sum_ba1 = r_a1b_norm + r_ba1_me_norm;                                //1 add

        double repulsion_x = rx_ab / r_ab_norm + rx_ab_mex / r_ab_me_norm; //1 add, 2 div
        double repulsion_x_b = -rx_ab / r_ab_norm + rx_ba_mex / r_ba_me_norm; //1 add, 2 div
        double repulsion_x_a1b = rx_a1b / r_a1b_norm + rx_a1b_mex / r_a1b_me_norm; //1 add, 2 div
        double repulsion_x_ba1 = -rx_a1b / r_a1b_norm + rx_ba1_mex / r_ba1_me_norm; //1 add, 2 div

        double repulsion_y = ry_ab / r_ab_norm + ry_ab_mey / r_ab_me_norm; //1 add, 2 div
        double repulsion_y_b = -ry_ab / r_ab_norm + ry_ba_mey / r_ba_me_norm; //1 add, 2 div
        double repulsion_y_a1b = ry_a1b / r_a1b_norm + ry_a1b_mey / r_a1b_me_norm; //1 add, 2 div
        double repulsion_y_ba1 = -ry_a1b / r_a1b_norm + ry_ba1_mey / r_ba1_me_norm; //1 add, 2 div


        double b_a = sqrt(norm_sum_a * norm_sum_a - delta_b * delta_b) * 0.5; //1 add, 3 mult, 1 sqrt
        double b_b = sqrt(norm_sum_b * norm_sum_b - delta_a * delta_a) * 0.5; //1 add, 3 mult, 1 sqrt
        double b_a1b = sqrt(norm_sum_a1b * norm_sum_a1b - delta_b * delta_b) * 0.5; //1 add, 3 mult, 1 sqrt
        double b_ba1 = sqrt(norm_sum_ba1 * norm_sum_ba1 - delta_a_1 * delta_a_1) * 0.5; //1 add, 3 mult, 1 sqrt

        double common_factor_a = exp_fast(-b_a * inv_sigma) * norm_sum_a * DIV_FACTOR / b_a; //2 mult, 2 div, 1 exp
        double common_factor_b = exp_fast(-b_b * inv_sigma) * norm_sum_b * DIV_FACTOR / b_b; //2 mult, 2 div, 1 exp
        double common_factor_a1b = exp_fast(-b_a1b * inv_sigma) * norm_sum_a1b * DIV_FACTOR / b_a1b; //2 mult, 2 div, 1 exp
        double common_factor_ba1 = exp_fast(-b_ba1 * inv_sigma) * norm_sum_ba1 * DIV_FACTOR / b_ba1; //2 mult, 2 div, 1 exp

        repulsion_x *= common_factor_a; //1 mult
        repulsion_x_b *= common_factor_b; //1 mult
        repulsion_x_a1b *= common_factor_a1b; //1 mult
        repulsion_x_ba1 *= common_factor_ba1; //1 mult

        repulsion_y *= common_factor_a; //1 mult
        repulsion_y_b *= common_factor_b; //1 mult
        repulsion_y_a1b *= common_factor_a1b; //1 mult
        repulsion_y_ba1 *= common_factor_ba1; //1 mult

        double check = ex_a * repulsion_x + ey_a * repulsion_y;                                             //1 add, 2 mult
        double check_b = ex_b * repulsion_x_b + ey_b * repulsion_y_b;                                             //1 add, 2 mult
        double check_a1b = ex_a_1 * repulsion_x_a1b + ey_a_1 * repulsion_y_a1b;                                             //1 add, 2 mult
        double check_ba1 = ex_b * repulsion_x_ba1 + ey_b * repulsion_y_ba1;                                             //1 add, 2 mult

        double threshold = sqrt(repulsion_x * repulsion_x + repulsion_y * repulsion_y) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
        double threshold_b = sqrt(repulsion_x_b * repulsion_x_b + repulsion_y_b * repulsion_y_b) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
        double threshold_a1b = sqrt(repulsion_x_a1b * repulsion_x_a1b + repulsion_y_a1b * repulsion_y_a1b) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
        double threshold_ba1 = sqrt(repulsion_x_ba1 * repulsion_x_ba1 + repulsion_y_ba1 * repulsion_y_ba1) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt

        double w = -check >= threshold ? 1 : INFLUENCE;
        double w_b = -check_b >= threshold_b ? 1 : INFLUENCE;
        double w_a1b = -check_a1b >= threshold_a1b ? 1 : INFLUENCE;
        double w_ba1 = -check_ba1 >= threshold_ba1 ? 1 : INFLUENCE;

        social_force_x += repulsion_x * w; // 1 add => 1 flop
        social_force_y += repulsion_y * w; // 1 add => 1 flop
        social_force[IndexX(j)] += repulsion_x_b * w_b;
        social_force[IndexY(j,n)] += repulsion_y_b * w_b;
        social_force_x_1 += repulsion_x_a1b * w_a1b; // 1 add => 1 flop
        social_force_y_1 += repulsion_y_a1b * w_a1b; // 1 add => 1 flop
        social_force[IndexX(j)] += repulsion_x_ba1 * w_ba1;
        social_force[IndexY(j,n)] += repulsion_y_ba1 * w_ba1;
      }

      /************************************************/
      //UPDATE ACCELERATION TERM
      /************************************************/
      // get actual velocity, desired direction, desired speed
      // compute velocity difference
      double v_delta_x = desired_speed_value * ex_a; // 1 mul, 1 flop
      double v_delta_y = desired_speed_value * ey_a; // 1 mul, 1 flop
      v_delta_x -= actual_velocity_x;                // 1 add, 1 flop
      v_delta_y -= actual_velocity_y;                // 1 add, 1 flop

      // apply realxation time
      social_force_x += INV_RELAX_TIME * v_delta_x; // 1 mul => 1 flops
      social_force_y += INV_RELAX_TIME * v_delta_y; // 1 mul => 1 flops

      /************************************************/


      // get actual velocity, desired direction, desired speed
      // compute velocity difference
      double v_delta_x_1 = desired_speed_value_1 * ex_a_1; // 1 mul, 1 flop
      double v_delta_y_1 = desired_speed_value_1 * ey_a_1; // 1 mul, 1 flop
      v_delta_x_1 -= actual_velocity_x_1;                // 1 add, 1 flop
      v_delta_y_1 -= actual_velocity_y_1;                // 1 add, 1 flop

      // apply realxation time
      social_force_x_1 += INV_RELAX_TIME * v_delta_x_1; // 1 mul => 1 flops
      social_force_y_1 += INV_RELAX_TIME * v_delta_y_1; // 1 mul => 1 flops

      /************************************************/
      // UPDATE POSITION
      /************************************************/
      //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
      double prefered_velocity_x = actual_velocity_x + social_force_x * TIMESTEP; // 1 add, 1 mult => 2 flops
      double prefered_velocity_y = actual_velocity_y + social_force_y * TIMESTEP; // 1 add, 1 mult => 2 flops

      //compute the norm of the preferd velocity
      double x_sq_plus_y_sq = (prefered_velocity_x * prefered_velocity_x) + (prefered_velocity_y * prefered_velocity_y); // 1 add, 2 mults => 3 flops
      double norm_value = sqrt(x_sq_plus_y_sq);                                                                          // 1 sqrt => 1 flops

      //fromula 12 in the paper --> compute control_value according to norm
      double control_value = norm_value > max_speed ? (max_speed / norm_value) : 1.0; // 1 div => 1 flops

      //apply control value
      prefered_velocity_x *= control_value; // 1 mul, 1 flop
      prefered_velocity_y *= control_value; // 1 mul, 1 flop

      //update position
      rx_a += prefered_velocity_x * TIMESTEP; // 1 add, 1 mul => 2 flops
      ry_a += prefered_velocity_y * TIMESTEP; // 1 add, 1 mul => 2 flops

      //*************************************************************************************

      //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
      double prefered_velocity_x_1 = actual_velocity_x_1 + social_force_x_1 * TIMESTEP; // 1 add, 1 mult => 2 flops
      double prefered_velocity_y_1 = actual_velocity_y_1 + social_force_y_1 * TIMESTEP; // 1 add, 1 mult => 2 flops

      //compute the norm of the preferd velocity
      double x_sq_plus_y_sq_1 = (prefered_velocity_x_1 * prefered_velocity_x_1) + (prefered_velocity_y_1 * prefered_velocity_y_1); // 1 add, 2 mults => 3 flops
      double norm_value_1 = sqrt(x_sq_plus_y_sq_1);                                                                          // 1 sqrt => 1 flops

      //fromula 12 in the paper --> compute control_value according to norm
      double control_value_1 = norm_value_1 > max_speed_1 ? (max_speed_1 / norm_value_1) : 1.0; // 1 div => 1 flops

      //apply control value
      prefered_velocity_x_1 *= control_value_1; // 1 mul, 1 flop
      prefered_velocity_y_1 *= control_value_1; // 1 mul, 1 flop

      //update position
      rx_a_1 += prefered_velocity_x_1 * TIMESTEP; // 1 add, 1 mul => 2 flops
      ry_a_1 += prefered_velocity_y_1 * TIMESTEP; // 1 add, 1 mul => 2 flops

      /************************************************/
      //UPDATE DESIRED DIRECTION
      /************************************************/

      // compute differences
      double delta_x = target_x - rx_a; // 1 add => 1 flop
      double delta_y = target_y - ry_a; // 1 add => 1 flop

      // normalization constant
      double d = delta_x * delta_x + delta_y * delta_y; // 1 add, 2 mult => 3 flops
      double normalizer = sqrt(d);                      // 1 sqrt => 1 flop

      // save desired_direction and position //save speed value, desire direction, actual_velocity
      desired_direction_temp[IndexX(i)] = delta_x / normalizer;     // 1 div => 1 flop
      desired_direction_temp[IndexY(i,n)] = delta_y / normalizer; // 1 div => 1 flop
      position_temp[IndexX(i)] = rx_a;
      position_temp[IndexY(i,n)] = ry_a;
      speed_temp[i] = control_value * norm_value * TIMESTEP; // 1 mul, 1 flop
      actual_velocity[IndexX(i)] = prefered_velocity_x;
      actual_velocity[IndexY(i,n)] = prefered_velocity_y;
      social_force[IndexX(i)] = 0.0;
      social_force[IndexY(i,n)] = 0.0;

      /************************************************/

      // compute differences
      double delta_x_1 = target_x_1 - rx_a_1; // 1 add => 1 flop
      double delta_y_1 = target_y_1 - ry_a_1; // 1 add => 1 flop

      // normalization constant
      double d_1 = delta_x_1 * delta_x_1 + delta_y_1 * delta_y_1; // 1 add, 2 mult => 3 flops
      double normalizer_1 = sqrt(d_1);                      // 1 sqrt => 1 flop

      // save desired_direction and position //save speed value, desire direction, actual_velocity
      desired_direction_temp[IndexX(i+1)] = delta_x_1 / normalizer_1;     // 1 div => 1 flop
      desired_direction_temp[IndexY(i+1,n)] = delta_y_1 / normalizer_1; // 1 div => 1 flop
      position_temp[IndexX(i+1)] = rx_a_1;
      position_temp[IndexY(i+1,n)] = ry_a_1;
      speed_temp[(i+1)] = control_value_1 * norm_value_1 * TIMESTEP; // 1 mul, 1 flop
      actual_velocity[IndexX(i+1)] = prefered_velocity_x_1;
      actual_velocity[IndexY(i+1,n)] = prefered_velocity_y_1;
      social_force[IndexX(i+1)] = 0.0;
      social_force[IndexY(i+1,n)] = 0.0;
    }

    // iterate over every person for step % 2 == 1
    for (int i = 0; i < n; i=i+2)
    {
      /************************************************/
      // load values
      /************************************************/
      double social_force_x = social_force[IndexX(i)];// 0.0;
      double social_force_y = social_force[IndexY(i,n)];// 0.0;
      double rx_a = position_temp[IndexX(i)];
      double ry_a = position_temp[IndexY(i,n)];
      double ex_a = desired_direction_temp[IndexX(i)];
      double ey_a = desired_direction_temp[IndexY(i,n)];
      double actual_velocity_x = actual_velocity[IndexX(i)];
      double actual_velocity_y = actual_velocity[IndexY(i,n)];
      double max_speed = desired_max_speed[i];
      double desired_speed_value = desired_speed[i];
      double target_x = final_destination[IndexX(i)];
      double target_y = final_destination[IndexY(i,n)];

      /************************************************/

      double social_force_x_1 = social_force[IndexX(i+1)];// 0.0;
      double social_force_y_1 = social_force[IndexY(i+1,n)];// 0.0;
      double rx_a_1 = position_temp[IndexX(i+1)];
      double ry_a_1 = position_temp[IndexY(i+1,n)];
      double ex_a_1 = desired_direction_temp[IndexX(i+1)];
      double ey_a_1 = desired_direction_temp[IndexY(i+1,n)];
      double actual_velocity_x_1 = actual_velocity[IndexX(i+1)];
      double actual_velocity_y_1 = actual_velocity[IndexY(i+1,n)];
      double max_speed_1 = desired_max_speed[(i+1)];
      double desired_speed_value_1 = desired_speed[(i+1)];
      double target_x_1 = final_destination[IndexX(i+1)];
      double target_y_1 = final_destination[IndexY(i+1,n)];

       /************************************************/
       // UPDATE BORDER REPULSION TERM
       /************************************************/
       //unrolled borders since there are only 2 in the scenario
      
        
      double ry_aB_0 = ry_a - border_0; //1 add
      double ry_aB_1 = ry_a - border_1; //1 add

      double r_aB_0_norm = ry_aB_0 > 0 ? ry_aB_0 : -ry_aB_0;
      double r_aB_1_norm = ry_aB_1 > 0 ? ry_aB_1 : -ry_aB_1;

      double shared_expression_0 = exp_fast((-r_aB_0_norm) * INV_R) * UTIMESR / r_aB_0_norm; //1 exp, 2 mult, 1 div
      double shared_expression_1 = exp_fast((-r_aB_1_norm) * INV_R) * UTIMESR / r_aB_1_norm; //1 exp, 2 mult, 1 div

      double repulsion_0 = shared_expression_0 * ry_aB_0; //1 mult
      double repulsion_1 = shared_expression_1 * ry_aB_1; //1 mult

      social_force_y += repulsion_0 + repulsion_1; //1 add
      
      /************************************************/

      double ry_a_1B_0 = ry_a_1 - border_0; //1 add
      double ry_a_1B_1 = ry_a_1 - border_1; //1 add

      double r_a_1B_0_norm = ry_a_1B_0 > 0 ? ry_a_1B_0 : -ry_a_1B_0;
      double r_a_1B_1_norm = ry_a_1B_1 > 0 ? ry_a_1B_1 : -ry_a_1B_1;

      double shared_expression_0_1 = exp_fast((-r_a_1B_0_norm) * INV_R) * UTIMESR / r_a_1B_0_norm; //1 exp, 2 mult, 1 div
      double shared_expression_1_1 = exp_fast((-r_a_1B_1_norm) * INV_R) * UTIMESR / r_a_1B_1_norm; //1 exp, 2 mult, 1 div

      double repulsion_0_1 = shared_expression_0_1 * ry_a_1B_0; //1 mult
      double repulsion_1_1 = shared_expression_1_1 * ry_a_1B_1; //1 mult

      social_force_y_1 += repulsion_0_1 + repulsion_1_1; //1 add

      /************************************************/
      // UPDATE PEOPLE REPULSION TERM
      /************************************************/

      // solve common repulsion of alpha0 and alpha1
      double delta_a = speed_temp[i];
      double delta_a_1 = speed_temp[i+1];

      double rx_a0a1 = rx_a - rx_a_1;
      double ry_a0a1 = ry_a - ry_a_1;
      double r_a0a1_norm = sqrt(rx_a0a1 * rx_a0a1 + ry_a0a1 * ry_a0a1);

      // from now on compute repulsion term directly for both alpha and beta
      //me stands for "minus e"
      double rx_a0a1_mex = rx_a0a1 - delta_a_1 * ex_a_1; //1 add, 1 mult
      double rx_a1a0_mex = -rx_a0a1 - delta_a * ex_a; //1 add, 1 mult

      double ry_a0a1_mey = ry_a0a1 - delta_a_1 * ey_a_1; //1 add, 1 mult
      double ry_a1a0_mey = -ry_a0a1 - delta_a * ey_a; //1 add, 1 mult

      double r_a0a1_me_norm = sqrt(rx_a0a1_mex * rx_a0a1_mex + ry_a0a1_mey * ry_a0a1_mey); //1 add, 2 mult, 1 sqrt
      double r_a1a0_me_norm = sqrt(rx_a1a0_mex * rx_a1a0_mex + ry_a1a0_mey * ry_a1a0_mey); //1 add, 2 mult, 1 sqrt

      double norm_sum_a = r_a0a1_norm + r_a0a1_me_norm;                                //1 add
      double norm_sum_a1 = r_a0a1_norm + r_a1a0_me_norm;                                //1 add

      double repulsion_x = rx_a0a1 / r_a0a1_norm + rx_a0a1_mex / r_a0a1_me_norm; //1 add, 2 div
      double repulsion_x_a1 = -rx_a0a1 / r_a0a1_norm + rx_a1a0_mex / r_a1a0_me_norm; //1 add, 2 div

      double repulsion_y = ry_a0a1 / r_a0a1_norm + ry_a0a1_mey / r_a0a1_me_norm; //1 add, 2 div
      double repulsion_y_a1 = -ry_a0a1 / r_a0a1_norm + ry_a1a0_mey / r_a1a0_me_norm; //1 add, 2 div

      double b_a = sqrt(norm_sum_a * norm_sum_a - delta_a_1 * delta_a_1) * 0.5; //1 add, 3 mult, 1 sqrt
      double b_a1 = sqrt(norm_sum_a1 * norm_sum_a1 - delta_a * delta_a) * 0.5; //1 add, 3 mult, 1 sqrt

      double common_factor_a = exp_fast(-b_a * inv_sigma) * norm_sum_a * DIV_FACTOR / b_a; //2 mult, 2 div, 1 exp
      double common_factor_a1 = exp_fast(-b_a1 * inv_sigma) * norm_sum_a1 * DIV_FACTOR / b_a1; //2 mult, 2 div, 1 exp

      repulsion_x *= common_factor_a; //1 mult
      repulsion_x_a1 *= common_factor_a1; //1 mult

      repulsion_y *= common_factor_a; //1 mult
      repulsion_y_a1 *= common_factor_a1; //1 mult

      double check = ex_a * repulsion_x + ey_a * repulsion_y;                                             //1 add, 2 mult
      double check_a1 = ex_a_1 * repulsion_x_a1 + ey_a_1 * repulsion_y_a1;                                             //1 add, 2 mult

      double threshold = sqrt(repulsion_x * repulsion_x + repulsion_y * repulsion_y) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
      double threshold_a1 = sqrt(repulsion_x_a1 * repulsion_x_a1 + repulsion_y_a1 * repulsion_y_a1) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt

      double w = -check >= threshold ? 1 : INFLUENCE;
      double w_a1 = -check_a1 >= threshold_a1 ? 1 : INFLUENCE;

      social_force_x += repulsion_x * w; // 1 add => 1 flop
      social_force_y += repulsion_y * w; // 1 add => 1 flop
      social_force_x_1 += repulsion_x_a1 * w_a1;
      social_force_y_1 += repulsion_y_a1 * w_a1;


      //iterate over all people
      
      for (int j = i+2; j < n; j++)
      {
        double rx_ab = rx_a - position_temp[IndexX(j)];
        double ry_ab = ry_a - position_temp[IndexY(j,n)];
        double ex_b = desired_direction_temp[IndexX(j)];
        double ey_b = desired_direction_temp[IndexY(j,n)];
        double delta_b = speed_temp[j];
        double r_ab_norm = sqrt(rx_ab * rx_ab + ry_ab * ry_ab);

        double rx_a1b = rx_a_1 - position_temp[IndexX(j)];
        double ry_a1b = ry_a_1 - position_temp[IndexY(j,n)];
        double r_a1b_norm = sqrt(rx_a1b * rx_a1b + ry_a1b * ry_a1b);

        // from now on compute repulsion term directly for both alpha and beta
        //me stands for "minus e"
        double rx_ab_mex = rx_ab - delta_b * ex_b; //1 add, 1 mult
        double rx_ba_mex = -rx_ab - delta_a * ex_a; //1 add, 1 mult
        double rx_a1b_mex = rx_a1b - delta_b * ex_b; //1 add, 1 mult
        double rx_ba1_mex = -rx_a1b - delta_a_1 * ex_a_1; //1 add, 1 mult

        double ry_ab_mey = ry_ab - delta_b * ey_b; //1 add, 1 mult
        double ry_ba_mey = -ry_ab - delta_a * ey_a; //1 add, 1 mult
        double ry_a1b_mey = ry_a1b - delta_b * ey_b; //1 add, 1 mult
        double ry_ba1_mey = -ry_a1b - delta_a_1 * ey_a_1; //1 add, 1 mult

        double r_ab_me_norm = sqrt(rx_ab_mex * rx_ab_mex + ry_ab_mey * ry_ab_mey); //1 add, 2 mult, 1 sqrt
        double r_ba_me_norm = sqrt(rx_ba_mex * rx_ba_mex + ry_ba_mey * ry_ba_mey); //1 add, 2 mult, 1 sqrt
        double r_a1b_me_norm = sqrt(rx_a1b_mex * rx_a1b_mex + ry_a1b_mey * ry_a1b_mey); //1 add, 2 mult, 1 sqrt
        double r_ba1_me_norm = sqrt(rx_ba1_mex * rx_ba1_mex + ry_ba1_mey * ry_ba1_mey); //1 add, 2 mult, 1 sqrt

        double norm_sum_a = r_ab_norm + r_ab_me_norm;                                //1 add
        double norm_sum_b = r_ab_norm + r_ba_me_norm;                                //1 add
        double norm_sum_a1b = r_a1b_norm + r_a1b_me_norm;                                //1 add
        double norm_sum_ba1 = r_a1b_norm + r_ba1_me_norm;                                //1 add

        double repulsion_x = rx_ab / r_ab_norm + rx_ab_mex / r_ab_me_norm; //1 add, 2 div
        double repulsion_x_b = -rx_ab / r_ab_norm + rx_ba_mex / r_ba_me_norm; //1 add, 2 div
        double repulsion_x_a1b = rx_a1b / r_a1b_norm + rx_a1b_mex / r_a1b_me_norm; //1 add, 2 div
        double repulsion_x_ba1 = -rx_a1b / r_a1b_norm + rx_ba1_mex / r_ba1_me_norm; //1 add, 2 div

        double repulsion_y = ry_ab / r_ab_norm + ry_ab_mey / r_ab_me_norm; //1 add, 2 div
        double repulsion_y_b = -ry_ab / r_ab_norm + ry_ba_mey / r_ba_me_norm; //1 add, 2 div
        double repulsion_y_a1b = ry_a1b / r_a1b_norm + ry_a1b_mey / r_a1b_me_norm; //1 add, 2 div
        double repulsion_y_ba1 = -ry_a1b / r_a1b_norm + ry_ba1_mey / r_ba1_me_norm; //1 add, 2 div


        double b_a = sqrt(norm_sum_a * norm_sum_a - delta_b * delta_b) * 0.5; //1 add, 3 mult, 1 sqrt
        double b_b = sqrt(norm_sum_b * norm_sum_b - delta_a * delta_a) * 0.5; //1 add, 3 mult, 1 sqrt
        double b_a1b = sqrt(norm_sum_a1b * norm_sum_a1b - delta_b * delta_b) * 0.5; //1 add, 3 mult, 1 sqrt
        double b_ba1 = sqrt(norm_sum_ba1 * norm_sum_ba1 - delta_a_1 * delta_a_1) * 0.5; //1 add, 3 mult, 1 sqrt

        double common_factor_a = exp_fast(-b_a * inv_sigma) * norm_sum_a * DIV_FACTOR / b_a; //2 mult, 2 div, 1 exp
        double common_factor_b = exp_fast(-b_b * inv_sigma) * norm_sum_b * DIV_FACTOR / b_b; //2 mult, 2 div, 1 exp
        double common_factor_a1b = exp_fast(-b_a1b * inv_sigma) * norm_sum_a1b * DIV_FACTOR / b_a1b; //2 mult, 2 div, 1 exp
        double common_factor_ba1 = exp_fast(-b_ba1 * inv_sigma) * norm_sum_ba1 * DIV_FACTOR / b_ba1; //2 mult, 2 div, 1 exp

        repulsion_x *= common_factor_a; //1 mult
        repulsion_x_b *= common_factor_b; //1 mult
        repulsion_x_a1b *= common_factor_a1b; //1 mult
        repulsion_x_ba1 *= common_factor_ba1; //1 mult

        repulsion_y *= common_factor_a; //1 mult
        repulsion_y_b *= common_factor_b; //1 mult
        repulsion_y_a1b *= common_factor_a1b; //1 mult
        repulsion_y_ba1 *= common_factor_ba1; //1 mult

        double check = ex_a * repulsion_x + ey_a * repulsion_y;                                             //1 add, 2 mult
        double check_b = ex_b * repulsion_x_b + ey_b * repulsion_y_b;                                             //1 add, 2 mult
        double check_a1b = ex_a_1 * repulsion_x_a1b + ey_a_1 * repulsion_y_a1b;                                             //1 add, 2 mult
        double check_ba1 = ex_b * repulsion_x_ba1 + ey_b * repulsion_y_ba1;                                             //1 add, 2 mult

        double threshold = sqrt(repulsion_x * repulsion_x + repulsion_y * repulsion_y) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
        double threshold_b = sqrt(repulsion_x_b * repulsion_x_b + repulsion_y_b * repulsion_y_b) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
        double threshold_a1b = sqrt(repulsion_x_a1b * repulsion_x_a1b + repulsion_y_a1b * repulsion_y_a1b) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
        double threshold_ba1 = sqrt(repulsion_x_ba1 * repulsion_x_ba1 + repulsion_y_ba1 * repulsion_y_ba1) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt

        double w = -check >= threshold ? 1 : INFLUENCE;
        double w_b = -check_b >= threshold_b ? 1 : INFLUENCE;
        double w_a1b = -check_a1b >= threshold_a1b ? 1 : INFLUENCE;
        double w_ba1 = -check_ba1 >= threshold_ba1 ? 1 : INFLUENCE;

        social_force_x += repulsion_x * w; // 1 add => 1 flop
        social_force_y += repulsion_y * w; // 1 add => 1 flop
        social_force[IndexX(j)] += repulsion_x_b * w_b;
        social_force[IndexY(j,n)] += repulsion_y_b * w_b;
        social_force_x_1 += repulsion_x_a1b * w_a1b; // 1 add => 1 flop
        social_force_y_1 += repulsion_y_a1b * w_a1b; // 1 add => 1 flop
        social_force[IndexX(j)] += repulsion_x_ba1 * w_ba1;
        social_force[IndexY(j,n)] += repulsion_y_ba1 * w_ba1;
      }

      /************************************************/
      //UPDATE ACCELERATION TERM
      /************************************************/
      // get actual velocity, desired direction, desired speed
      // compute velocity difference
      double v_delta_x = desired_speed_value * ex_a; // 1 mul, 1 flop
      double v_delta_y = desired_speed_value * ey_a; // 1 mul, 1 flop
      v_delta_x -= actual_velocity_x;                // 1 add, 1 flop
      v_delta_y -= actual_velocity_y;                // 1 add, 1 flop

      // apply realxation time
      social_force_x += INV_RELAX_TIME * v_delta_x; // 1 mul => 1 flops
      social_force_y += INV_RELAX_TIME * v_delta_y; // 1 mul => 1 flops

      /************************************************/


      // get actual velocity, desired direction, desired speed
      // compute velocity difference
      double v_delta_x_1 = desired_speed_value_1 * ex_a_1; // 1 mul, 1 flop
      double v_delta_y_1 = desired_speed_value_1 * ey_a_1; // 1 mul, 1 flop
      v_delta_x_1 -= actual_velocity_x_1;                // 1 add, 1 flop
      v_delta_y_1 -= actual_velocity_y_1;                // 1 add, 1 flop

      // apply realxation time
      social_force_x_1 += INV_RELAX_TIME * v_delta_x_1; // 1 mul => 1 flops
      social_force_y_1 += INV_RELAX_TIME * v_delta_y_1; // 1 mul => 1 flops

      /************************************************/
      // UPDATE POSITION
      /************************************************/
      //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
      double prefered_velocity_x = actual_velocity_x + social_force_x * TIMESTEP; // 1 add, 1 mult => 2 flops
      double prefered_velocity_y = actual_velocity_y + social_force_y * TIMESTEP; // 1 add, 1 mult => 2 flops

      //compute the norm of the preferd velocity
      double x_sq_plus_y_sq = (prefered_velocity_x * prefered_velocity_x) + (prefered_velocity_y * prefered_velocity_y); // 1 add, 2 mults => 3 flops
      double norm_value = sqrt(x_sq_plus_y_sq);                                                                          // 1 sqrt => 1 flops

      //fromula 12 in the paper --> compute control_value according to norm
      double control_value = norm_value > max_speed ? (max_speed / norm_value) : 1.0; // 1 div => 1 flops

      //apply control value
      prefered_velocity_x *= control_value; // 1 mul, 1 flop
      prefered_velocity_y *= control_value; // 1 mul, 1 flop

      //update position
      rx_a += prefered_velocity_x * TIMESTEP; // 1 add, 1 mul => 2 flops
      ry_a += prefered_velocity_y * TIMESTEP; // 1 add, 1 mul => 2 flops

      //*************************************************************************************

      //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
      double prefered_velocity_x_1 = actual_velocity_x_1 + social_force_x_1 * TIMESTEP; // 1 add, 1 mult => 2 flops
      double prefered_velocity_y_1 = actual_velocity_y_1 + social_force_y_1 * TIMESTEP; // 1 add, 1 mult => 2 flops

      //compute the norm of the preferd velocity
      double x_sq_plus_y_sq_1 = (prefered_velocity_x_1 * prefered_velocity_x_1) + (prefered_velocity_y_1 * prefered_velocity_y_1); // 1 add, 2 mults => 3 flops
      double norm_value_1 = sqrt(x_sq_plus_y_sq_1);                                                                          // 1 sqrt => 1 flops

      //fromula 12 in the paper --> compute control_value according to norm
      double control_value_1 = norm_value_1 > max_speed_1 ? (max_speed_1 / norm_value_1) : 1.0; // 1 div => 1 flops

      //apply control value
      prefered_velocity_x_1 *= control_value_1; // 1 mul, 1 flop
      prefered_velocity_y_1 *= control_value_1; // 1 mul, 1 flop

      //update position
      rx_a_1 += prefered_velocity_x_1 * TIMESTEP; // 1 add, 1 mul => 2 flops
      ry_a_1 += prefered_velocity_y_1 * TIMESTEP; // 1 add, 1 mul => 2 flops

      /************************************************/
      //UPDATE DESIRED DIRECTION
      /************************************************/

      // compute differences
      double delta_x = target_x - rx_a; // 1 add => 1 flop
      double delta_y = target_y - ry_a; // 1 add => 1 flop

      // normalization constant
      double d = delta_x * delta_x + delta_y * delta_y; // 1 add, 2 mult => 3 flops
      double normalizer = sqrt(d);                      // 1 sqrt => 1 flop

      // save desired_direction and position //save speed value, desire direction, actual_velocity
      desired_direction[IndexX(i)] = delta_x / normalizer;     // 1 div => 1 flop
      desired_direction[IndexY(i,n)] = delta_y / normalizer; // 1 div => 1 flop
      position[IndexX(i)] = rx_a;
      position[IndexY(i,n)] = ry_a;
      speed[i] = control_value * norm_value; // 1 mul, 1 flop
      actual_velocity[IndexX(i)] = prefered_velocity_x;
      actual_velocity[IndexY(i,n)] = prefered_velocity_y;
      social_force[IndexX(i)] = 0.0;
      social_force[IndexY(i,n)] = 0.0;

      /************************************************/

      // compute differences
      double delta_x_1 = target_x_1 - rx_a_1; // 1 add => 1 flop
      double delta_y_1 = target_y_1 - ry_a_1; // 1 add => 1 flop

      // normalization constant
      double d_1 = delta_x_1 * delta_x_1 + delta_y_1 * delta_y_1; // 1 add, 2 mult => 3 flops
      double normalizer_1 = sqrt(d_1);                      // 1 sqrt => 1 flop

      // save desired_direction and position //save speed value, desire direction, actual_velocity
      desired_direction[IndexX(i+1)] = delta_x_1 / normalizer_1;     // 1 div => 1 flop
      desired_direction[IndexY(i+1,n)] = delta_y_1 / normalizer_1; // 1 div => 1 flop
      position[IndexX(i+1)] = rx_a_1;
      position[IndexY(i+1,n)] = ry_a_1;
      speed[(i+1)] = control_value_1 * norm_value_1; // 1 mul, 1 flop
      actual_velocity[IndexX(i+1)] = prefered_velocity_x_1;
      actual_velocity[IndexY(i+1,n)] = prefered_velocity_y_1;
      social_force[IndexX(i+1)] = 0.0;
      social_force[IndexY(i+1,n)] = 0.0;
    }
  }
}