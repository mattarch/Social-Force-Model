/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <immintrin.h>

#include "social_force_model_vectorize_3.h"
#include "../social_force.h"
#include "../testing.h"
#include "../utility.h"

extern char filename_global[80];

__m256d exp_fast_vec_3(__m256d x, __m256d one, __m256d exp_constant)
{
  x = _mm256_fmadd_pd(x, exp_constant, one);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  x = _mm256_mul_pd(x, x);
  return x;
}

/*
  This function updates the desired direction for all people.
  This function corresponds to formula (1) from the paper.

	Cost:  adds: n * 3
				mults: n * 2
				 divs: n * 2
		  	sqrts: n
				Flops: n * 8

  Assumptions: There is only one final destination per person.
  Parameters:
              position: (n,2) : array of 2d position of people
     final_destination: (n,2) : array with 2d coordinate of the final destinations of people
     desired_direction: (n,2) : array of 2d unit vectors pointing from a person's current position 
                                towards the corresponging final_destination
                     n: number of people
*/
void update_desired_direction_vectorize_3(double *position, double *final_destination, double *desired_direction, int n)
{
  __m256d current_x;
  __m256d current_y;
  __m256d target_x;
  __m256d target_y;
  __m256d delta_x;
  __m256d delta_y;
  __m256d delta_x_2;
  __m256d delta_y_2;
  __m256d normalizer;

  __m256d result_x;
  __m256d result_y;

  __m256d one = _mm256_set1_pd(1);

  // iterate over all persons and update desired_direction
  for (int i = 0; i < n - 3; i += 4)
  {
    // get current position and target
    current_x = _mm256_load_pd(position + i); // now xy positions for two persons in register
    target_x = _mm256_load_pd(final_destination + i);
    current_y = _mm256_load_pd(position + n + i); // now xy positions for two persons in register
    target_y = _mm256_load_pd(final_destination + n + i);

    // compute differences
    delta_x = target_x - current_x;
    delta_y = target_y - current_y;

    // compute norm
    delta_x_2 = _mm256_mul_pd(delta_x, delta_x); // square each entry
    delta_y_2 = _mm256_mul_pd(delta_y, delta_y);
    normalizer = _mm256_sqrt_pd(_mm256_add_pd(delta_x_2, delta_y_2));

    normalizer = _mm256_div_pd(one, normalizer);

    result_x = _mm256_mul_pd(delta_x, normalizer);
    result_y = _mm256_mul_pd(delta_y, normalizer);

    _mm256_store_pd(desired_direction + i, result_x);
    _mm256_store_pd(desired_direction + n + i, result_y);
  }

  /*
  // iterate over all persons and update desired_direction
  for (int i = 0; i < n; i++)
  {
    // get current position and target
    double current_x = position[IndexX(i)];
    double current_y = position[IndexY(i,n)];
    double target_x = final_destination[IndexX(i)];
    double target_y = final_destination[IndexY(i,n)];

    // compute differences
    double delta_x = target_x - current_x; // 1 add => 1 flop
    double delta_y = target_y - current_y; // 1 add => 1 flop

    // normalization constant
    double d = delta_x * delta_x + delta_y * delta_y; // 1 add, 2 mult => 3 flops
    double normalizer = sqrt(d);                      // 1 sqrt => 1 flop

    // update desired_direction
    desired_direction[IndexX(i)] = delta_x / normalizer;     // 1 div => 1 flop
    desired_direction[IndexY(i,n)] = delta_y / normalizer; // 1 div => 1 flop
  }
  */
}

/*
  This function updates the acceleration term for all people.
  This function is part of formula (2) from the paper.

  FLOPS = n * (2 adds, 4 mults)

  Assumptions: - The RELAX_TIME macro is never 0.
               - actual_velocity needs to be up to date, 
               - desired_direction needs to be up to date
               ATTENTION: if not actual_velocity and desired_direction are not up to date
                          first call compute_actual_velocity before update_desired_direction
  Parameters:   
     desired_direction: (n,2) : array of 2d unit vectors pointing from a person's current position 
                                towards the corresponging final_destination
     acceleration_term: (n,2) : array of x- and y-acceleration for every person
       actual_velocity: (n,2) : array of 2d velocity vectors for every person
                                actual_velocity = actual_speed * desired_direction
                     n: number of people
*/
void update_acceleration_term_vectorize_3(double *desired_direction, double *acceleration_term, double *actual_velocity, double *desired_speed, int n)
{
  //!ATTENTION: function compute_actual_velocity and update_desired_direction have to be called befor this function in this order

  __m256d actual_velocity_x;
  __m256d actual_velocity_y;
  __m256d desired_direction_x;
  __m256d desired_direction_y;
  __m256d desired_speed_value;
  __m256d v_delta_x;
  __m256d v_delta_y;
  __m256d inv_relax_time_vec = _mm256_set1_pd(INV_RELAX_TIME);

  for (int i = 0; i < n - 3; i += 4)
  {
    actual_velocity_x = _mm256_load_pd(actual_velocity + i);
    actual_velocity_y = _mm256_load_pd(actual_velocity + n + i);
    desired_direction_x = _mm256_load_pd(desired_direction + i);
    desired_direction_y = _mm256_load_pd(desired_direction + n + i);
    desired_speed_value = _mm256_load_pd(desired_speed + i);

    v_delta_x = _mm256_mul_pd(desired_speed_value, desired_direction_x);
    v_delta_y = _mm256_mul_pd(desired_speed_value, desired_direction_y);
    v_delta_x = _mm256_sub_pd(v_delta_x, actual_velocity_x);
    v_delta_y = _mm256_sub_pd(v_delta_y, actual_velocity_y);

    v_delta_x = _mm256_mul_pd(v_delta_x, inv_relax_time_vec);
    v_delta_y = _mm256_mul_pd(v_delta_y, inv_relax_time_vec);

    _mm256_store_pd(acceleration_term + i, v_delta_x);
    _mm256_store_pd(acceleration_term + n + i, v_delta_y);
  }

  /*
  // compute the new acceleration terms for every person
  // iterate over every person
  for (int i = 0; i < n; i++)
  {
    // get actual velocity, desired direction, desired speed
    double actual_velocity_x = actual_velocity[IndexX(i)];
    double actual_velocity_y = actual_velocity[IndexY(i, n)];
    double desired_direction_x = desired_direction[IndexX(i)];
    double desired_direction_y = desired_direction[IndexY(i, n)];
    double desired_speed_value = desired_speed[i];

    // compute velocity difference
    double v_delta_x = desired_speed_value * desired_direction_x; // 1 mul, 1 flop
    double v_delta_y = desired_speed_value * desired_direction_y; // 1 mul, 1 flop
    v_delta_x -= actual_velocity_x;                               // 1 add, 1 flop
    v_delta_y -= actual_velocity_y;                               // 1 add, 1 flop

    // apply realxation time
    acceleration_term[IndexX(i)] = INV_RELAX_TIME * v_delta_x;    // 1 mul => 1 flops
    acceleration_term[IndexY(i, n)] = INV_RELAX_TIME * v_delta_y; // 1 mul => 1 flops
  }
  */
}

/*
  This function updates the repulsion between every pair of people in the 
  set wrt the relative position.
  This function corresponds to formulae (4), (7) and (8) from the paper.

  FLOPS = (n^2 - n) * (12 add, 20 mult, 7 div, 4 sqrt, 1 exp)
  Assumptions: two different people can not be in the same spot at the same time
  Parameters: 
                     position: (n,2) : array of 2d position of people
            desired_direction: (n,2) : array of 2d unit vectors pointing from a person's current position 
                                       towards the corresponging final_destination
                 actual_speed: (n,1) : array of the actual speed for every person
        people_repulsion_term: (2n,2n) : matrix containing the force of repulsion between person i and j
                            n: number of people
*/
void update_people_repulsion_term_vectorize_3(double *position, double *desired_direction, double *actual_speed, double *people_repulsion_term, int n)
{
  __m256d position_i_x;
  __m256d position_i_y;

  __m256d position_j_x;
  __m256d position_j_y;

  __m256d r_ab_x;
  __m256d r_ab_y;
  __m256d e_a_x;
  __m256d e_a_y;
  __m256d e_b_x;
  __m256d e_b_y;

  __m256d vb;
  __m256d delta_b;

  __m256d r_ab_2_x;
  __m256d r_ab_2_y;
  __m256d r_ab_norm;

  __m256d r_ab_me_x;
  __m256d r_ab_me_y;

  __m256d r_ab_me_2_x;
  __m256d r_ab_me_2_y;
  __m256d r_ab_me_norm;

  __m256d norm_sum;
  __m256d repulsion_x;
  __m256d repulsion_y;

  __m256d norm_sum_2, delta_b_2;
  __m256d b;
  __m256d exp;
  __m256d common_factor;
  __m256d check;
  __m256d check_x;
  __m256d check_y;

  __m256d repulsion_2_x;
  __m256d repulsion_2_y;
  __m256d repulsion_norm;
  __m256d threshold;
  __m256d w;
  __m256d mask;

  __m256d timestep_vec = _mm256_set1_pd(-TIMESTEP);
  __m256d minus_sigma_inv_vec = _mm256_set1_pd(-1.0/SIGMA);
  __m256d div_factor_vec = _mm256_set1_pd(DIV_FACTOR);
  __m256d projection_factor_vec = _mm256_set1_pd(PROJECTION_FACTOR);
  __m256d influencer_vec = _mm256_set1_pd(INFLUENCE);

  __m256d one = _mm256_set1_pd(1);
  __m256d half_vec = _mm256_set1_pd(0.5);
  __m256d minus1_vec = _mm256_set1_pd(-1);
  __m256d eps = _mm256_set1_pd(1e-12);
  __m256d exp_constant = _mm256_set1_pd(0.00006103515); // 1 / 16384

  for (int i = 0; i < n; i++)
  {
    position_i_x = _mm256_broadcast_sd(position + i);
    position_i_y = _mm256_broadcast_sd(position + n + i);
    e_a_x = _mm256_broadcast_sd(desired_direction + i);
    e_a_y = _mm256_broadcast_sd(desired_direction + n + i);

    for (int j = 0; j < n - 3; j += 4)
    {
      position_j_x = _mm256_load_pd(position + j);
      position_j_y = _mm256_load_pd(position + n + j);

      r_ab_x = _mm256_sub_pd(position_i_x, position_j_x);
      r_ab_y = _mm256_sub_pd(position_i_y, position_j_y);

      e_b_x = _mm256_load_pd(desired_direction + j);
      e_b_y = _mm256_load_pd(desired_direction + n + j);

      vb = _mm256_load_pd(actual_speed + j);

      delta_b = _mm256_mul_pd(vb, timestep_vec);

      // compute norm r_ab
      r_ab_2_x = _mm256_mul_pd(r_ab_x, r_ab_x);
      r_ab_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_y, r_ab_y, r_ab_2_x));

      r_ab_me_x = _mm256_fmadd_pd(delta_b, e_b_x, r_ab_x);
      r_ab_me_y = _mm256_fmadd_pd(delta_b, e_b_y, r_ab_y);

      // compute norm r_ab_me
      r_ab_me_2_x = _mm256_mul_pd(r_ab_me_x, r_ab_me_x);
      r_ab_me_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

      // sum up norms
      norm_sum = _mm256_add_pd(r_ab_norm, r_ab_me_norm);

      // take care of i==j
      mask = _mm256_cmp_pd(r_ab_norm, eps, _CMP_GE_OQ);
      r_ab_norm = _mm256_blendv_pd(one, r_ab_norm, mask);
      mask = _mm256_cmp_pd(r_ab_me_norm, eps, _CMP_GE_OQ);
      r_ab_me_norm = _mm256_blendv_pd(one, r_ab_me_norm, mask);

      // invert norm to save divs
      r_ab_norm = _mm256_div_pd(one, r_ab_norm);
      r_ab_me_norm = _mm256_div_pd(one, r_ab_me_norm);

      repulsion_x = _mm256_mul_pd(r_ab_x, r_ab_norm);
      repulsion_x = _mm256_fmadd_pd(r_ab_me_x, r_ab_me_norm, repulsion_x);

      repulsion_y = _mm256_mul_pd(r_ab_y, r_ab_norm);
      repulsion_y = _mm256_fmadd_pd(r_ab_me_y, r_ab_me_norm, repulsion_y);

      delta_b_2 = _mm256_mul_pd(delta_b, delta_b);
      b = _mm256_sqrt_pd(_mm256_fmsub_pd(norm_sum, norm_sum, delta_b_2));
      b = _mm256_mul_pd(b, half_vec);

      exp = _mm256_mul_pd(b, minus_sigma_inv_vec);
      exp = exp_fast_vec_3(exp, one, exp_constant);

      common_factor = _mm256_mul_pd(norm_sum, div_factor_vec);
      common_factor = _mm256_div_pd(common_factor, b);
      common_factor = _mm256_mul_pd(exp, common_factor);

      repulsion_x = _mm256_mul_pd(repulsion_x, common_factor);
      repulsion_y = _mm256_mul_pd(repulsion_y, common_factor);

      // compute norm r_ab
      repulsion_2_y = _mm256_mul_pd(repulsion_y, repulsion_y);
      threshold = _mm256_sqrt_pd(_mm256_fmadd_pd(repulsion_x,repulsion_x, repulsion_2_y));

      check_y = _mm256_mul_pd(e_a_y, repulsion_y);
      check = _mm256_fmadd_pd(e_a_x, repulsion_x, check_y);

      threshold = _mm256_mul_pd(threshold, projection_factor_vec);

      mask = _mm256_cmp_pd(_mm256_mul_pd(check, minus1_vec), threshold, _CMP_GE_OQ);

      w = _mm256_blendv_pd(influencer_vec, one, mask);

      repulsion_x = _mm256_mul_pd(w, repulsion_x);
      repulsion_y = _mm256_mul_pd(w, repulsion_y);

      _mm256_store_pd(people_repulsion_term + i * n + j, repulsion_x);
      _mm256_store_pd(people_repulsion_term + n * n + i * n + j, repulsion_y);
    }
  }

  for (int i = 0; i < n; i++)
  {
    people_repulsion_term[i * n + i] = 0;
    people_repulsion_term[n * n + i * n + i] = 0;
  }

  /*
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i == j)
        continue;
      double rx_ab = position[IndexX(i)] - position[IndexX(j)];       //1 add
      double ry_ab = position[IndexY(i, n)] - position[IndexY(j, n)]; //1 add
      double ex_a = desired_direction[IndexX(i)];
      double ey_a = desired_direction[IndexY(i, n)];
      double ex_b = desired_direction[IndexX(j)];
      double ey_b = desired_direction[IndexY(j, n)];
      double vb = actual_speed[j];
      double delta_b = vb * TIMESTEP; //1 mult

      double r_ab_norm = sqrt(rx_ab * rx_ab + ry_ab * ry_ab); //1 add, 2 mult, 1 sqrt

      //me stands for "minus e"
      double rx_ab_mex = rx_ab - delta_b * ex_b; //1 add, 1 mult
      double ry_ab_mey = ry_ab - delta_b * ey_b; //1 add, 1 mult

      double r_ab_me_norm = sqrt(rx_ab_mex * rx_ab_mex + ry_ab_mey * ry_ab_mey); //1 add, 2 mult, 1 sqrt
      double norm_sum = r_ab_norm + r_ab_me_norm;                                //1 add

      double repulsion_x = rx_ab / r_ab_norm + rx_ab_mex / r_ab_me_norm; //1 add, 2 div
      double repulsion_y = ry_ab / r_ab_norm + ry_ab_mey / r_ab_me_norm; //1 add, 2 div

      double b = sqrt(norm_sum * norm_sum - delta_b * delta_b) / 2; //1 add, 2 mult, 1 div, 1 sqrt

      double common_factor = exp_fast(-b / SIGMA) * norm_sum * DIV_FACTOR / b; //2 mult, 2 div, 1 exp

      repulsion_x *= common_factor; //1 mult
      repulsion_y *= common_factor; //1 mult

      double check = ex_a * repulsion_x + ey_a * repulsion_y;                                             //1 add, 2 mult
      double threshold = sqrt(repulsion_x * repulsion_x + repulsion_y * repulsion_y) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
      double w = -check >= threshold ? 1 : INFLUENCE;

      people_repulsion_term[IndexX_matrix(i, j, n)] = w * repulsion_x; //1 mult
      people_repulsion_term[IndexY_matrix(i, j, n)] = w * repulsion_y; //1 mult
    }
  }
  */
}

/*social_force_model_basic
  This function updates the repulsion between every person and every boarder.
  Here the border B is assumed to be a sidewalk.
  This function corresponds to formula (5) from the paper.

  Cost:  adds: n_borders * n * 1
				mults: n_borders * n * 3
				 divs: n_borders * n * 3
          exp: n_borders * n * 1
        FLOPS: n_borders * n * 8

  Assumptions: The border B is a straight sidewalk (walking direction east-west), sidewalk described by two borders, a northern and southern border
  Parameters:
                 position: (n,2) : array of 2d position of people
                  borders: (1,2) : array of borders for simple sidewalk scenario
                                   b[0] contains the northern border, b[1] contains the southern border of the sidewalk      
    border_repulsion_term: (n, n_borders): matrix containing the force of repulsion between pedestrain i and border j
                        n: number of people
                n_borders: number of borders
*/
void update_border_repulsion_term_vectorize_3(double *position, double *borders, double *border_repulsion_term, int n, int n_borders)
{

  __m256d border;
  __m256d r_a_y;

  __m256d r_aB_y;

  __m256d r_aB_minus_y;

  __m256d r_aB_norm;

  __m256d mask_y;

  __m256d common_factor;
  __m256d exp;

  __m256d zero = _mm256_set1_pd(0);
  __m256d one = _mm256_set1_pd(1);

  __m256d minus1 = _mm256_set1_pd(-1);

  __m256d r_vec = _mm256_set1_pd(R);
  __m256d u_alpha_b_vec = _mm256_set1_pd(U_ALPHA_B);

  __m256d exp_constant = _mm256_set1_pd(0.00006103515); // 1 / 16384

  for (int j = 0; j < 2; j++)
  {
    border = _mm256_broadcast_sd(borders + j);

    for (int i = 0; i < n - 3; i += 4)
    {
      r_a_y = _mm256_load_pd(position + n + i);

      r_aB_y = _mm256_sub_pd(r_a_y, border);

      r_aB_minus_y = _mm256_mul_pd(r_aB_y, minus1);

      mask_y = _mm256_cmp_pd(r_aB_y, zero, _CMP_GE_OQ);
      r_aB_norm = _mm256_blendv_pd(r_aB_minus_y, r_aB_y, mask_y);

      exp = _mm256_div_pd(r_aB_norm, r_vec);
      exp = _mm256_mul_pd(exp, minus1);
      exp = exp_fast_vec_3(exp, one, exp_constant);

      common_factor = _mm256_div_pd(u_alpha_b_vec, r_vec);
      common_factor = _mm256_div_pd(common_factor, r_aB_norm);
      common_factor = _mm256_mul_pd(exp, common_factor);

      common_factor = _mm256_mul_pd(r_aB_y, common_factor);

      _mm256_store_pd(border_repulsion_term + (n + j * 2 * n + i), common_factor);

    } // (1 add, 3 mult, 3 div, 1 exp) * n_borders
  }   // (1 add, 3 mult, 3 div, 1 exp) * n_borders * n

  /*
  for (int j = 0; j < n_borders; j++)
  {
    for (int i = 0; i < n; i++)
    {

      double rx_a = position[IndexX(i)];
      double ry_a = position[IndexY(i, n)];

      double rx_aB = 0.0;
      double ry_aB = ry_a - borders[j]; //1 add => 1 flop

      double r_aB_norm = ry_aB > 0 ? ry_aB : -ry_aB;

      double shared_expression = exp_fast((-r_aB_norm) / R) * U_ALPHA_B / R / r_aB_norm; // 1 exp, 3 div, 1 mult => 4 flops + 1 exp

      double repulsion_x = shared_expression * rx_aB; // 1 mult => 1 flop

      double repulsion_y = shared_expression * ry_aB; // 1 mult => 1 flop

      //       | xb0 | yb0 | xb1 | yb1 |
      border_repulsion_term[IndexX_border(i, j, n)] = repulsion_x;
      border_repulsion_term[IndexY_border(i, j, n)] = repulsion_y;

    } // (1 add, 3 mult, 3 div, 1 exp) * n_borders
  }   // (1 add, 3 mult, 3 div, 1 exp) * n_borders * n
  */
}

/*
  This function computes the social force for each person and stores the results in the array soacial_force
  This function corresponds to formula (9) of the paper.

	Cost:  adds: 2 * n * (n + n_borders social_force_model_basic
  Assumptions: The acceleration, people, and border terms are up to date.
  Parameters:       
             acceleration_term: (n,2) : array of x- and y-acceleration for every person
         people_repulsion_term: (n,n) : matrix containing the force of repulsion between person a and b
         border_repulsion_term: (n, n_borders): matrix containing the force of repulsion between pedestrain a and border b
                  social_force: (n,2) : array containing the social forces for every person
                             n: number of people
                     n_borders: number of borders
*/
void compute_social_force_vectorize_3(double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term, double *social_force, int n, int n_borders)
{

  __m256d social_force_x;
  __m256d social_force_y;

  __m256d border0_x;
  __m256d border0_y;
  __m256d border1_x;
  __m256d border1_y;
  __m256d unpacked;
  __m256d unpacked1;

  __m256d people_repulsion_x;
  __m256d people_repulsion_y;

  __m256d row_add0_x = _mm256_set1_pd(0);
  __m256d row_add1_x = _mm256_set1_pd(0);
  __m256d row_add2_x = _mm256_set1_pd(0);
  __m256d row_add3_x = _mm256_set1_pd(0);

  __m256d row_add0_y = _mm256_set1_pd(0);
  __m256d row_add1_y = _mm256_set1_pd(0);
  __m256d row_add2_y = _mm256_set1_pd(0);
  __m256d row_add3_y = _mm256_set1_pd(0);

  for (int p = 0; p < n - 3; p += 4)
  {

    row_add0_x = _mm256_set1_pd(0);
    row_add1_x = _mm256_set1_pd(0);
    row_add2_x = _mm256_set1_pd(0);
    row_add3_x = _mm256_set1_pd(0);

    row_add0_y = _mm256_set1_pd(0);
    row_add1_y = _mm256_set1_pd(0);
    row_add2_y = _mm256_set1_pd(0);
    row_add3_y = _mm256_set1_pd(0);

    social_force_x = _mm256_load_pd(acceleration_term + p);
    social_force_y = _mm256_load_pd(acceleration_term + n + p);

    border0_y = _mm256_load_pd(border_repulsion_term + n + p);
    border1_y = _mm256_load_pd(border_repulsion_term + 3 * n + p);

    social_force_y = _mm256_add_pd(social_force_y, border0_y);
    social_force_y = _mm256_add_pd(social_force_y, border1_y);

    // add repulsive terms toward other people
    for (int beta = 0; beta < n - 3; beta += 4)
    {
      row_add0_x = _mm256_add_pd(row_add0_x, _mm256_load_pd(people_repulsion_term + p * n + beta));
      row_add1_x = _mm256_add_pd(row_add1_x, _mm256_load_pd(people_repulsion_term + (p + 1) * n + beta));
      row_add2_x = _mm256_add_pd(row_add2_x, _mm256_load_pd(people_repulsion_term + (p + 2) * n + beta));
      row_add3_x = _mm256_add_pd(row_add3_x, _mm256_load_pd(people_repulsion_term + (p + 3) * n + beta));
    }

    row_add0_x = _mm256_hadd_pd(row_add0_x, _mm256_permute2f128_pd(row_add0_x, row_add0_x, 1));
    row_add1_x = _mm256_hadd_pd(row_add1_x, _mm256_permute2f128_pd(row_add1_x, row_add1_x, 1));
    row_add2_x = _mm256_hadd_pd(row_add2_x, _mm256_permute2f128_pd(row_add2_x, row_add2_x, 1));
    row_add3_x = _mm256_hadd_pd(row_add3_x, _mm256_permute2f128_pd(row_add3_x, row_add3_x, 1));

    row_add0_x = _mm256_hadd_pd(row_add0_x, row_add0_x);
    row_add1_x = _mm256_hadd_pd(row_add1_x, row_add1_x);
    row_add2_x = _mm256_hadd_pd(row_add2_x, row_add2_x);
    row_add3_x = _mm256_hadd_pd(row_add3_x, row_add3_x);

    unpacked = _mm256_unpacklo_pd(row_add0_x, row_add1_x);
    unpacked1 = _mm256_unpacklo_pd(row_add2_x, row_add3_x);

    unpacked = _mm256_blend_pd(unpacked, unpacked1, 12);

    social_force_x = _mm256_add_pd(unpacked, social_force_x);
    _mm256_store_pd(social_force + p, social_force_x);

    // add repulsive terms toward other people
    for (int beta = 0; beta < n - 3; beta += 4)
    {
      row_add0_y = _mm256_add_pd(row_add0_y, _mm256_load_pd(people_repulsion_term + n * n + p * n + beta));
      row_add1_y = _mm256_add_pd(row_add1_y, _mm256_load_pd(people_repulsion_term + n * n + (p + 1) * n + beta));
      row_add2_y = _mm256_add_pd(row_add2_y, _mm256_load_pd(people_repulsion_term + n * n + (p + 2) * n + beta));
      row_add3_y = _mm256_add_pd(row_add3_y, _mm256_load_pd(people_repulsion_term + n * n + (p + 3) * n + beta));
    }

    row_add0_y = _mm256_hadd_pd(row_add0_y, _mm256_permute2f128_pd(row_add0_y, row_add0_y, 1));
    row_add1_y = _mm256_hadd_pd(row_add1_y, _mm256_permute2f128_pd(row_add1_y, row_add1_y, 1));
    row_add2_y = _mm256_hadd_pd(row_add2_y, _mm256_permute2f128_pd(row_add2_y, row_add2_y, 1));
    row_add3_y = _mm256_hadd_pd(row_add3_y, _mm256_permute2f128_pd(row_add3_y, row_add3_y, 1));

    row_add0_y = _mm256_hadd_pd(row_add0_y, row_add0_y);
    row_add1_y = _mm256_hadd_pd(row_add1_y, row_add1_y);
    row_add2_y = _mm256_hadd_pd(row_add2_y, row_add2_y);
    row_add3_y = _mm256_hadd_pd(row_add3_y, row_add3_y);

    unpacked = _mm256_unpacklo_pd(row_add0_y, row_add1_y);
    unpacked1 = _mm256_unpacklo_pd(row_add2_y, row_add3_y);

    unpacked = _mm256_blend_pd(unpacked, unpacked1, 12);

    social_force_y = _mm256_add_pd(unpacked, social_force_y);

    _mm256_store_pd(social_force + n + p, social_force_y);
  }

  /*
  // compute the social force for each person
  for (int p = 0; p < n; p++)
  {
    // acceleration term
    social_force[IndexX(p)] = acceleration_term[IndexX(p)];
    social_force[IndexY(p, n)] = acceleration_term[IndexY(p, n)];

    // add repulsive terms toward other people
    for (int beta = 0; beta < n; beta++)
    {
      // leave out term if beta = p
      if (beta == p)
      {
        continue;
      }

      // add repulsive term towards person beta
      social_force[IndexX(p)] += people_repulsion_term[IndexX_matrix(p, beta, n)];    // 1 add => 1 flop
      social_force[IndexY(p, n)] += people_repulsion_term[IndexY_matrix(p, beta, n)]; // 1 add => 1 flop
    }
  }

  // add repulsive terms of borders
  for (int b = 0; b < n_borders; b++)
  {
    for (int p = 0; p < n; p++)
    {
      social_force[IndexX(p)] += border_repulsion_term[IndexX_border(p, b, n)];    // 1 add => 1 flop
      social_force[IndexY(p, n)] += border_repulsion_term[IndexY_border(p, b, n)]; // 1 add => 1 flop
    }
  }
  */
}

/*
  This function computes the new velocity according to the social force and updates the position of every person.
  It implements formulas 10 to 12 in the paper.

  FLOPS = n * (5 adds, 9 mults, 3 divs, 1 sqrts)
        //This is the count if you always execute the if statement.
        
  Assumptions: The social force needs to be computed before calling this function.
  Parameters:
                   position: (n,2) : array of 2d position of people
          desired_direction: (n,2) : array of 2d unit vectors pointing from a person's current position 
                                    towards the corresponging final_destination
               actual_speed: (n,1) : array of the actual speed for every person
               social_force: (n,2) : array containing the social forces for every person
            actual_velocity: (n,2) : array of 2d velocity vectors for every person
                                     actual_velocity = actual_speed * desired_direction
                          n: number of people
*/
void update_position_vectorize_3(double *position, double *desired_direction, double *actual_speed, double *social_force, double *actual_velocity, double *desired_max_speed, int n)
{

  __m256d prefered_velocity_x;
  __m256d prefered_velocity_y;

  __m256d position_x;
  __m256d position_y;

  __m256d prefered_velocity_2_x;
  __m256d prefered_velocity_2_y;

  __m256d prefered_velocity_norm;
  __m256d max_speed;
  __m256d control_value;
  __m256d mask;
  __m256d social_force_x;
  __m256d social_force_y;

  __m128d actual_speed_x;
  __m128d actual_speed_y;

  __m256d actual_speed_xy_256;

  __m256d control_value_prefered_velocity_xy_norm;

  __m256d timestep_vec = _mm256_set1_pd(TIMESTEP);
  __m256d one = _mm256_set1_pd(1);

  for (int i = 0; i < n - 3; i += 4)
  {
    prefered_velocity_x = _mm256_load_pd(actual_velocity + i);
    prefered_velocity_y = _mm256_load_pd(actual_velocity + n + i);

    social_force_x = _mm256_load_pd(social_force + i);
    social_force_y = _mm256_load_pd(social_force + n + i);

    prefered_velocity_x = _mm256_fmadd_pd(social_force_x, timestep_vec, prefered_velocity_x);
    prefered_velocity_y = _mm256_fmadd_pd(social_force_y, timestep_vec, prefered_velocity_y);

    // norm
    prefered_velocity_2_x = _mm256_mul_pd(prefered_velocity_x, prefered_velocity_x);
    prefered_velocity_2_y = _mm256_mul_pd(prefered_velocity_y, prefered_velocity_y);

    prefered_velocity_norm = _mm256_sqrt_pd(_mm256_add_pd(prefered_velocity_2_x, prefered_velocity_2_y));

    max_speed = _mm256_load_pd(desired_max_speed + i);

    mask = _mm256_cmp_pd(prefered_velocity_norm, max_speed, _CMP_GT_OQ);
    control_value = _mm256_blendv_pd(one, _mm256_div_pd(max_speed, prefered_velocity_norm), mask);

    prefered_velocity_x = _mm256_mul_pd(prefered_velocity_x, control_value);
    prefered_velocity_y = _mm256_mul_pd(prefered_velocity_y, control_value);

    _mm256_store_pd(actual_velocity + i, prefered_velocity_x);
    _mm256_store_pd(actual_velocity + n + i, prefered_velocity_y);

    control_value_prefered_velocity_xy_norm = _mm256_mul_pd(control_value, prefered_velocity_norm);

    _mm256_store_pd(actual_speed + i, control_value_prefered_velocity_xy_norm);

    _mm256_store_pd(desired_direction + i, _mm256_div_pd(prefered_velocity_x, control_value_prefered_velocity_xy_norm));
    _mm256_store_pd(desired_direction + n + i, _mm256_div_pd(prefered_velocity_y, control_value_prefered_velocity_xy_norm));

    position_x = _mm256_load_pd(position + i);
    position_y = _mm256_load_pd(position + n + i);

    _mm256_store_pd(position + i, _mm256_fmadd_pd(prefered_velocity_x, timestep_vec, position_x));
    _mm256_store_pd(position + n + i, _mm256_fmadd_pd(prefered_velocity_y, timestep_vec, position_y));
  }

  /*
  double control_value;
  double norm_value;
  for (int i = 0; i < n; i++)
  {

    //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
    double prefered_velocity_x = actual_velocity[IndexX(i)] + social_force[IndexX(i)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
    double prefered_velocity_y = actual_velocity[IndexY(i, n)] + social_force[IndexY(i, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops

    //compute the norm of the preferd velocity
    double x_sq_plus_y_sq = (prefered_velocity_x * prefered_velocity_x) + (prefered_velocity_y * prefered_velocity_y); // 1 add, 2 mults => 3 flops
    norm_value = sqrt(x_sq_plus_y_sq);                                                                                 // 1 sqrt => 1 flops

    //fromula 12 in the paper --> compute control_value according to norm
    double max_speed = desired_max_speed[i];
    control_value = norm_value > max_speed ? (max_speed / norm_value) : 1.0; // 1 div => 1 flops

    //apply control value
    prefered_velocity_x *= control_value; // 1 mul, 1 flop
    prefered_velocity_y *= control_value; // 1 mul, 1 flop

    //update speed value, desire direction, actual_velocity
    actual_speed[i] = control_value * norm_value;                            // 1 mul, 1 flop
    desired_direction[IndexX(i)] = prefered_velocity_x / actual_speed[i];    // 1 div, 1 flop
    desired_direction[IndexY(i, n)] = prefered_velocity_y / actual_speed[i]; // 1 div, 1 flop
    actual_velocity[IndexX(i)] = prefered_velocity_x;
    actual_velocity[IndexY(i, n)] = prefered_velocity_y;
    //update position
    position[IndexX(i)] += prefered_velocity_x * TIMESTEP;    // 1 add, 1 mul => 2 flops
    position[IndexY(i, n)] += prefered_velocity_y * TIMESTEP; // 1 add, 1 mul => 2 flops
  }
  */
}

void simulation_basic_vectorize_3(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination, double *borders, double *actual_velocity, double *acceleration_term,
                                  double *people_repulsion_term, double *border_repulsion_term, double *social_force, double *desired_speed, double *desired_max_speed)
{
  // start simulation
  CONSOLE_PRINT(("Start simulation with %d persons\n", number_of_people));

  // simulate steps
  for (int step = 0; step < n_timesteps; step++)
  {
    // update variables
    update_desired_direction_vectorize_3(position, final_destination, desired_direction, number_of_people);
    update_acceleration_term_vectorize_3(desired_direction, acceleration_term, actual_velocity, desired_speed, number_of_people);
    update_people_repulsion_term_vectorize_3(position, desired_direction, speed, people_repulsion_term, number_of_people);
    update_border_repulsion_term_vectorize_3(position, borders, border_repulsion_term, number_of_people, N_BORDERS);
    compute_social_force_vectorize_3(acceleration_term, people_repulsion_term, border_repulsion_term, social_force, number_of_people, N_BORDERS);
    update_position_vectorize_3(position, desired_direction, speed, social_force, actual_velocity, desired_max_speed, number_of_people);
    CONSOLE_PRINT(("Finished iteration %d\n", (step + 1)));
  }

  CONSOLE_PRINT(("Simulation terminated\n"));
}
