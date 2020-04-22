/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#include <stdio.h>
#include <math.h>
#include <immintrin.h>

#include "social_force_model_vectorize_0.h"
#include "social_force.h"
#include "testing.h"
#include "utility.h"

extern char filename_global[80];
extern struct arguments arguments;

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
void update_desired_direction_vectorize_0(double *position, double *final_destination, double *desired_direction, int n)
{

  __m256d current_xy, target_xy, delta_xy, delta_xy_squared_turned, delta_xy_squared, d, normalizer;

  // iterate over all persons and update desired_direction
  for (int i = 0; i < 2 * n - 3; i += 4)
  {
    // get current position and target
    current_xy = _mm256_load_pd(position + i); // now xy positions for two persons in register
    target_xy = _mm256_load_pd(final_destination + i);

    // compute differences
    delta_xy = target_xy - current_xy;

    // compute norm
    delta_xy_squared = _mm256_mul_pd(delta_xy, delta_xy); // square each entry
    delta_xy_squared_turned = _mm256_permute_pd(delta_xy_squared, 0b0101); //
    d = _mm256_add_pd(delta_xy_squared, delta_xy_squared_turned);
    normalizer = _mm256_sqrt_pd(d);

    // update desired_direction
    delta_xy = _mm256_div_pd(delta_xy, normalizer);

    _mm256_store_pd(desired_direction + i, delta_xy);
  }

  /*

  for (int i = 0; i < n; i++)
  {

    double current_x = position[i * 2];
    double current_y = position[i * 2 + 1];
    double target_x = final_destination[i * 2];
    double target_y = final_destination[i * 2 + 1];

    // compute differences
    double delta_x = target_x - current_x; // 1 add => 1 flop
    double delta_y = target_y - current_y; // 1 add => 1 flop

    // normalization constant
    double d = delta_x * delta_x + delta_y * delta_y; // 1 add, 2 mult => 3 flops
    double normalizer = sqrt(d);                      // 1 sqrt => 1 flop

    // update desired_direction
    desired_direction[i * 2] = delta_x / normalizer;     // 1 div => 1 flop
    desired_direction[i * 2 + 1] = delta_y / normalizer; // 1 div => 1 flop
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
void update_acceleration_term_vectorize_0(double *desired_direction, double *acceleration_term, double *actual_velocity, double *desired_speed, int n)
{
  //!ATTENTION: function compute_actual_velocity and uupdate_desired_direction have to be called befor this function in this order

  // compute the new acceleration terms for every person
  // iterate over 2 persons at a time

  __m256d actual_velocity_xy, desired_direction_xy, desired_speed_value, desired_speed_permuted, v_delta_xy;
  __m256d inv_relax_time_vec = _mm256_set1_pd(INV_RELAX_TIME);
  for (int i = 0; i < 2 * n - 3; i += 4)
  {
    // get actual velocity, desired direction, desired speed
    actual_velocity_xy = _mm256_load_pd(actual_velocity + i);
    desired_direction_xy = _mm256_load_pd(desired_direction + i);
    desired_speed_value = _mm256_load_pd(desired_speed + (i / 2));
    desired_speed_permuted = _mm256_permute4x64_pd(desired_speed_value, 0b01010000);

    // compute velocity difference
    v_delta_xy = _mm256_mul_pd(desired_speed_permuted, desired_direction_xy);
    v_delta_xy = _mm256_sub_pd(v_delta_xy, actual_velocity_xy);

    // apply realxation time
    v_delta_xy = _mm256_mul_pd(v_delta_xy, inv_relax_time_vec);

    _mm256_store_pd(acceleration_term + i, v_delta_xy);
  }

  /*
  for (int i = 0; i < n; i++)
  {
    // get actual velocity, desired direction, desired speed
    double actual_velocity_x = actual_velocity[2 * i];
    double actual_velocity_y = actual_velocity[2 * i + 1];
    double desired_direction_x = desired_direction[2 * i];
    double desired_direction_y = desired_direction[2 * i + 1];
    double desired_speed_value = desired_speed[i];

    // compute velocity difference
    double v_delta_x = desired_speed_value * desired_direction_x; // 1 mul, 1 flop
    double v_delta_y = desired_speed_value * desired_direction_y; // 1 mul, 1 flop
    v_delta_x -= actual_velocity_x;                               // 1 add, 1 flop
    v_delta_y -= actual_velocity_y;                               // 1 add, 1 flop

    // apply realxation time
    acceleration_term[2 * i] = INV_RELAX_TIME * v_delta_x;     // 1 mul => 1 flops
    acceleration_term[2 * i + 1] = INV_RELAX_TIME * v_delta_y; // 1 mul => 1 flops
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
void update_people_repulsion_term_vectorize_0(double *position, double *desired_direction, double *actual_speed, double *people_repulsion_term, int n)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i == j)
        continue;
      double rx_ab = position[i * 2] - position[j * 2];         //1 add
      double ry_ab = position[i * 2 + 1] - position[j * 2 + 1]; //1 add
      double ex_a = desired_direction[i * 2];
      double ey_a = desired_direction[i * 2 + 1];
      double ex_b = desired_direction[j * 2];
      double ey_b = desired_direction[j * 2 + 1];
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

      people_repulsion_term[i * (2 * n) + 2 * j] = w * repulsion_x;     //1 mult
      people_repulsion_term[i * (2 * n) + 2 * j + 1] = w * repulsion_y; //1 mult
    }
  }
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
void update_border_repulsion_term_vectorize_0(double *position, double *borders, double *border_repulsion_term, int n, int n_borders)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n_borders; j++)
    {
      double rx_a = position[i * 2];
      double ry_a = position[i * 2 + 1];

      double rx_aB = 0.0;
      double ry_aB = ry_a - borders[j]; //1 add => 1 flop

      double r_aB_norm = ry_aB > 0 ? ry_aB : -ry_aB;

      double shared_expression = exp_fast((-r_aB_norm) / R) * U_ALPHA_B / R / r_aB_norm; // 1 exp, 3 div, 1 mult => 4 flops + 1 exp

      double repulsion_x = shared_expression * rx_aB; // 1 mult => 1 flop

      double repulsion_y = shared_expression * ry_aB; // 1 mult => 1 flop

      border_repulsion_term[i * (2 * n_borders) + 2 * j] = repulsion_x;
      border_repulsion_term[i * (2 * n_borders) + 2 * j + 1] = repulsion_y;
    } // (1 add, 3 mult, 3 div, 1 exp) * n_borders
  }   // (1 add, 3 mult, 3 div, 1 exp) * n_borders * n
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
void compute_social_force_vectorize_0(double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term, double *social_force, int n, int n_borders)
{
  // compute the social force for each person
  for (int p = 0; p < n; p++)
  {
    // acceleration term
    social_force[2 * p] = acceleration_term[2 * p];
    social_force[2 * p + 1] = acceleration_term[2 * p + 1];

    // add repulsive terms toward other people
    for (int beta = 0; beta < n; beta++)
    {
      // leave out term if beta = p
      if (beta == p)
      {
        continue;
      }

      // add repulsive term towards person beta
      social_force[2 * p] += people_repulsion_term[p * (2 * n) + 2 * beta];         // 1 add => 1 flop
      social_force[2 * p + 1] += people_repulsion_term[p * (2 * n) + 2 * beta + 1]; // 1 add => 1 flop
    }

    // add repulsive terms of borders
    for (int b = 0; b < n_borders; b++)
    {
      social_force[2 * p] += border_repulsion_term[p * (2 * n_borders) + 2 * b];         // 1 add => 1 flop
      social_force[2 * p + 1] += border_repulsion_term[p * (2 * n_borders) + 2 * b + 1]; // 1 add => 1 flop
    }
  }
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
void update_position_vectorize_0(double *position, double *desired_direction, double *actual_speed, double *social_force, double *actual_velocity, double *desired_max_speed, int n)
{
  double control_value;
  double norm_value;
  for (int i = 0; i < n; i++)
  {

    //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
    double prefered_velocity_x = actual_velocity[2 * i] + social_force[2 * i] * TIMESTEP;         // 1 add, 1 mult => 2 flops
    double prefered_velocity_y = actual_velocity[2 * i + 1] + social_force[2 * i + 1] * TIMESTEP; // 1 add, 1 mult => 2 flops

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
    actual_speed[i] = control_value * norm_value;                         // 1 mul, 1 flop
    desired_direction[i * 2] = prefered_velocity_x / actual_speed[i];     // 1 div, 1 flop
    desired_direction[i * 2 + 1] = prefered_velocity_y / actual_speed[i]; // 1 div, 1 flop
    actual_velocity[2 * i] = prefered_velocity_x;
    actual_velocity[2 * i + 1] = prefered_velocity_y;
    //update position
    position[i * 2] += prefered_velocity_x * TIMESTEP;     // 1 add, 1 mul => 2 flops
    position[i * 2 + 1] += prefered_velocity_y * TIMESTEP; // 1 add, 1 mul => 2 flops
  }
}

void simulation_basic_vectorize_0(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination, double *borders, double *actual_velocity, double *acceleration_term,
                                  double *people_repulsion_term, double *border_repulsion_term, double *social_force, double *desired_speed, double *desired_max_speed)
{
  // start simulation
  CONSOLE_PRINT(("Start simulation with %d persons\n", number_of_people));

  // simulate steps
  for (int step = 0; step < n_timesteps; step++)
  {
    // update variables
    update_desired_direction_vectorize_0(position, final_destination, desired_direction, number_of_people);
    update_acceleration_term_vectorize_0(desired_direction, acceleration_term, actual_velocity, desired_speed, number_of_people);
    update_people_repulsion_term_vectorize_0(position, desired_direction, speed, people_repulsion_term, number_of_people);
    update_border_repulsion_term_vectorize_0(position, borders, border_repulsion_term, number_of_people, N_BORDERS);
    compute_social_force_vectorize_0(acceleration_term, people_repulsion_term, border_repulsion_term, social_force, number_of_people, N_BORDERS);
    update_position_vectorize_0(position, desired_direction, speed, social_force, actual_velocity, desired_max_speed, number_of_people);
    CONSOLE_PRINT(("Finished iteration %d\n", (step + 1)));
  }

  CONSOLE_PRINT(("Simulation terminated\n"));
}
