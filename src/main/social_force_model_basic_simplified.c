/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#include <stdio.h>
#include <math.h>

#include "social_force_model_basic_simplified.h"
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
void update_desired_direction_simplified(float *position, float *final_destination, float *desired_direction, int n)
{
  // iterate over all persons and update desired_direction
  for (int i = 0; i < n; i++)
  {
    // get current position and target
    float current_x = position[IndexX(i)];
    float current_y = position[IndexY(i,n)];
    float target_x = final_destination[IndexX(i)];
    float target_y = final_destination[IndexY(i,n)];

    // compute differences
    float delta_x = target_x - current_x; // 1 add => 1 flop
    float delta_y = target_y - current_y; // 1 add => 1 flop

    // normalization constant
    float d = delta_x * delta_x + delta_y * delta_y; // 1 add, 2 mult => 3 flops
    float normalizer = sqrt(d);                      // 1 sqrt => 1 flop

    // update desired_direction
    desired_direction[IndexX(i)] = delta_x / normalizer;     // 1 div => 1 flop
    desired_direction[IndexY(i,n)] = delta_y / normalizer; // 1 div => 1 flop
  }
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
void update_acceleration_term_simplified(float *desired_direction, float *acceleration_term, float *actual_velocity, float *desired_speed, int n)
{
  //!ATTENTION: function compute_actual_velocity and uupdate_desired_direction have to be called befor this function in this order

  // compute the new acceleration terms for every person
  // iterate over every person
  for (int i = 0; i < n; i++)
  {
    // get actual velocity, desired direction, desired speed
    float actual_velocity_x = actual_velocity[IndexX(i)];
    float actual_velocity_y = actual_velocity[IndexY(i,n)];
    float desired_direction_x = desired_direction[IndexX(i)];
    float desired_direction_y = desired_direction[IndexY(i,n)];
    float desired_speed_value = desired_speed[i];

    // compute velocity difference
    float v_delta_x = desired_speed_value * desired_direction_x; // 1 mul, 1 flop
    float v_delta_y = desired_speed_value * desired_direction_y; // 1 mul, 1 flop
    v_delta_x -= actual_velocity_x;                               // 1 add, 1 flop
    v_delta_y -= actual_velocity_y;                               // 1 add, 1 flop

    // apply realxation time
    acceleration_term[IndexX(i)] = INV_RELAX_TIME * v_delta_x;     // 1 mul => 1 flops
    acceleration_term[IndexY(i,n)] = INV_RELAX_TIME * v_delta_y; // 1 mul => 1 flops
  }
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
void update_people_repulsion_term_simplified(float *position, float *desired_direction, float *actual_speed, float *people_repulsion_term, int n)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i == j)
        continue;
      float rx_ab = position[IndexX(i)] - position[IndexX(j)];         //1 add
      float ry_ab = position[IndexY(i,n)] - position[IndexY(j,n)]; //1 add
      float ex_a = desired_direction[IndexX(i)];
      float ey_a = desired_direction[IndexY(i,n)];
      float ex_b = desired_direction[IndexX(j)];
      float ey_b = desired_direction[IndexY(j,n)];
      float vb = actual_speed[j];
      float delta_b = vb * TIMESTEP; //1 mult

      float r_ab_norm = sqrt(rx_ab * rx_ab + ry_ab * ry_ab); //1 add, 2 mult, 1 sqrt

      //me stands for "minus e"
      float rx_ab_mex = rx_ab - delta_b * ex_b; //1 add, 1 mult
      float ry_ab_mey = ry_ab - delta_b * ey_b; //1 add, 1 mult

      float r_ab_me_norm = sqrt(rx_ab_mex * rx_ab_mex + ry_ab_mey * ry_ab_mey); //1 add, 2 mult, 1 sqrt
      float norm_sum = r_ab_norm + r_ab_me_norm;                                //1 add

      float repulsion_x = rx_ab / r_ab_norm + rx_ab_mex / r_ab_me_norm; //1 add, 2 div
      float repulsion_y = ry_ab / r_ab_norm + ry_ab_mey / r_ab_me_norm; //1 add, 2 div

      float b = sqrt(norm_sum * norm_sum - delta_b * delta_b) / 2; //1 add, 2 mult, 1 div, 1 sqrt

      float common_factor = exp_fast(-b / SIGMA) * norm_sum * DIV_FACTOR / b; //2 mult, 2 div, 1 exp

      repulsion_x *= common_factor; //1 mult
      repulsion_y *= common_factor; //1 mult

      float check = ex_a * repulsion_x + ey_a * repulsion_y;                                             //1 add, 2 mult
      float threshold = sqrt(repulsion_x * repulsion_x + repulsion_y * repulsion_y) * PROJECTION_FACTOR; //1 add, 3 mult, 1 sqrt
      float w = -check >= threshold ? 1 : INFLUENCE;

      people_repulsion_term[IndexX_matrix(i,j,n)] = w * repulsion_x;         //1 mult
      people_repulsion_term[IndexY_matrix(i,j,n)] = w * repulsion_y; //1 mult
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
void update_border_repulsion_term_simplified(float *position, float *borders, float *border_repulsion_term, int n, int n_borders)
{
  for (int j = 0; j < n_borders; j++)
  {
    for (int i = 0; i < n; i++)
    {

      float rx_a = position[IndexX(i)];
      float ry_a = position[IndexY(i,n)];

      float rx_aB = 0.0;
      float ry_aB = ry_a - borders[j]; //1 add => 1 flop

      float r_aB_norm = ry_aB > 0 ? ry_aB : -ry_aB;

      float shared_expression = exp_fast((-r_aB_norm) / R) * U_ALPHA_B / R / r_aB_norm; // 1 exp, 3 div, 1 mult => 4 flops + 1 exp

      float repulsion_x = shared_expression * rx_aB; // 1 mult => 1 flop

      float repulsion_y = shared_expression * ry_aB; // 1 mult => 1 flop

      //       | xb0 | yb0 | xb1 | yb1 |
      border_repulsion_term[IndexX_border(i,j,n)] = repulsion_x;
      border_repulsion_term[IndexY_border(i,j,n)] = repulsion_y;

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
void compute_social_force_simplified(float *acceleration_term, float *people_repulsion_term, float *border_repulsion_term, float *social_force, int n, int n_borders)
{
  // compute the social force for each person
  for (int p = 0; p < n; p++)
  {
    // acceleration term
    social_force[IndexX(p)] = acceleration_term[IndexX(p)];
    social_force[IndexY(p,n)] = acceleration_term[IndexY(p,n)];

    // add repulsive terms toward other people
    for (int beta = 0; beta < n; beta++)
    {
      // leave out term if beta = p
      if (beta == p)
      {
        continue;
      }

      // add repulsive term towards person beta
      social_force[IndexX(p)] += people_repulsion_term[IndexX_matrix(p,beta,n)];             // 1 add => 1 flop
      social_force[IndexY(p,n)] += people_repulsion_term[IndexY_matrix(p,beta,n)]; // 1 add => 1 flop
    }
  }

  // add repulsive terms of borders
  for (int b = 0; b < n_borders; b++)
  {
    for (int p = 0; p < n; p++)
    {
      social_force[IndexX(p)] += border_repulsion_term[IndexX_border(p,b,n)];         // 1 add => 1 flop
      social_force[IndexY(p,n)] += border_repulsion_term[IndexY_border(p,b,n)]; // 1 add => 1 flop
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
void update_position_simplified(float *position, float *desired_direction, float *actual_speed, float *social_force, float *actual_velocity, float *desired_max_speed, int n)
{
  float control_value;
  float norm_value;
  for (int i = 0; i < n; i++)
  {

    //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
    float prefered_velocity_x = actual_velocity[IndexX(i)] + social_force[IndexX(i)] * TIMESTEP;         // 1 add, 1 mult => 2 flops
    float prefered_velocity_y = actual_velocity[IndexY(i,n)] + social_force[IndexY(i,n)] * TIMESTEP; // 1 add, 1 mult => 2 flops

    //compute the norm of the preferd velocity
    float x_sq_plus_y_sq = (prefered_velocity_x * prefered_velocity_x) + (prefered_velocity_y * prefered_velocity_y); // 1 add, 2 mults => 3 flops
    norm_value = sqrt(x_sq_plus_y_sq);                                                                                 // 1 sqrt => 1 flops

    //fromula 12 in the paper --> compute control_value according to norm
    float max_speed = desired_max_speed[i];
    control_value = norm_value > max_speed ? (max_speed / norm_value) : 1.0; // 1 div => 1 flops

    //apply control value
    prefered_velocity_x *= control_value; // 1 mul, 1 flop
    prefered_velocity_y *= control_value; // 1 mul, 1 flop

    //update speed value, desire direction, actual_velocity
    actual_speed[i] = control_value * norm_value;                     // 1 mul, 1 flop
    desired_direction[IndexX(i)] = prefered_velocity_x / actual_speed[i];     // 1 div, 1 flop
    desired_direction[IndexY(i,n)] = prefered_velocity_y / actual_speed[i]; // 1 div, 1 flop
    actual_velocity[IndexX(i)] = prefered_velocity_x;
    actual_velocity[IndexY(i,n)] = prefered_velocity_y;
    //update position
    position[IndexX(i)] += prefered_velocity_x * TIMESTEP;     // 1 add, 1 mul => 2 flops
    position[IndexY(i,n)] += prefered_velocity_y * TIMESTEP; // 1 add, 1 mul => 2 flops
  }
}

void simulation_basic_simplified(int number_of_people, int n_timesteps, float *position, float *speed, float *desired_direction, float *final_destination, float *borders, float *actual_velocity, float *acceleration_term,
                                 float *people_repulsion_term, float *border_repulsion_term, float *social_force, float *desired_speed, float *desired_max_speed)
{
  // start simulation
  CONSOLE_PRINT(("Start simulation with %d persons\n", number_of_people));

  // simulate steps
  for (int step = 0; step < n_timesteps; step++)
  {
    // update variables
    update_desired_direction_simplified(position, final_destination, desired_direction, number_of_people);
    update_acceleration_term_simplified(desired_direction, acceleration_term, actual_velocity, desired_speed, number_of_people);
    update_people_repulsion_term_simplified(position, desired_direction, speed, people_repulsion_term, number_of_people);
    update_border_repulsion_term_simplified(position, borders, border_repulsion_term, number_of_people, N_BORDERS);
    compute_social_force_simplified(acceleration_term, people_repulsion_term, border_repulsion_term, social_force, number_of_people, N_BORDERS);
    update_position_simplified(position, desired_direction, speed, social_force, actual_velocity, desired_max_speed, number_of_people);
    CONSOLE_PRINT(("Finished iteration %d\n", (step + 1)));
  }

  CONSOLE_PRINT(("Simulation terminated\n"));
}
/*
void test_simulation_basic_simplified(int number_of_people, int n_timesteps, float *position, float *speed, float *desired_direction, float *final_destination, float *borders, float *actual_velocity, float *acceleration_term,
                                      float *people_repulsion_term, float *border_repulsion_term, float *social_force, float *desired_speed, float *desired_max_speed)
{

  // start simulation
  printf("Start simulation with %d persons\n", number_of_people);

  int ntimesteps = arguments.visual ? n_timesteps : NTESTS_FINITE_DIFFERENCES;
  // simulate steps
  for (int step = 0; step < ntimesteps; step++)
  {
    // update variables
    update_desired_direction_simplified(position, final_destination, desired_direction, number_of_people);
    update_acceleration_term_simplified(desired_direction, acceleration_term, actual_velocity, desired_speed, number_of_people);

    update_people_repulsion_term_simplified(position, desired_direction, speed, people_repulsion_term, number_of_people);
    update_border_repulsion_term_simplified(position, borders, border_repulsion_term, number_of_people, N_BORDERS);

    test_people_repulsion_with_FD(people_repulsion_term, number_of_people, position, desired_direction, speed);
    test_border_repulsion_with_FD(border_repulsion_term, position, borders, N_BORDERS, number_of_people);

    compute_social_force_simplified(acceleration_term, people_repulsion_term, border_repulsion_term, social_force, number_of_people, N_BORDERS);
    update_position_simplified(position, desired_direction, speed, social_force, actual_velocity, desired_max_speed, number_of_people);

    //printf("Finished iteration %d\n", (step + 1));
  }

  printf("Simulation terminated\n");
}
*/