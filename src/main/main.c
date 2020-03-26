/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

// import stuff
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// debug Makro
// #define DEBUG_MODE

// default values
#define AVG_SPEED 1.34        // in basic scenario this is the desired speed of every person
#define MAX_SPEED 1.742       // maximum speed allowed for every person
#define RELAX_TIME 0.5        // in [seconds]; used for smooth acceleration
#define WALK_WAY_LENGTH 50.0  // in [meters]; walkway dimension in x direction
#define WALK_WAY_WIDTH 10.0   // in [meters]; walkway dimension in y direction; also distance between borders in basic scenario
#define NUMBER_OF_PEOPLE 2    // number of people in the simulation
#define N_FEATURES 7          // number of featutres stored in the People matrix for every person
#define N_BORDERS 2           // number of biarders used in the scenario
#define TIMESTEP 2            // in [seconds]; timestep for simulation
#define N_TIMESTEP 1          // number of timesteps being simulated

// parameters model PAGE 8
#define V_ALPHA_BETA 2.1   // in m^{2}s^{-2}
#define SIGMA 0.3          // in m
#define U_ALPHA_B 10.0     // in m^{2}s^{-2}
#define R 0.2              // in m
#define DELTA_T 2.0        // in s
#define PSI 1.75           // in radians
#define INFLUENCE 0.5      // pure (?)       
// add rest of functions




// function prototypes
void allocate_memory(double **m, int n);
void initialize_people(double *people, int n);
void initialize_borders(double *borders, int n);
void initialize_vector_to_zero(double *m, int n);
void destroy(double *m);
void update_desired_direction(double *People, int n);
void update_acceleration_term(double *People, double *acceleration_terms, double *actual_velocity, int n);
void compute_actual_velocity(double *People, double *actual_velocity, int n);
void update_people_repulsion_term(double *People, double *Repulsion_term, int n);
void update_border_repulsion_term(double *People, double* borders, double *border_repulsion_term, int n, int n_borders);
void compute_social_force(double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term, double *social_force, int n, int n_borders);
void update_position(double *People, double *social_force, double *prefered_velocity, int n);
void run_simulation();


// main function
// int main()
// {
//   printf("Test\n");
  
//   return 0;
// }

// implementation of functions
//------------------------------------------------------------------------------------------
/*
  This function allocates memory for an array of size n filled with doubles.
   after calling this function *m points to the first element of the array.

  Assumptions:  array will be accessed by doubles.

  Parameters: m: pointer to a pointer to a double
              n: size of the array
*/
void allocate_memory(double **m, int n)
{
  *m = (double *)(malloc(n * sizeof(double)));
}

/*
  This function fills the People array with reasonable starting values.

  Assumptions: there is at least one person in the peoples array.
                sidewalk scenario --> half of the people start from left other half from right.
 Parameters:   People: array of people
                     n: number of people
*/
void initialize_people(double *People, int n)
{
  for(int i = 0;i<n;i++)
  {
    // initialize values independant of starting point and target point
    People[i*N_FEATURES +1] = rand()*WALK_WAY_WIDTH/RAND_MAX; // starting position y coordinate
    People[i*N_FEATURES +2] = 0.0; // starting velocity
    People[i*N_FEATURES +4] = 0.0; // starting value for direct_y
    People[i*N_FEATURES +6] = rand()*WALK_WAY_WIDTH/RAND_MAX; // target y coordinate

    if(i%2) // initialize this person to walk from left to right
    {
      People[i*N_FEATURES ] = 0.0; // starting position x coordinate
      People[i*N_FEATURES +3] = 1.0; // starting value for direct_x      
      People[i*N_FEATURES +5] = WALK_WAY_LENGTH; // target x coordinate
    }
    else   // initialize this person to walk from right to left
    {
      People[i*N_FEATURES ] = WALK_WAY_LENGTH; // starting position x coordinate
      People[i*N_FEATURES +3] = -1.0; // starting value for direct_x      
      People[i*N_FEATURES +5] = 0.0; // target x coordinate  
    }
  }
}

/*
  This function fills the borders array with reasonable starting values.

  Assumptions: 
                sidewalk scenario --> two horizontal borders, bottom one on height 0 top one on height d
 Parameters:   doreders: array of borders
                     n: number of borders
*/
void initialize_borders(double *borders, int n)
{
  // check for the basic scenario
  if(n!=2)
  {
    printf("There are more than 2 borders, borders get initialized to 0.");
    for(int i =0;i<n;i++)
    {
      borders[i] = 0.0;
    }
  }
  else
  {
    borders[0] = WALK_WAY_WIDTH;
    borders[1] = 0.0;
  }
  
}

/*
  This function fills the  array with zeros as starting values.

  Assumptions: none
 Parameters:   m: array
               n: length of m
*/
void initialize_vector_to_zero(double *m, int n)
{
  for(int i = 0; i<n;i++)
  {
    m[i] = 0.0;
  }
}
/*
  This function frees the memory allocated in location m

  Assumptions:  *m was recieved form a previous malloc or calloc function

  Parameters: m: initial address of memory to be freed
*/
void destroy(double *m)
{
  free(m);
}

/*
  This function updates the desired direction for all Persons in the People array.
  This function corresponds to formula (1) from the paper

  Assumptions: only one edge in the direction polygon (use target People array)

  Parameters:   People: array of people
                     n: number of people
*/
void update_desired_direction(double *People, int n)
{
  // iterate over all persons and update direction
  for(int i=0; i<n; i++)
  {
    // get current position and target
    double current_x = People[N_FEATURES*i];
    double current_y = People[N_FEATURES*i + 1];
    double target_x = People[N_FEATURES*i + 5];
    double target_y = People[N_FEATURES*i + 6];
    
    // compute differences
    double delta_x  = target_x - current_x;
    double delta_y = target_y - current_y;

    // normalization constant
    double d = delta_x * delta_x + delta_y * delta_y;
    double normalizer = sqrt(d);

    // update direction
    People[N_FEATURES*i + 3] = delta_x / normalizer;
    People[N_FEATURES*i + 4] = delta_y / normalizer;
  }
}

/*
  This function computes the actual velocity for all Persons in the People array and stors it in the actual velocity array.
  This function is part of formula (2) from the paper

  Assumptions: none

  Parameters:   People: array of people
      actual_velocity: array of size n X 2 for every person
                     n: number of people
*/
void compute_actual_velocity(double *People, double *actual_velocity, int n)
{
  // compute actual velocity for every person  
  // iterate over all people
  for(int i=0; i<n; i++)
  {
    actual_velocity[2*i] = People[N_FEATURES*i + 2]*People[N_FEATURES*i + 3];
    actual_velocity[2*i + 1] = People[N_FEATURES*i + 2]*People[N_FEATURES*i + 4];
  }
}

/*
  This function updates the acceleration term for all Persons in the People array.
  This function is part of formula (2) from the paper

  Assumptions: The RELAX_TIME is never 0.
               The function compute_actual_velocity and update_desired_direction was called in advance in THIS order!!!

  Parameters:   People: array of people
    acceleration_terms: array of x- and y-acceleration for every person
      actual_velocity: array with the actual velocity of every person
                     n: number of people
*/
void update_acceleration_term(double *People, double *acceleration_terms, double *actual_velocity, int n)
{
  //!ATTENTION: function compute_actual_velocity and uupdate_desired_direction have to be called befor this function in this order

  // compute the new acceleration terms for every person
  // iterate over every person
  for(int i=0; i<n; i++)
  {
    // compute velocity difference
    acceleration_terms[2*i] = AVG_SPEED*People[N_FEATURES*i + 3] - actual_velocity[2*i];
    acceleration_terms[2*i + 1] = AVG_SPEED*People[N_FEATURES*i + 4] - actual_velocity[2*i + 1];

    // apply realxation time
    acceleration_terms[2*i] = (1/RELAX_TIME)*acceleration_terms[2*i];
    acceleration_terms[2*i + 1] = (1/RELAX_TIME)*acceleration_terms[2*i + 1];
  }
}


/*
TODO: ARE WE SURE THE DESIRED DIRECTION (e) IS EQUAL TO THE REAL DIRECTION?
*/

/*
  This function updates the repulsion between every pair of people in the 
  set wrt the relative position.
  This function corresponds to formulae (4), (7) and (8) from the paper

  Assumptions: none

  Parameters:   People: array of people
        repulsion_term: matrix containing the force of repulsion between a and b
                     n: number of people
*/
void update_people_repulsion_term(double *People, double *Repulsion_term, int n)
{
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      if(i==j)continue;
      double rx_ab                = People[i * N_FEATURES] - People[j * N_FEATURES];
      double ry_ab                = People[i * N_FEATURES + 1] - People[j * N_FEATURES + 1];
      double ex_a                 = People[i * N_FEATURES + 3];
      double ey_a                 = People[i * N_FEATURES + 4];

      double r_ab_norm            = sqrt(rx_ab * rx_ab + ry_ab * ry_ab);                      //(1)

      //me stands for "minus e" 
      double rx_ab_mex            = rx_ab 
                                      - People[j * N_FEATURES + 2] * DELTA_T * People[j * N_FEATURES + 3];
      double ry_ab_mey            = ry_ab 
                                      - People[j * N_FEATURES + 2] * DELTA_T * People[j * N_FEATURES + 4];

      double r_ab_me_norm         = sqrt(rx_ab_mex * rx_ab_mex + ry_ab_mey * ry_ab_mey);      //(2)

      double repulsion_x          = rx_ab / r_ab_norm + rx_ab_mex / r_ab_me_norm;
      double repulsion_y          = ry_ab / r_ab_norm + ry_ab_mey / r_ab_me_norm;

      double b                    = sqrt((r_ab_norm + r_ab_me_norm) * (r_ab_norm + r_ab_me_norm) 
                                      - ((People[j * N_FEATURES + 2] * DELTA_T) * (People[j * N_FEATURES + 2] * DELTA_T))) / 2;
      
      repulsion_x                *= exp(-b / SIGMA) * (r_ab_norm + r_ab_me_norm);
      repulsion_x                *= V_ALPHA_BETA / 4.0 / SIGMA / b;

      repulsion_y                *= exp(-b / SIGMA) * (r_ab_norm + r_ab_me_norm);
      repulsion_y                *= V_ALPHA_BETA / 4.0 / SIGMA / b;

      double check                = ex_a * (-1.0 * repulsion_x) + ey_a * (-1.0 * repulsion_y);
      double threshold            = sqrt(repulsion_x * repulsion_x + repulsion_y * repulsion_y) * cos(PSI);
      double w                    = check >= threshold ? 1 : INFLUENCE;

      Repulsion_term[i * (2*n) + 2 * j]     = w * repulsion_x;
      Repulsion_term[i * (2*n) + 2 * j + 1] = w * repulsion_y;
    }
  }
}


/*
  This function updates the repulsion between every person and every boarder.
  Here the border B is assumed to be a sidewalk.
  This function corresponds to formula (5) from the paper

  Assumptions: Border B is a straight sidewalk (walking direction east-west), sidewalk described by two borders, a northern and southern border

  Parameters:   People: array of people
               borders: b[0] contains the northern border, b[1] contains the southern border of the sidewalk      
        repulsion_term: matrix containing the force of repulsion between pedestrain i and border j
                     n: number of people
             n_borders: number of borders
*/
void update_border_repulsion_term(double *People, double* borders, double *border_repulsion_term, int n, int n_borders)
{
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n_borders; j++)
    {
      double rx_a                 = People[i * N_FEATURES];
      double ry_a                 = People[i * N_FEATURES + 1];

      double rx_aB                = 0.0;
      double ry_aB                = ry_a -  borders[j];

      double r_aB_norm            = fabs(ry_aB);

      double repulsion_x          =  exp(-r_aB_norm / R) * (rx_aB/r_aB_norm);
      repulsion_x                 *= U_ALPHA_B / (double) R;
      
      double repulsion_y          =  exp(-r_aB_norm / R) * (ry_aB/r_aB_norm);
      repulsion_y                 *= U_ALPHA_B / (double) R;

      border_repulsion_term[i * n + 2 * j]     = repulsion_x;
      border_repulsion_term[i * n + 2 * j + 1] = repulsion_y;
    }
  }
}

/*
  This function computes the social force for each person and stores the results in the array soacial_force
  This function corresponds to formula (9) of the paper

  Assumption: acceleration, people, and border terms have been computed before

  input:      acceleration_term: array with acceleration terms
          people_repulsion_term: array with people repulsion terms
          border_repulsion_term: array with border repulsion terms
                    social_fore: array containing the social forces
                              n: number of people
                      n_borders: number of borders

*/
void compute_social_force(double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term, double *social_force, int n, int n_borders)
{
  // compute the social force for each person
  for(int p = 0; p < n; p++)
  {
    // acceleration term
    social_force[2 * p]     = acceleration_term[2 * p];
    social_force[2 * p + 1] = acceleration_term[2 * p + 1];

    // add repulsive terms toward other people
    for(int beta = 0;  beta < n; beta++)
    {
      // leave out term if beta = p
      if(beta == p){
        continue;
      }

      // add repulsive term towards person beta
      social_force[2 * p]     += people_repulsion_term[p * n + 2 * beta];
      social_force[2 * p + 1] += people_repulsion_term[p * n + 2 * beta + 1];
    }

    // add repulsive terms of borders
    for(int b = 0; b < n_borders; b++)
    {
      social_force[2 * p]     += border_repulsion_term[p * n + 2 * b];
      social_force[2 * p + 1] += border_repulsion_term[p * n + 2 * b + 1];
    }
  }
}

/*
  This function computes the new velocity according to the social force and updates the position of every person
  It implements formulas 10 to 12 in the paper
  Assumptions: sozial force is computed before calling this fuinction

  input: People = array of people
        sozial_force = the sozial force vectors for every person
        prefered_velocity = just a 2*n array for storing the prefered velocity
        n = number of people
*/
void update_position(double *People, double *social_force, double *prefered_velocity, int n)
{
  double control_value;
  double norm_value;
  for(int i = 0; i<n;i++){
    control_value             = 1.0;
    //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
    prefered_velocity[2*i]    = People[N_FEATURES*i +2]*People[N_FEATURES*i +3] + social_force[2*i]*TIMESTEP;
    prefered_velocity[2*i+1]  = People[N_FEATURES*i +2]*People[N_FEATURES*i +4] + social_force[2*i+1]*TIMESTEP;

    //compute the norm of the preferd velocity
    norm_value                = sqrt(pow(prefered_velocity[2*i],2) + pow(prefered_velocity[2*i + 1],2));

    //fromula 12 in the paper --> compute control_value according to norm 
    if(norm_value > MAX_SPEED)
      {
        control_value         = MAX_SPEED/norm_value;
      }
    
    //apply control value
    prefered_velocity[2*i]    *=control_value;
    prefered_velocity[2*i +1] *=control_value;

    //update speed term in People matrix --> this is the new speed
    People[N_FEATURES*i + 2]  = sqrt(pow(prefered_velocity[2*i],2) + pow(prefered_velocity[2*i + 1],2));

    //update position QUESTION: should I use the computed velocity or should I use the updated speed times the desired direction? 
    // they might not be equal because the social force term is not included in the desired direction of movement
    People[N_FEATURES*i]      += prefered_velocity[2*i]*TIMESTEP;
    People[N_FEATURES*i + 1]  += prefered_velocity[2*i+1]*TIMESTEP;
  }

}
/*
  This function runs the simulation 
*/
void run_simulation()
{
  // allocate memory
  double *People                  = (double *) calloc(NUMBER_OF_PEOPLE * N_FEATURES, sizeof(double));
  double *borders                 = (double *) calloc(N_BORDERS, sizeof(double));
  double *actual_velocity         = (double *) calloc(NUMBER_OF_PEOPLE * 2, sizeof(double));
  double *acceleration_term       = (double *) calloc(NUMBER_OF_PEOPLE * 2, sizeof(double));
  double *people_repulsion_term   = (double *) calloc(NUMBER_OF_PEOPLE * NUMBER_OF_PEOPLE * 2, sizeof(double));
  double *border_repulsion_term   = (double *) calloc(NUMBER_OF_PEOPLE * N_BORDERS * 2, sizeof(double));
  double *social_force            = (double *) calloc(NUMBER_OF_PEOPLE * 2, sizeof(double));
  double *prefered_velocity       = (double *) calloc(NUMBER_OF_PEOPLE * 2, sizeof(double));

  // check if calloc worked correctly
  if (People == NULL || borders == NULL || actual_velocity == NULL || acceleration_term == NULL 
  || people_repulsion_term == NULL || border_repulsion_term == NULL || social_force == NULL)
  {
    printf("Error: calloc failed\n");
    return;
  }

  // initialize arrays
  initialize_people(People, NUMBER_OF_PEOPLE);
  initialize_borders(borders, N_BORDERS);

  // start simulation
  printf("Start simulation with %d persons\n", NUMBER_OF_PEOPLE);

  // simulate steps
  for(int step = 0; step < N_TIMESTEP; step++)
  {
    // update variables
    compute_actual_velocity(People, actual_velocity, NUMBER_OF_PEOPLE);
    update_desired_direction(People, NUMBER_OF_PEOPLE);
    update_acceleration_term(People, acceleration_term, actual_velocity, NUMBER_OF_PEOPLE);
    update_people_repulsion_term(People, people_repulsion_term, NUMBER_OF_PEOPLE);
    update_border_repulsion_term(People, borders, border_repulsion_term, NUMBER_OF_PEOPLE, N_BORDERS);
    compute_social_force(acceleration_term, people_repulsion_term, border_repulsion_term, social_force, NUMBER_OF_PEOPLE, N_BORDERS);
    update_position(People, social_force, prefered_velocity, NUMBER_OF_PEOPLE);

    printf("Finished iteration %d\n", (step+1));
  }

  printf("Simulation terminated\n");
}