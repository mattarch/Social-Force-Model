#ifndef FOO_H_   /* Include guard */
#define FOO_H_

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

void allocate_memory(double **m, int n);
void initialize_people(double *people, int n);
void initialize_borders(double *borders, int n);
void initialize_vector_to_zero(double *m, int n);
void destroy(double *m);
void update_direction_of_motion(double *People, int n);
void update_acceleration_term(double *People, double *acceleration_terms, double *actual_velocity, int n);
void compute_actual_velocity(double *People, double *actual_velocity, int n);
void update_people_repulsion_term(double *People, double *Repulsion_term, int n);
void update_border_repulsion_term(double *People, double* borders, double *border_repulsion_term, int n, int n_borders);
void compute_social_force(double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term, double *social_force, int n, int n_borders);
void update_position(double *People, double *social_force, double *prefered_velocity, int n);
void run_simulation();

#endif