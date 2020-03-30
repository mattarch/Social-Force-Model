#ifndef FOO_H_   /* Include guard */
#define FOO_H_

#ifdef DEBUG
# define CONSOLE_PRINT(x) printf x
#else
# define CONSOLE_PRINT(x) do {} while (0)
#endif

#include <time.h>
#include <stdio.h>

// default values
#define AVG_SPEED 1.34        // in basic scenario this is the desired speed of every person
#define MAX_SPEED 1.742       // maximum speed allowed for every person
#define RELAX_TIME 0.5        // in [seconds]; used for smooth acceleration
#define WALK_WAY_LENGTH 50.0  // in [meters]; walkway dimension in x direction
#define WALK_WAY_WIDTH 10.0   // in [meters]; walkway dimension in y direction; also distance between borders in basic scenario
#define NUMBER_OF_PEOPLE 2    // number of people in the simulation
#define N_FEATURES 7          // number of featutres stored in the People matrix for every person
#define N_BORDERS 2           // number of biarders used in the scenario
#define TIMESTEP 2.0            // in [seconds]; timestep for simulation
#define N_TIMESTEP 30         // number of timesteps being simulated

// parameters model PAGE 8
#define V_ALPHA_BETA 2.1   // in m^{2}s^{-2}
#define SIGMA 0.3          // in m
#define U_ALPHA_B 10.0     // in m^{2}s^{-2}
#define R 0.2              // in m
#define DELTA_T 2.0        // in s
#define PSI 1.75           // in radians
#define INFLUENCE 0.5      // pure (?)       
// add rest of functions



static char filename_global[40];

void initialize_people(double *people, int n);
void initialize_borders(double *borders, int n);
void update_desired_direction(double *People, int n);
void update_acceleration_term(double *People, double *acceleration_terms, double *actual_velocity, int n);
void compute_actual_velocity(double *People, double *actual_velocity, int n);
void update_people_repulsion_term(double *People, double *Repulsion_term, int n);
void update_border_repulsion_term(double *People, double* borders, double *border_repulsion_term, int n, int n_borders);
void compute_social_force(double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term, double *social_force, int n, int n_borders);
void update_position(double *People, double *social_force, double *prefered_velocity, int n);
void run_simulation();

/* function defined in the header file itself */
void output_to_file_initial_state(char*,double*,int,int,int);
void output_to_file_persons(char*,double*,int,int,int);
void output_to_file_constants(char*);
void get_filename();

/*
  This function outputs the initial state of the Person matrix to the given filename
  Parameters: filename: name of the .txt file
                people: array of people
                     n: number of people
            n_features: number of features per person
*/
void output_to_file_initial_state(char* filename, double* people, int n, int n_features, int n_timestep)
{
    FILE *fptr;

    fptr = fopen(filename,"a"); // -a option: data is appended to end of the file; If the file does not exist, it will be created.

    if(!fptr)
    {
        perror("Error reading file in output_to_file_initial_state");
    }

    output_to_file_constants(filename);

    output_to_file_persons(filename, people, n,n_features, n_timestep);
 
}

/*
  This function outputs the defined constants to the given filename
  Parameters: filename: name of the .txt file
*/
void output_to_file_constants(char* filename)
{
    FILE *fptr;

    fptr = fopen(filename,"a"); // -a option: data is appended to end of the file; If the file does not exist, it will be created.

    if(!fptr)
    {
        perror("Error reading file in output_to_file_constants");
        return;
    }

    fprintf(fptr,"17\n"); // output number of Variables
    fprintf(fptr,"AVG_SPEED %f\n",AVG_SPEED);
    fprintf(fptr,"MAX_SPEED %f\n",MAX_SPEED);
    fprintf(fptr,"RELAX_TIME %f\n",RELAX_TIME);
    fprintf(fptr,"WALK_WAY_LENGTH %f\n",WALK_WAY_LENGTH);
    fprintf(fptr,"WALK_WAY_WIDTH %f\n",WALK_WAY_WIDTH);
    fprintf(fptr,"NUMBER_OF_PEOPLE %d\n",(int)NUMBER_OF_PEOPLE);
    fprintf(fptr,"N_FEATURES %d\n",(int)N_FEATURES);
    fprintf(fptr,"N_BORDERS %d\n",(int)N_BORDERS);
    fprintf(fptr,"TIMESTEP %f\n",TIMESTEP);
    fprintf(fptr,"N_TIMESTEP %d\n",(int)N_TIMESTEP);

    fprintf(fptr,"V_ALPHA_BETA %f\n",V_ALPHA_BETA);
    fprintf(fptr,"SIGMA %f\n",SIGMA);
    fprintf(fptr,"U_ALPHA_B %f\n",U_ALPHA_B);
    fprintf(fptr,"R %f\n",R);
    fprintf(fptr,"DELTA_T %f\n",DELTA_T);
    fprintf(fptr,"PSI %f\n",PSI);
    fprintf(fptr,"INFLUENCE %f\n",INFLUENCE);

    fclose(fptr);
}


/*
  This function output the Person matrix to the given filename (needed after initialization)
  Parameters: filename: name of the .txt file
                people: array of people
                     n: number of people
            n_features: number of features per person
*/
void output_to_file_persons(char* filename, double* people, int n, int n_features, int n_timestep)
{
    FILE *fptr;

    fptr = fopen(filename,"a"); // -a option: data is appended to end of the file; If the file does not exist, it will be created.

    if(!fptr)
    {
        perror("Error reading file in output_to_file_persons");
    }

    for(int i = 0; i < n; i++)
    {
    for(int j = 0; j < n_features; j++)
    {
        double x = people[i * n_features + j];
        fprintf(fptr,"%f ",x);
    }
    fprintf(fptr,"\n");
    }

    fclose(fptr);
}
/*
  This function generates a filename with the current timestamp, accessible through filename_global
*/
void get_filename()
{
    struct tm *timenow;
    time_t now = time(NULL);
    timenow = gmtime(&now);

    strftime(filename_global, sizeof(filename_global), "../../test/basic_%Y-%m-%d_%H-%M-%S", timenow);
}

#endif