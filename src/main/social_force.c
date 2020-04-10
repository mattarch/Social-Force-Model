/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

// system includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// constants
#include "social_force.h"
#include "social_force_model_basic.h"

// testing
#include "testing.h"

// parsing arguments
#include "parse_args.h"

// benachmarking
#include "tsc_x86.h"

// utility
#include "utility.h"

// globals
struct arguments arguments = {0};
char filename_global[40];

// main function
int main(int argc, char *argv[])
{

    // default values for arguments
    arguments.n_people = 300;
    arguments.n_timesteps = 300;
    arguments.test = false;
    arguments.walkway_width = 4;
    arguments.walkway_length = 50;

    // parse arguments
    parse_args(argc, argv, &arguments);

    /* start with simulation stuff */

    int number_of_people = arguments.n_people;
    int n_timesteps = arguments.n_timesteps;
    // allocate memory
    double *position = (double *)calloc(number_of_people * 2, sizeof(double));
    double *speed = (double *)calloc(number_of_people, sizeof(double));
    double *desired_direction = (double *)calloc(number_of_people * 2, sizeof(double));
    double *final_destination = (double *)calloc(number_of_people * 2, sizeof(double));
    double *borders = (double *)calloc(N_BORDERS, sizeof(double));
    double *actual_velocity = (double *)calloc(number_of_people * 2, sizeof(double));
    double *acceleration_term = (double *)calloc(number_of_people * 2, sizeof(double));
    double *people_repulsion_term = (double *)calloc(number_of_people * number_of_people * 2, sizeof(double));
    double *border_repulsion_term = (double *)calloc(number_of_people * N_BORDERS * 2, sizeof(double));
    double *social_force = (double *)calloc(number_of_people * 2, sizeof(double));
    double *desired_speed = (double *)calloc(number_of_people, sizeof(double));
    // check if calloc worked correctly
    if (position == NULL || speed == NULL || desired_direction == NULL || final_destination == NULL || borders == NULL || actual_velocity == NULL || acceleration_term == NULL || people_repulsion_term == NULL || border_repulsion_term == NULL || social_force == NULL || desired_speed == NULL)
    {
        printf("Error: calloc failed\n");
        return 1;
    }

    // initialize arrays
    initialize_people(position, desired_direction, final_destination, desired_speed, number_of_people);
    initialize_borders(borders, N_BORDERS);

    if (arguments.test)
    {
        run_tests();
        test_simulation_basic(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed);
    }
    else
    {
        simulation_basic(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed);
    }

    return 0;
}

// implementation of functions
//------------------------------------------------------------------------------------------

/*
  This function initializes the arrays associated with people with reasonable starting values.

  Assumptions: There is at least one person.
               In the sidewalk scenario half of the people start from left the other half from right.
  Parameters: 
              position: (n,2) : array of 2d position of people
     desired_direction: (n,2) : array of 2d unit vectors pointing from a person's current position 
                                towards the corresponging final_destination
     final_destination: (n,2) : array with 2d coordinate of the final destinations of people
                     n: number of people
*/
void initialize_people(double *position, double *desired_direction, double *final_destination, double *desired_speed, int n)
{
    for (int i = 0; i < n; i++)
    {
        // initialize values independant of starting point and target point
        position[i * 2 + 1] = rand() * arguments.walkway_width / RAND_MAX;          // starting position y coordinate
        desired_direction[i * 2 + 1] = 0.0;                                         // starting value for direct_y
        final_destination[i * 2 + 1] = rand() * arguments.walkway_width / RAND_MAX; // target y coordinate
        desired_speed[i] = sampleNormal(0.0676, AVG_SPEED);

        if (i % 2) // initialize this person to walk from left to right
        {
            position[i * 2] = 0.0 - rand() * arguments.walkway_length / RAND_MAX; // starting position x coordinate
            desired_direction[i * 2] = 1.0;                                       // starting value for direct_x
            final_destination[i * 2] = arguments.walkway_length + 10;             // target x coordinate
        }
        else // initialize this person to walk from right to left
        {
            position[i * 2] = arguments.walkway_length + rand() * arguments.walkway_length / RAND_MAX; // starting position x coordinate
            desired_direction[i * 2] = -1.0;                                                           // starting value for direct_x
            final_destination[i * 2] = 0.0 - 10;                                                       // target x coordinate
        }
    }
}

/*
  This function initializes the borders array with reasonable starting values.

  Assumptions: sidewalk scenario --> two horizontal borders, bottom one on height 0 top one on height arguments.walkway_width.
  Parameters:
               borders: (n_borders,2) : array of borders
             n_borders: number of borders
*/
void initialize_borders(double *borders, int n_borders)
{
    // check for the basic scenario
    if (n_borders != 2)
    {
        printf("There are more than 2 borders, borders get initialized to 0.");
        for (int i = 0; i < n_borders; i++)
        {
            borders[i] = 0.0;
        }
    }
    else
    {
        borders[0] = arguments.walkway_width;
        borders[1] = 0.0;
    }
}
