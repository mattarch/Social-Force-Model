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
#include <string.h>

// parsing arguments
#include "parse_args.h"

// constants
#include "social_force.h"
#include "social_force_model_basic.h"

// testing
#include "testing.h"

// benachmarking
#include "tsc_x86.h"

// utility
#include "utility.h"

// globals
struct arguments arguments = {0};
char filename_global[40];

//cost of each operation
const int add_cost = 1;
const int mult_cost = 1;
const int div_cost = 1;
const int exp_cost = 1;
const int sqrt_cost = 1;
const int fabs_cost = 1;

// main function
int main(int argc, char *argv[])
{

    sim_t *sim_list[50];
    //simulations and tests list
    sim_func test_functions_list[50];

    int sim_counter = 0;

    int test_func_counter = 0;

    // default values for arguments
    arguments.n_people = 300;
    arguments.n_timesteps = 300;
    arguments.test = false;
    arguments.walkway_width = 4;
    arguments.walkway_length = 50;
    arguments.benchmark = "all";
    arguments.visual = false;
    arguments.filename = "\n";

    // parse arguments
    parse_args(argc, argv, &arguments);

    add_implementations(sim_list, &sim_counter, test_functions_list, &test_func_counter);
    if (arguments.test || arguments.visual)
    {
        run_tests(sim_list, sim_counter);
        for (int i = 0; i < test_func_counter; i++)
        {
            sim_func current = test_functions_list[i];
            run_sim_test(current,arguments);
        }
    }
    else
    {
        if (strstr(arguments.benchmark, "all") != 0)
        {
            for (int i = 0; i < sim_counter; i++)
            {
                sim_t sim = *sim_list[i];
                run_bench(sim);
            }
        }
        else
        {
            for (int i = 0; i < sim_counter; i++)
            {
                sim_t sim = *sim_list[i];
                if (strstr(arguments.benchmark, sim.name) != 0)
                {
                    run_bench(sim);
                }
            }
        }
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

/*
*   Function used to add the different implementations of the simulation to benchmark.
*/
void add_implementations(sim_t **sim_list, int *sim_counter, sim_func *test_functions_list, int *test_func_counter)
{
    add_function(sim_list, sim_counter, simulation_basic, "basic");
    add_test_function(test_functions_list, test_simulation_basic, test_func_counter);
}

/*
*   Function used to sort the cycles.
*/
int compare(const void *a, const void *b)
{
    return (*(double *)a - *(double *)b);
}

/*
*   This function benchmarks all the implementation stored into the list.
*/
void run_bench(sim_t sim)
{
    char *name = sim.name;
    sim_func f = sim.f;
    double cycles = 0.;
    long num_runs = 5;
    double multiplier = 1;
    myInt64 start, end;

    int number_of_people = arguments.n_people;
    int n_timesteps = arguments.n_timesteps;
    // allocate memory

    long long unsigned int flops = arguments.n_timesteps * compute_flops(number_of_people);

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
        return;
    }

    // initialize arrays
    initialize_people(position, desired_direction, final_destination, desired_speed, number_of_people);
    initialize_borders(borders, N_BORDERS);

    // Warm-up phase
    do
    {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++)
        {
            f(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    double *cycles_list = malloc(sizeof(double) * REP);

    initialize_people(position, desired_direction, final_destination, desired_speed, number_of_people);
    initialize_borders(borders, N_BORDERS);
    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++)
    {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i)
        {
            f(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycles_list[j] = cycles;
    }
    total_cycles /= REP;

    free(position);
    free(speed);
    free(desired_direction);
    free(final_destination);
    free(borders);
    free(actual_velocity);
    free(acceleration_term);
    free(people_repulsion_term);
    free(border_repulsion_term);
    free(social_force);
    free(desired_speed);
    qsort(cycles_list, REP, sizeof(double), compare);
    cycles = cycles_list[REP / 2]; //total_cycles;
    free(cycles_list);

    printf("%s %d %llu %f %.8f\n", name, number_of_people, flops, cycles, flops / cycles);
}



/*
*   Function that returns the number of flops. Computed used wxMaxima.
*/
long long unsigned int compute_flops(int number_of_people)
{
    return number_of_people * (fabs_cost * N_BORDERS + 2 * sqrt_cost * N_BORDERS +
                               6 * div_cost * N_BORDERS + 4 * mult_cost * N_BORDERS + 3 * add_cost * N_BORDERS +
                               2 * exp_cost * number_of_people + 2 * sqrt_cost * number_of_people +
                               13 * div_cost * number_of_people + 15 * mult_cost * number_of_people +
                               16 * add_cost * number_of_people - 2 * exp_cost - sqrt_cost - 9 * div_cost - 7 * mult_cost - 9 * add_cost);
}

/*
*   Appends a simulation function to the list
*/
void add_function(sim_t **sim_list, int *sim_counter, sim_func f, char *name)
{
    (*sim_counter) = *sim_counter + 1;
    sim_t *sim = (sim_t *)malloc(sizeof(sim_t));
    sim->f = f;
    sim->name = name;
    sim_list[*sim_counter - 1] = sim;
}

/*
*   Appends a test simulation function to the list
*/
void add_test_function(sim_func *test_functions_list, sim_func f, int *test_func_counter)
{
    (*test_func_counter) = *test_func_counter + 1;
    test_functions_list[*test_func_counter - 1] = f;
}