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

// fix aligenment
#include "aligned_free.h"
#include "aligned_malloc.h"

// parsing arguments
#include "parse_args.h"

// constants
#include "social_force.h"
#include "social_force_model_basic.h"
#include "social_force_model_basic_simplified.h"

// vectorized versions
#include "vectorize/social_force_model_vectorize_1.h"
#include "vectorize/social_force_model_vectorize_2.h"
#include "vectorize/social_force_model_vectorize_3.h"
#include "vectorize/social_force_model_vectorize_4.h"
#include "vectorize/social_force_model_vectorize_5.h"

#include "vectorize/social_force_model_vectorize_2_5_1.h"

// testing
#include "testing.h"

// benachmarking
#include "tsc_x86.h"

// utility
#include "utility.h"

// globals
struct arguments arguments = {0};
char filename_global[80];

int seed = 0;

#define IS_FLOAT 0
#define IS_DOUBLE 1

//cost of each operation
const int add_cost = 1;
const int mult_cost = 1;
const int div_cost = 1;
const int exp_cost = 1;
const int sqrt_cost = 1;
const int fabs_cost = 1;
const int fast_exp_cost = 16;

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
        for (int j = 0; j < 1; j++)
        {
            seed = j;
            run_tests(sim_list, sim_counter);
        }
        for (int i = 0; i < test_func_counter; i++)
        {
            sim_func current = test_functions_list[i];
            run_finite_differences(current, arguments);
        }
    }
    else
    {
        if (contains_substring(arguments.benchmark, "all"))
        {
            for (int i = 0; i < sim_counter; i++)
            {

                sim_t sim = *sim_list[i];
                if (sim.is_double)
                {
                    run_bench_double(sim);
                }
                else
                {
                    run_bench_float(sim);
                }
            }
        }
        else
        {
            for (int i = 0; i < sim_counter; i++)
            {
                sim_t sim = *sim_list[i];

                if (contains_substring(arguments.benchmark, sim.name))
                {
                    if (sim.is_double)
                    {
                        run_bench_double(sim);
                    }
                    else
                    {
                        run_bench_float(sim);
                    }
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
void initialize_people(float *position, float *desired_direction, float *final_destination, float *desired_speed, int n)
{

    srand(seed);

    for (int i = 0; i < n; i++)
    {
        // initialize values independant of starting point and target point
        double current_y = rand() * arguments.walkway_width / RAND_MAX; // starting position y coordinate
        double target_y = rand() * arguments.walkway_width / RAND_MAX;  // target y coordinate
        desired_speed[i] = sampleNormal(0.0676, AVG_SPEED);
        double current_x, target_x;
        if (i % 2) // initialize this person to walk from left to right
        {
            current_x = 0.0 - rand() * arguments.walkway_length / RAND_MAX;                                 // starting position x coordinate
            target_x = arguments.walkway_length + 10 + rand() * (arguments.walkway_length - 10) / RAND_MAX; // target x coordinate
        }
        else // initialize this person to walk from right to left
        {
            current_x = arguments.walkway_length + rand() * arguments.walkway_length / RAND_MAX; // starting position x coordinate
            target_x = -10 - rand() * (arguments.walkway_length - 10) / RAND_MAX;                // target x coordinate
        }

        // compute differences
        double delta_x = target_x - current_x; // 1 add => 1 flop
        double delta_y = target_y - current_y; // 1 add => 1 flop

        // normalization constant
        double d = delta_x * delta_x + delta_y * delta_y; // 1 add, 2 mult => 3 flops
        double normalizer = sqrt(d);                      // 1 sqrt => 1 flop

        // update desired_direction
        desired_direction[IndexX(i)] = delta_x / normalizer;    // 1 div => 1 flop
        desired_direction[IndexY(i, n)] = delta_y / normalizer; // 1 div => 1 flop

        // write position and target to vector
        position[IndexX(i)] = current_x;
        position[IndexY(i, n)] = current_y;
        final_destination[IndexX(i)] = target_x;
        final_destination[IndexY(i, n)] = target_y;
    }
}
/*
  This function computes the max speed value for every person.

  Assumptions: There is at least one person.
               
  Parameters: desired_speed (n,1) = vector of desired speed values
              desired_max_speed (n,1) = vector of max speed values
                     n: number of people
*/
void compute_max_speed(float *desired_speed, float *desired_max_speed, int n)
{
    for (int i = 0; i < n; i++)
    {
        desired_max_speed[i] = desired_speed[i] * 1.3;
    }
}

/*
  This function initializes the borders array with reasonable starting values.

  Assumptions: sidewalk scenario --> two horizontal borders, bottom one on height 0 top one on height arguments.walkway_width.
  Parameters:
               borders: (n_borders,2) : array of borders
             n_borders: number of borders
*/
void initialize_borders(float *borders, int n_borders)
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
    //add_function(sim_list, sim_counter, simulation_basic, compute_basic_flops, "basic");
    // add_function(sim_list, sim_counter, simulation_basic_simplified, compute_simplified_flops, "simplified");

    //add_function(sim_list, sim_counter, simulation_basic_vectorize_1, compute_simplified_flops, "vectorize_1");
    //add_function(sim_list, sim_counter, simulation_basic_vectorize_2_5_1, compute_simplified_flops, "vectorize_2_5_1");
    //add_function(sim_list, sim_counter, simulation_basic_vectorize_2, compute_simplified_flops, "vectorize_2");
    //add_function(sim_list, sim_counter, simulation_basic_vectorize_3, compute_simplified_flops, "vectorize_3");
    add_function(sim_list, sim_counter, simulation_basic_vectorize_4, compute_simplified_flops, IS_FLOAT, "vectorize_4");
    add_function(sim_list, sim_counter, simulation_basic_vectorize_5, compute_simplified_flops, IS_FLOAT, "vectorize_5");

    //add_test_function(test_functions_list, test_simulation_basic, test_func_counter);
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
void run_bench_float(sim_t sim)
{
    char *name = sim.name;
    sim_func f = sim.f;
    flops_func flops_f = sim.flops_f;
    double cycles = 0.;
    long num_runs = 5;
    double multiplier = 1;
    myInt64 start, end;

    int number_of_people = arguments.n_people;
    int n_timesteps = arguments.n_timesteps;
    // allocate memory

    long long unsigned int flops = arguments.n_timesteps * flops_f(number_of_people);

    float *position = (float *)aligned_malloc(number_of_people * 2 * sizeof(float), 32);
    set_zero(position, number_of_people * 2);
    float *speed = (float *)aligned_malloc(number_of_people * sizeof(float), 32);
    set_zero(speed, number_of_people);
    float *desired_direction = (float *)aligned_malloc(number_of_people * 2 * sizeof(float), 32);
    set_zero(desired_direction, number_of_people * 2);
    float *final_destination = (float *)aligned_malloc(number_of_people * 2 * sizeof(float), 32);
    set_zero(final_destination, number_of_people * 2);
    float *borders = (float *)aligned_malloc(N_BORDERS * sizeof(float), 32);
    set_zero(borders, N_BORDERS);
    float *actual_velocity = (float *)aligned_malloc(number_of_people * 2 * sizeof(float), 32);
    set_zero(actual_velocity, number_of_people * 2);
    float *acceleration_term = (float *)aligned_malloc(number_of_people * 2 * sizeof(float), 32);
    set_zero(acceleration_term, number_of_people * 2);
    float *people_repulsion_term = (float *)aligned_malloc(number_of_people * 8 * 2 * sizeof(float), 32);
    set_zero(people_repulsion_term, 8 * number_of_people * 2);
    float *border_repulsion_term = (float *)aligned_malloc(number_of_people * N_BORDERS * 2 * sizeof(float), 32);
    set_zero(border_repulsion_term, number_of_people * N_BORDERS * 2);
    float *social_force = (float *)aligned_malloc(number_of_people * 2 * sizeof(float), 32);
    set_zero(social_force, number_of_people * 2);
    float *desired_speed = (float *)aligned_malloc(number_of_people * sizeof(float), 32);
    set_zero(desired_speed, number_of_people);
    float *desired_max_speed = (float *)aligned_malloc(number_of_people * sizeof(float), 32);
    set_zero(desired_max_speed, number_of_people);
    // check if calloc worked correctly
    if (position == NULL || speed == NULL || desired_direction == NULL || final_destination == NULL || borders == NULL || actual_velocity == NULL || acceleration_term == NULL || people_repulsion_term == NULL || border_repulsion_term == NULL || social_force == NULL || desired_speed == NULL || desired_max_speed == NULL)
    {
        printf("Error: calloc failed\n");
        return;
    }

    // initialize arrays
    initialize_people(position, desired_direction, final_destination, desired_speed, number_of_people);
    compute_max_speed(desired_speed, desired_max_speed, number_of_people);
    initialize_borders(borders, N_BORDERS);

    // Warm-up phase
    do
    {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++)
        {
            f(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed, desired_max_speed);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    double *cycles_list = malloc(sizeof(double) * REP);

    initialize_people(position, desired_direction, final_destination, desired_speed, number_of_people);
    compute_max_speed(desired_speed, desired_max_speed, number_of_people);
    initialize_borders(borders, N_BORDERS);
    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++)
    {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i)
        {
            f(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed, desired_max_speed);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycles_list[j] = cycles;
    }
    total_cycles /= REP;

    aligned_free(position);
    aligned_free(speed);
    aligned_free(desired_direction);
    aligned_free(final_destination);
    aligned_free(borders);
    aligned_free(actual_velocity);
    aligned_free(acceleration_term);
    aligned_free(people_repulsion_term);
    aligned_free(border_repulsion_term);
    aligned_free(social_force);
    aligned_free(desired_speed);
    aligned_free(desired_max_speed);
    qsort(cycles_list, REP, sizeof(double), compare);
    cycles = cycles_list[REP / 2]; //total_cycles;
    free(cycles_list);

    printf("%s %d %llu %f %.8f\n", name, number_of_people, flops, cycles, flops / cycles);
}

/*
*   This function benchmarks all the implementation stored into the list.
*/
void run_bench_double(sim_t sim)
{
}

/*
*   Function that returns the number of flops. Computed useing wxMaxima.
*/
long long unsigned int compute_basic_flops(int number_of_people)
{
    int n = number_of_people;
    int nb = N_BORDERS;
    int a, b, c, d, e, f;
    a = add_cost;
    b = mult_cost;
    c = div_cost;
    d = sqrt_cost;
    e = exp_cost;
    f = fabs_cost;
    long long unsigned int flops = n * (3 * a + 2 * b + 2 * c * d) +
                                   n * 2 * b +
                                   n * (2 * a + 4 * b + 2 * c) +
                                   (n * n - n) * (15 * a + 22 * b + 13 * c + 4 * d + 2 * e) +
                                   (nb * n) * (a + 4 * b + 6 * c + 2 * e + f) +
                                   n * (n + nb) * 2 * a +
                                   n * (a * 6 + 12 * b + 3 * c + 2 * d);
    return flops;
}

/*
*   Function that returns the number of flops. Computed useing wxMaxima.
*/
long long unsigned int compute_simplified_flops(int number_of_people)
{
    int n = number_of_people;
    int nb = N_BORDERS;
    int a, b, c, d, fe;
    a = add_cost;
    b = mult_cost;
    c = div_cost;
    d = sqrt_cost;
    fe = fast_exp_cost;
    long long unsigned int flops = n * (3 * a + 2 * b + 2 * c * d) +
                                   n * (2 * a + 4 * b) +
                                   (n * n - n) * (12 * a + 20 * b + 7 * c + 4 * d + fe) +
                                   (nb * n) * (a + 3 * b + 3 * c + fe) +
                                   n * (n + nb) * 2 * a +
                                   n * (5 * a + 9 * b + 3 * c + d);
    return flops;
}

/*
*   Appends a simulation function to the list
*/
void add_function(sim_t **sim_list, int *sim_counter, sim_func f, flops_func flops_f, int is_double, char *name)
{
    (*sim_counter) = *sim_counter + 1;
    sim_t *sim = (sim_t *)malloc(sizeof(sim_t));
    sim->f = f;
    sim->flops_f = flops_f;
    sim->is_double = is_double;
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