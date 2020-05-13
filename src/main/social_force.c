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

// constants and basic versions
#include "social_force.h"
#include "social_force_model_basic_simplified.h"
#include "social_force_model_basic_simplified_double.h"

// vectorized versions
#include "vectorize/social_force_model_vectorize_1.h"
#include "vectorize/social_force_model_vectorize_2.h"
#include "vectorize/social_force_model_vectorize_3.h"
#include "vectorize/social_force_model_vectorize_4.h"
#include "vectorize/social_force_model_vectorize_1_double.h"
#include "vectorize/social_force_model_vectorize_2_double.h"
#include "vectorize/social_force_model_vectorize_4_double.h"



#include "vectorize/social_force_model_vectorize_2_5_1.h"

// standard C versions
#include "stdc_opt/social_force_model_stdc_optv_2_5_1.h"
#include "stdc_opt/social_force_model_stdc_optv_2_4.h"

// testing
#include "testing_float/test_sets_float.h"
#include "testing_float/compare_simulations_float.h"

#include "testing_double/test_sets_double.h"
#include "testing_double/test_gradients_double.h"
#include "testing_double/compare_simulations_double.h"

// benachmarking
#include "tsc_x86.h"

// utility
#include "utility.h"

// globals
struct arguments arguments = {0};
char filename_global[80];

int seed = 0;

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
    //simulations and tests list
    sim_t *sim_list_float[50];
    sim_t *sim_list_double[50];

    sim_t *test_functions_list[50];

    int sim_counter_float = 0;
    int sim_counter_double = 0;

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

    add_implementations_float(sim_list_float, &sim_counter_float, test_functions_list, &test_func_counter);
    add_implementations_double(sim_list_double, &sim_counter_double, test_functions_list, &test_func_counter);

    if (arguments.test || arguments.visual)
    {
        //seed = j;
        run_tests(sim_list_float, sim_counter_float);
        run_tests_double(sim_list_double, sim_counter_double);

        for (int i = 0; i < test_func_counter; i++)
        {
            if (test_functions_list[i]->is_double)
            {
                sim_func_double current = test_functions_list[i]->f_double;
                run_finite_differences_double(current, arguments);
            }
        }
    }
    else
    {
        /* run benchmarks */

        if (contains_substring(arguments.benchmark, "all"))
        {
            for (int i = 0; i < sim_counter_float; i++)
            {
                sim_t sim = *sim_list_float[i];
                run_bench_float(sim);
            }
            for (int i = 0; i < sim_counter_double; i++)
            {
                sim_t sim = *sim_list_double[i];
                run_bench_double(sim);
            }
        }
        else
        {
            /* run benchmarks only for some of the implementations  */

            for (int i = 0; i < sim_counter_float; i++)
            {
                sim_t sim = *sim_list_float[i];

                if (contains_substring(arguments.benchmark, sim.name))
                {
                    run_bench_float(sim);
                }
            }
            for (int i = 0; i < sim_counter_double; i++)
            {
                sim_t sim = *sim_list_double[i];

                if (contains_substring(arguments.benchmark, sim.name))
                {
                    run_bench_double(sim);
                }
            }
        }
    }
    return 0;
}

/*
*   Function used to add the different implementations of the simulation to benchmark.
*/
void add_implementations_float(sim_t **sim_list, int *sim_counter, sim_t **test_functions_list, int *test_func_counter)
{
    // add_function(sim_list, sim_counter, simulation_basic_simplified, NULL, compute_simplified_flops, IS_FLOAT,compute_operational_intensity_0, "simplified_float");

    add_function(sim_list, sim_counter, simulation_basic_vectorize_1, NULL, compute_simplified_flops, IS_FLOAT, compute_operational_intensity_0, "vectorize_1");
    add_function(sim_list, sim_counter, simulation_basic_vectorize_2, NULL, compute_simplified_flops, IS_FLOAT, compute_operational_intensity_0, "vectorize_2");
    add_function(sim_list, sim_counter, simulation_basic_vectorize_3, NULL, compute_simplified_flops, IS_FLOAT, compute_operational_intensity_0, "vectorize_3");
    add_function(sim_list, sim_counter, simulation_basic_vectorize_4, NULL, compute_simplified_flops, IS_FLOAT, compute_operational_intensity_0, "vectorize_4");
    add_function(sim_list, sim_counter, simulation_basic_vectorize_2_5_1, NULL, compute_simplified_flops, IS_FLOAT,compute_operational_intensity_251, "vectorize_2_5_1");

    add_test_function(test_functions_list, test_simulation_basic_simplified, NULL, IS_FLOAT, test_func_counter);
}

void add_implementations_double(sim_t **sim_list, int *sim_counter, sim_t **test_functions_list, int *test_func_counter)
{
    add_function(sim_list, sim_counter, NULL, simulation_basic_simplified_double, compute_simplified_flops, IS_DOUBLE,compute_operational_intensity_0, "simplified_double");
    add_function(sim_list, sim_counter, NULL, simulation_basic_vectorize_1_double, compute_simplified_flops, IS_DOUBLE, compute_operational_intensity_0, "vectorize_1_double");
    add_function(sim_list, sim_counter, NULL, simulation_basic_vectorize_2_double, compute_simplified_flops, IS_DOUBLE, compute_operational_intensity_0, "vectorize_2_double");
    add_function(sim_list, sim_counter, NULL, simulation_basic_vectorize_4_double, compute_simplified_flops, IS_DOUBLE, compute_operational_intensity_0, "vectorize_4_double");

    add_function(sim_list, sim_counter, NULL, simulation_basic_optv_2_5_1, compute_simplified_flops, IS_DOUBLE,compute_operational_intensity_251, "stdc_optv_2_5_1_double");
    // add_function(sim_list, sim_counter, NULL, simulation_basic_optv_2_4, compute_simplified_flops, IS_DOUBLE, compute_operational_intensity_0, "stdc_optv_2_4_double");

    add_test_function(test_functions_list, NULL, test_simulation_basic_simplified_double, IS_DOUBLE, test_func_counter);
}

/*
*   Appends a simulation function to the list
*/
void add_function(sim_t **sim_list, int *sim_counter, sim_func f, sim_func_double f_double, flops_func flops_f, int is_double, op_int_func op_f, char *name)
{
    (*sim_counter) = *sim_counter + 1;
    sim_t *sim = (sim_t *)malloc(sizeof(sim_t));
    sim->f = f;
    sim->f_double = f_double;
    sim->flops_f = flops_f;
    sim->is_double = is_double;
    sim->op_f = op_f;
    sim->name = name;
    sim_list[*sim_counter - 1] = sim;
}

/*
*   Appends a test simulation function to the list
*/
void add_test_function(sim_t **test_functions_list, sim_func f, sim_func_double f_double, int is_double, int *test_func_counter)
{
    (*test_func_counter) = *test_func_counter + 1;
    sim_t *sim = (sim_t *)malloc(sizeof(sim_t));
    sim->f = f;
    sim->f_double = f_double;
    sim->is_double = is_double;
    test_functions_list[*test_func_counter - 1] = sim;
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
void initialize_people_float(float *position, float *desired_direction, float *final_destination, float *desired_speed, int n)
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

void initialize_people_double(double *position, double *desired_direction, double *final_destination, double *desired_speed, int n)
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
void initialize_max_speed_float(float *desired_speed, float *desired_max_speed, int n)
{
    for (int i = 0; i < n; i++)
    {
        desired_max_speed[i] = desired_speed[i] * 1.3;
    }
}

void initialize_max_speed_double(double *desired_speed, double *desired_max_speed, int n)
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
void initialize_borders_float(float *borders, int n_borders)
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

void initialize_borders_double(double *borders, int n_borders)
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
*   This function benchmarks all the implementation stored into the list.
*/
void run_bench_float(sim_t sim)
{
    char *name = sim.name;
    sim_func f = sim.f;
    flops_func flops_f = sim.flops_f;
    op_int_func op_f = sim.op_f;
    double cycles = 0.;
    double multiplier = 1;
    myInt64 start, end;

    int number_of_people = arguments.n_people;
    int n_timesteps = 1;
    double op_intensity = op_f(number_of_people);
    // allocate memory

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
    float *people_repulsion_term;
    if (contains_substring(name, "simplified") || contains_substring(name, "vectorize_1") || contains_substring(name, "vectorize_2"))
    {
        people_repulsion_term = (float *)aligned_malloc(number_of_people * number_of_people * 2 * sizeof(float), 32);
        set_zero(people_repulsion_term, number_of_people * number_of_people * 2);
    }
    else
    {
        people_repulsion_term = (float *)aligned_malloc(number_of_people * 8 * 2 * sizeof(float), 32);
        set_zero(people_repulsion_term, 8 * number_of_people * 2);
    }
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
    initialize_people_float(position, desired_direction, final_destination, desired_speed, number_of_people);
    initialize_max_speed_float(desired_speed, desired_max_speed, number_of_people);
    initialize_borders_float(borders, N_BORDERS);

    // Warm-up phase
    do
    {
        n_timesteps = n_timesteps * multiplier;
        start = start_tsc();
        f(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed, desired_max_speed);
        end = stop_tsc(start);

        cycles = (double)end;
        // printf("%Lf\t%llu\n", cycles, end);
        multiplier = (CYCLES_REQUIRED) / (cycles);
    } while (multiplier > 2);

    long long unsigned flops = n_timesteps * flops_f(number_of_people);
    double *cycles_list = malloc(sizeof(double) * REP);

    initialize_people_float(position, desired_direction, final_destination, desired_speed, number_of_people);
    initialize_max_speed_float(desired_speed, desired_max_speed, number_of_people);
    initialize_borders_float(borders, N_BORDERS);
    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++)
    {
        start = start_tsc();
        f(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed, desired_max_speed);
        end = stop_tsc(start);

        cycles = ((double)end);
        // printf("%llu\t%Lf\t%d\n", end, cycles, n_timesteps);
        // total_cycles += cycles;

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
    double performance = (double)flops / cycles;
    printf("%s %d %llu %.2f %.3f %.2f\n", name, number_of_people, flops, cycles, performance, op_intensity);
}

/*
*   This function benchmarks all the implementation stored into the list.
*/
void run_bench_double(sim_t sim)
{
    char *name = sim.name;
    sim_func_double f = sim.f_double;
    flops_func flops_f = sim.flops_f;
    op_int_func op_f = sim.op_f;
    char *fun_name = sim.name;
    double cycles = 0.;
    double multiplier = 1;
    myInt64 start, end;

    int number_of_people = arguments.n_people;
    int n_timesteps = 1;
    double op_intensity = op_f(number_of_people);
    // allocate memory

    double *position = (double *)aligned_malloc(number_of_people * 2 * sizeof(double), 32);
    set_zero_double(position, number_of_people * 2);
    double *speed = (double *)aligned_malloc(number_of_people * sizeof(double), 32);
    set_zero_double(speed, number_of_people);
    double *desired_direction = (double *)aligned_malloc(number_of_people * 2 * sizeof(double), 32);
    set_zero_double(desired_direction, number_of_people * 2);
    double *final_destination = (double *)aligned_malloc(number_of_people * 2 * sizeof(double), 32);
    set_zero_double(final_destination, number_of_people * 2);
    double *borders = (double *)aligned_malloc(N_BORDERS * sizeof(double), 32);
    set_zero_double(borders, N_BORDERS);
    double *actual_velocity = (double *)aligned_malloc(number_of_people * 2 * sizeof(double), 32);
    set_zero_double(actual_velocity, number_of_people * 2);
    double *acceleration_term = (double *)aligned_malloc(number_of_people * 2 * sizeof(double), 32);
    set_zero_double(acceleration_term, number_of_people * 2);
    double *people_repulsion_term;
    if (contains_substring(name, "simplified") || contains_substring(name, "vectorize_2") || contains_substring(name, "vectorize_3"))
    {
        people_repulsion_term = (double *)aligned_malloc(number_of_people * number_of_people * 2 * sizeof(double), 32);
        set_zero_double(people_repulsion_term, number_of_people * number_of_people * 2);
    }
    else
    {
        people_repulsion_term = (double *)aligned_malloc(number_of_people * 8 * 2 * sizeof(double), 32);
        set_zero_double(people_repulsion_term, number_of_people * 8 * 2);
    }
    double *border_repulsion_term = (double *)aligned_malloc(number_of_people * N_BORDERS * 2 * sizeof(double), 32);
    set_zero_double(border_repulsion_term, number_of_people * N_BORDERS * 2);
    double *social_force = (double *)aligned_malloc(number_of_people * 2 * sizeof(double), 32);
    set_zero_double(social_force, number_of_people * 2);
    double *desired_speed = (double *)aligned_malloc(number_of_people * sizeof(double), 32);
    set_zero_double(desired_speed, number_of_people);
    double *desired_max_speed = (double *)aligned_malloc(number_of_people * sizeof(double), 32);
    set_zero_double(desired_max_speed, number_of_people);
    // check if calloc worked correctly
    if (position == NULL || speed == NULL || desired_direction == NULL || final_destination == NULL || borders == NULL || actual_velocity == NULL || acceleration_term == NULL || people_repulsion_term == NULL || border_repulsion_term == NULL || social_force == NULL || desired_speed == NULL || desired_max_speed == NULL)
    {
        printf("Error: calloc failed\n");
        return;
    }

    // initialize arrays
    initialize_people_double(position, desired_direction, final_destination, desired_speed, number_of_people);
    initialize_max_speed_double(desired_speed, desired_max_speed, number_of_people);
    initialize_borders_double(borders, N_BORDERS);

    // Warm-up phase
    do
    {
        n_timesteps = n_timesteps * multiplier;
        start = start_tsc();
        f(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed, desired_max_speed);
        end = stop_tsc(start);

        cycles = (double)end;
        // printf("%Lf\t%llu\n", cycles, end);
        multiplier = (CYCLES_REQUIRED) / (cycles);
    } while (multiplier > 2);

    long long unsigned flops = n_timesteps * flops_f(number_of_people);
    double *cycles_list = malloc(sizeof(double) * REP);

    initialize_people_double(position, desired_direction, final_destination, desired_speed, number_of_people);
    initialize_max_speed_double(desired_speed, desired_max_speed, number_of_people);
    initialize_borders_double(borders, N_BORDERS);
    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++)
    {

        start = start_tsc();
        f(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed, desired_max_speed);
        end = stop_tsc(start);

        cycles = ((double)end);
        // printf("%llu\t%Lf\t%d\n", end, cycles, n_timesteps);
        // total_cycles += cycles;

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
    double performance = (double)flops / cycles;
    printf("%s %d %llu %.2f %.3f %.2f\n", name, number_of_people, flops, cycles, performance, op_intensity);
}

/*
*   Function used to sort the cycles.
*/
int compare(const void *a, const void *b)
{
    return (*(double *)a - *(double *)b);
}

double compute_operational_intensity_251(int number_of_people)
{
    return (52 * (double)number_of_people + 25) / 104;
}

double compute_operational_intensity_0(int n)
{
    return 0.0;
}

/*
*   Function that returns the number of flops. Computed useing wxMaxima.
*/
long long unsigned compute_basic_flops(int number_of_people)
{
    long long unsigned n = number_of_people;
    long long unsigned nb = N_BORDERS;
    long long unsigned a, b, c, d, e, f;
    a = add_cost;
    b = mult_cost;
    c = div_cost;
    d = sqrt_cost;
    e = exp_cost;
    f = fabs_cost;
    long long unsigned flops = n * (3 * a + 2 * b + 2 * c * d) +
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
long long unsigned compute_simplified_flops(int number_of_people)
{
    long long unsigned n = number_of_people;
    long long unsigned nb = N_BORDERS;
    long long unsigned a, b, c, d, fe;
    a = add_cost;
    b = mult_cost;
    c = div_cost;
    d = sqrt_cost;
    fe = fast_exp_cost;
    
    long long unsigned flops = n * (3 * a + 2 * b + 2 * c + d) +
                               n * (2 * a + 4 * b) +
                               (n * n - n) * (12 * a + 20 * b + 7 * c + 4 * d + fe) +
                               (nb * n) * (a + 3 * b + 3 * c + fe) +
                               n * (n + nb) * 2 * a +
                               n * (5 * a + 9 * b + 3 * c + d);
    return flops;
}
