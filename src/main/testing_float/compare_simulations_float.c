/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// fix aligenment
#include "aligned_free.h"
#include "aligned_malloc.h"

#include "parse_args.h"
#include "social_force.h"
#include "compare_simulations_float.h"
#include "test_sets_float.h"

#include "utility.h"

extern int seed;

void print_errors_to_file(int number_of_people, float *error_vector)
{
    FILE *filePtr;
    char file_[50] = "error_vector_seed_0_simplified-v1-15-";
    char snum[10];
    sprintf(snum, "%d.txt", seed);
    strcat(file_, snum);

    filePtr = fopen(file_, "w");

    for (int i = 0; i < number_of_people; i++)
    {
        fprintf(filePtr, "%lf\n", error_vector[i]);
    }
    fclose(filePtr);
}

/*
*   Function used to check the square distance between the expected result of a testcase
*   and the result of the function testes.
*   This function returns 1 if the distance is greater than the threshold, 0 otherwise.
*/
int check_absolute_distance_vec(float *expected, float *res, float *error_vector, int n)
{
    int check = 0;
    float err;
    for (int i = 0; i < n; i++)
    {
        // err = (expected[i] - res[i]) * (expected[i] - res[i]);
        err = fabs(expected[i] - res[i]);
        error_vector[i] = err;
        if (err > EPS || isnan(err))
        {
            //printf("expected: %f, got %f, abs error %f\n", expected[i], res[i], err);
            check = 1;
        }
    }

    return check;
}


/*
* Runs all the simulation and compare the result to the basic version
* Returns 0 if OK, returns 1 if there's an error.
*/
int compare_simulations(sim_t **sim_list, int sim_counter)
{
    int error_check = 0;
    int number_of_people = 80;
    int n_timesteps = 1;
    int n_test_timesteps = 300;
    float *oracle_position, *oracle_speed, *oracle_desired_direction, *oracle_final_destination,
        *oracle_borders, *oracle_actual_velocity, *oracle_acceleration_term, *oracle_people_repulsion_term,
        *oracle_border_repulsion_term, *oracle_social_force, *oracle_desired_speed, *oracle_desired_max_speed;
    float *current_position, *current_speed, *current_desired_direction, *current_final_destination,
        *current_borders, *current_actual_velocity, *current_acceleration_term, *current_people_repulsion_term,
        *current_border_repulsion_term, *current_social_force, *current_desired_speed, *current_desired_max_speed;

    float *error_vector = (float *)malloc(number_of_people * 2 * sizeof(float));

    float *starting_position = (float *)aligned_malloc(number_of_people * 2 * sizeof(float), 32);
    set_zero(starting_position, number_of_people * 2);
    float *starting_desired_direction = (float *)aligned_malloc(number_of_people * 2 * sizeof(float), 32);
    set_zero(starting_desired_direction, number_of_people * 2);
    float *starting_final_destination = (float *)aligned_malloc(number_of_people * 2 * sizeof(float), 32);
    set_zero(starting_final_destination, number_of_people * 2);
    float *starting_borders = (float *)aligned_malloc(N_BORDERS * sizeof(float), 32);
    set_zero(starting_borders, N_BORDERS);
    float *starting_desired_speed = (float *)aligned_malloc(number_of_people * sizeof(float), 32);
    set_zero(starting_desired_speed, number_of_people);
    float *starting_desired_max_speed = (float *)aligned_malloc(number_of_people * sizeof(float), 32);
    set_zero(starting_desired_max_speed, number_of_people);
    // // check if calloc worked correctly
    if (starting_position == NULL || starting_desired_direction == NULL || starting_final_destination == NULL || starting_borders == NULL || starting_desired_speed == NULL || starting_desired_max_speed == NULL)
    {
        printf("Error: calloc failed\n");
        return 1;
    }

    // initialize starting position
    initialize_people_float(starting_position, starting_desired_direction, starting_final_destination, starting_desired_speed, number_of_people);
    initialize_max_speed_float(starting_desired_speed, starting_desired_max_speed, number_of_people);
    initialize_borders_float(starting_borders, N_BORDERS);
    for (int i = 1; i < sim_counter; i++)
    {

        copy_init(starting_position, starting_desired_direction, starting_final_destination, starting_borders, starting_desired_speed, starting_desired_max_speed,
                  &oracle_position, &oracle_desired_direction, &oracle_final_destination, &oracle_borders, &oracle_desired_speed, &oracle_desired_max_speed, number_of_people);
        allocate_arrays(&oracle_speed, &oracle_actual_velocity, &oracle_acceleration_term, &oracle_people_repulsion_term,
                        &oracle_border_repulsion_term, &oracle_social_force, number_of_people);

        copy_init(starting_position, starting_desired_direction, starting_final_destination, starting_borders, starting_desired_speed, starting_desired_max_speed,
                  &current_position, &current_desired_direction, &current_final_destination, &current_borders, &current_desired_speed, &current_desired_max_speed, number_of_people);
        allocate_arrays(&current_speed, &current_actual_velocity, &current_acceleration_term, &current_people_repulsion_term,
                        &current_border_repulsion_term, &current_social_force, number_of_people);

        for (int j = 0; j < n_test_timesteps; j++)
        {
            error_check = 0;
            sim_func f = sim_list[0]->f;
            f(number_of_people, n_timesteps, oracle_position, oracle_speed, oracle_desired_direction, oracle_final_destination,
              oracle_borders, oracle_actual_velocity, oracle_acceleration_term,
              oracle_people_repulsion_term, oracle_border_repulsion_term, oracle_social_force, oracle_desired_speed, oracle_desired_max_speed);
            //test all implementations

            int check = 0;

            f = sim_list[i]->f;
            f(number_of_people, n_timesteps, current_position, current_speed, current_desired_direction, current_final_destination,
              current_borders, current_actual_velocity, current_acceleration_term,
              current_people_repulsion_term, current_border_repulsion_term, current_social_force, current_desired_speed, current_desired_max_speed);

            check += check_absolute_distance(oracle_position, current_position, number_of_people, 1);
            //check += check_absolute_distance(oracle_speed, current_speed, number_of_people, 0);
            //check += check_absolute_distance(oracle_desired_direction, current_desired_direction, number_of_people, 1);
            //check += check_absolute_distance(oracle_acceleration_term, current_acceleration_term, number_of_people, 1);
            //check += check_absolute_distance(oracle_people_repulsion_term, current_people_repulsion_term, number_of_people, 2);

            //check += check_absolute_distance(oracle_border_repulsion_term, current_border_repulsion_term, number_of_people, 3);

            if (check)
            {
                printf("ERROR: implementation %s differs in iteration %d from the base\n", sim_list[i]->name, j);
                error_check = 1;
            }
            else
            {
                printf("%s iteration %d CORRECT!\n", sim_list[i]->name, j);
            }

            copy_state(oracle_position, oracle_desired_direction, oracle_final_destination, oracle_borders, oracle_desired_speed, oracle_desired_max_speed,
                       &current_position, &current_desired_direction, &current_final_destination, &current_borders, &current_desired_speed, &current_desired_max_speed, number_of_people);
        }

        //check_absolute_distance_vec(oracle_position, current_position, error_vector, 2 * number_of_people);

        //print_errors_to_file(2 * number_of_people, error_vector);

        free_all(11, &current_position, &current_speed, &current_desired_direction, &current_final_destination, &current_borders, &current_actual_velocity,
                 &current_acceleration_term, &current_people_repulsion_term, &current_border_repulsion_term, &current_social_force, &current_desired_speed, &current_desired_max_speed);

        free_all(11, &oracle_position, &oracle_speed, &oracle_desired_direction, &oracle_final_destination, &oracle_borders, &oracle_actual_velocity,
                 &oracle_acceleration_term, &oracle_people_repulsion_term, &oracle_border_repulsion_term, &oracle_social_force, &oracle_desired_speed, &oracle_desired_max_speed);
        if (error_check)
        {
            printf("ERROR: implementation %s differs from the base\n", sim_list[i]->name);
        }
    }

    return error_check;
}

/*
*   Function used to check the absolute distance between the expected result of the simplified version
*   and the result of the function tested.
*   This function returns 1 if any distance is greater than the threshold, 0 otherwise.
*/
int check_absolute_distance(float *expected, float *res, int n, int case_n)
{
    int wrong = 0;
    float acc = 0.0;
    if (case_n == 0)
    {
        for (int i = 0; i < n; i++)
        {
            acc = fabs(expected[i] - res[i]);
            if (isnan(acc) || acc > EPS_NEW)
            {
                wrong = 1;
                printf("%lf %lf %lf %d \n", expected[i], res[i], acc, i);
            }
        }
    }
    else if (case_n == 1)
    {
        for (int i = 0; i < n; i++)
        {
            acc = fabs(expected[IndexX(i)] - res[IndexX(i)]);

            if (isnan(acc) || acc > EPS_NEW)
            {
                wrong = 1;
                printf("x: %lf %lf %lf %d \n", expected[IndexX(i)], res[IndexX(i)], acc, IndexX(i));
            }
        }

        for (int i = 0; i < n; i++)
        {
            acc = fabs(expected[IndexY(i, n)] - res[IndexY(i, n)]);

            if (isnan(acc) || acc > EPS_NEW)
            {
                wrong = 1;
                printf("y: %lf %lf %lf %d \n", expected[IndexY(i, n)], res[IndexY(i, n)], acc, IndexY(i, n));
            }
        }
    }
    else if (case_n == 2)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {

                acc = fabs(expected[IndexX_matrix(i, j, n)] - res[IndexX_matrix(i, j, n)]);

                if (isnan(acc) || acc > EPS_NEW)
                {
                    wrong = 1;
                    printf("x: %lf %lf %d %d\n", expected[IndexX_matrix(i, j, n)], res[IndexX_matrix(i, j, n)], i, j);
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                acc = fabs(expected[IndexY_matrix(i, j, n)] - res[IndexY_matrix(i, j, n)]);
                if (isnan(acc) || acc > EPS_NEW)
                {
                    wrong = 1;
                    printf("y: %lf %lf %d %d\n", expected[IndexY_matrix(i, j, n)], res[IndexY_matrix(i, j, n)], i, j);
                }
            }
        }
    }
    else if (case_n == 3)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 2; j++)
            {

                acc = fabs(expected[IndexX_border(i, j, n)] - res[IndexX_border(i, j, n)]);
                if (isnan(acc) || acc > EPS_NEW)
                {
                    wrong = 1;
                    printf("x: %lf %lf %d %d\n", expected[IndexX_border(i, j, n)], res[IndexX_border(i, j, n)], i, j);
                }
            }
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                acc = fabs(expected[IndexY_border(i, j, n)] - res[IndexY_border(i, j, n)]);
                if (isnan(acc) || acc > EPS_NEW)
                {
                    wrong = 1;
                    printf("y: %lf %lf %d %d\n", expected[IndexY_border(i, j, n)], res[IndexY_border(i, j, n)], i, j);
                }
            }
        }
    }
    return wrong;
}

void copy_init(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
               float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n)
{
    *pos = (float *)aligned_malloc(n * 2 * sizeof(float), 32);
    *dir = (float *)aligned_malloc(n * 2 * sizeof(float), 32);
    *fdes = (float *)aligned_malloc(n * 2 * sizeof(float), 32);
    *bor = (float *)aligned_malloc(N_BORDERS * sizeof(float), 32);
    *spe = (float *)aligned_malloc(n * sizeof(float), 32);
    *mspe = (float *)aligned_malloc(n * sizeof(float), 32);

    memcpy(*pos, s_pos, n * 2 * sizeof(float));   //
    memcpy(*dir, s_dir, n * 2 * sizeof(float));   //
    memcpy(*fdes, s_fdes, n * 2 * sizeof(float)); //
    memcpy(*bor, s_bor, N_BORDERS * sizeof(float));
    memcpy(*spe, s_spe, n * sizeof(float));
    memcpy(*mspe, s_mspe, n * sizeof(float));
}

void copy_state(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
                float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n)
{
    memcpy(*pos, s_pos, n * 2 * sizeof(float));   //
    memcpy(*dir, s_dir, n * 2 * sizeof(float));   //
    memcpy(*fdes, s_fdes, n * 2 * sizeof(float)); //
    memcpy(*bor, s_bor, N_BORDERS * sizeof(float));
    memcpy(*spe, s_spe, n * sizeof(float));
    memcpy(*mspe, s_mspe, n * sizeof(float));
}

void allocate_arrays(float **spe, float **vel, float **acc, float **prep, float **brep,
                     float **frc, int n)
{
    *spe = (float *)aligned_malloc(n * sizeof(float), 32);
    set_zero(*spe, n);
    *vel = (float *)aligned_malloc(n * 2 * sizeof(float), 32);
    set_zero(*vel, n * 2);
    *acc = (float *)aligned_malloc(n * 2 * sizeof(float), 32);
    set_zero(*acc, n * 2);
    *prep = (float *)aligned_malloc(n * n * 2 * sizeof(float), 32);
    set_zero(*prep, n * n * 2);
    *brep = (float *)aligned_malloc(n * N_BORDERS * 2 * sizeof(float), 32);
    set_zero(*brep, n * N_BORDERS * 2);
    *frc = (float *)aligned_malloc(n * 2 * sizeof(float), 32);
    set_zero(*frc, n * 2);
}
