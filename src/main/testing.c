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
#include "testing.h"
#include "test_sets.h"
#include "social_force_model_basic.h"
#include "utility.h"

testcase_t **testcases[N_TESTS];

int counter[2 * N_TESTS];

void (**direction_ptr_v)(float *, float *, float *, int);
void (**acceleration_ptr_v)(float *, float *, float *, float *, int);
void (**social_ptr_v)(float *, float *, float *, float *, int, int);
void (**pos_ptr_v)(float *, float *, float *, float *, float *, float *, int);

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
int check_square_distance_vec(float *expected, float *res, float *error_vector, int n)
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
*   Function that adds all the testcases to try to the list
*/
void add_testcases()
{
    //direction testcases
    add_direction_testcase("basic direction test", direction_position0, direction_fdest0, direction_expected0, direction_n0);

    //acceleration testcases
    add_acceleration_testcase("acceleration_test_from_start_straight_and_with_angle", acceleration_direction0, acceleration_vel0, acceleration_desired_speed, acceleration_expected0, acceleration_n0);
    add_acceleration_testcase("acceleration_test_with_velocity_straight_and_with_angle", acceleration_direction1, acceleration_vel1, acceleration_desired_speed, acceleration_expected1, acceleration_n1);
    add_acceleration_testcase("acceleration_test_to_0_straight_and_with_angle", acceleration_direction2, acceleration_vel2, acceleration_desired_speed, acceleration_expected2, acceleration_n2);
    add_acceleration_testcase("deacceleration_test_straight_and_with_angle", acceleration_direction3, acceleration_vel3, acceleration_desired_speed, acceleration_expected3, acceleration_n3);

    //social force
    add_compute_social_force_testcase("basic test", social_acc0, social_prep0, social_brep0, social_expected0, social_n0, social_nb0);

    //position
    add_position_testcase("position_test_from_origin_0_speed_unit_social_force", position_pos0, position_dir, position_speed, position_force0, position_vel0, position_desired_speed, position_expected0, position_n0);
    add_position_testcase("position_test_from_origin_MAX_speed_unit_social_force", position_pos1, position_dir, position_speed, position_force1, position_vel1, position_desired_speed, position_expected1, position_n1);
}

/*
*   function that adds different implementations to test for correctness
*/
void add_function_implementations()
{
    //update direction implementations
    add_direction_implementation(update_desired_direction);

    //update acceleration implementations
    add_acceleration_implementation(update_acceleration_term);

    //compute social force implementations
    add_social_implementation(compute_social_force);

    //update position implementations
    add_pos_implementation(update_position);
}

int run_tests(sim_t **sim_list, int sim_counter)
{
    //init auxilliary parameters
    for (int i = 0; i < N_TESTS * 2; i++)
    {
        counter[i] = 0;
    }

    //add tests
    add_testcases();
    add_function_implementations();

    //run tests
    int error = run_testcases();
    int errorsim = compare_simulations(sim_list, sim_counter);
    return error + errorsim;
}

/*
* Runs all the function in the test list. These functions are used to check that the
* gradient computed analytically are correct by comparing the result to a finite differences
* implementation of the gradients.
*/
void run_finite_differences(sim_func f, struct arguments arguments)
{
    int number_of_people = arguments.n_people;
    int n_timesteps = arguments.n_timesteps;
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
    float *people_repulsion_term = (float *)aligned_malloc(number_of_people * number_of_people * 2 * sizeof(float), 32);
    set_zero(people_repulsion_term, number_of_people * number_of_people * 2);
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
    f(number_of_people, n_timesteps, position, speed, desired_direction, final_destination, borders, actual_velocity, acceleration_term, people_repulsion_term, border_repulsion_term, social_force, desired_speed, desired_max_speed);
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
    initialize_people(starting_position, starting_desired_direction, starting_final_destination, starting_desired_speed, number_of_people);
    compute_max_speed(starting_desired_speed, starting_desired_max_speed, number_of_people);
    initialize_borders(starting_borders, N_BORDERS);
    for (int i = 1; i < sim_counter; i++)
    {

        if (strstr("vectorize_1", sim_list[i]->name))
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

                check += check_square_distance(oracle_position, current_position, 2 * number_of_people, 0);
                check += check_square_distance(oracle_speed, current_speed, number_of_people, 0);
                check += check_square_distance(oracle_desired_direction, current_desired_direction, 2 * number_of_people, 0);
                check += check_square_distance(oracle_acceleration_term, current_acceleration_term, 2 * number_of_people, 0);
                check += check_square_distance(oracle_people_repulsion_term, current_people_repulsion_term, 2 * number_of_people * number_of_people, 0);

                //check += check_square_distance(oracle_border_repulsion_term, current_border_repulsion_term, 2 * 2 * number_of_people, 0);

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
        }
        else
        {

            copy_init(starting_position, starting_desired_direction, starting_final_destination, starting_borders, starting_desired_speed, starting_desired_max_speed,
                      &oracle_position, &oracle_desired_direction, &oracle_final_destination, &oracle_borders, &oracle_desired_speed, &oracle_desired_max_speed, number_of_people);
            allocate_arrays(&oracle_speed, &oracle_actual_velocity, &oracle_acceleration_term, &oracle_people_repulsion_term,
                            &oracle_border_repulsion_term, &oracle_social_force, number_of_people);

            copy_init_new(starting_position, starting_desired_direction, starting_final_destination, starting_borders, starting_desired_speed, starting_desired_max_speed,
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

                check += check_square_distance(oracle_position, current_position, number_of_people, 1);
                //check += check_square_distance(oracle_speed, current_speed, number_of_people, 0);
                //check += check_square_distance(oracle_desired_direction, current_desired_direction, number_of_people, 1);
                //check += check_square_distance(oracle_acceleration_term, current_acceleration_term, number_of_people, 1);
                //check += check_square_distance(oracle_people_repulsion_term, current_people_repulsion_term, number_of_people, 2);

                //check += check_square_distance(oracle_border_repulsion_term, current_border_repulsion_term, number_of_people, 3);

                if (check)
                {
                    printf("ERROR: implementation %s differs in iteration %d from the base\n", sim_list[i]->name, j);
                    error_check = 1;
                }
                else
                {
                    printf("%s iteration %d CORRECT!\n", sim_list[i]->name, j);
                }

                copy_state_new(oracle_position, oracle_desired_direction, oracle_final_destination, oracle_borders, oracle_desired_speed, oracle_desired_max_speed,
                               &current_position, &current_desired_direction, &current_final_destination, &current_borders, &current_desired_speed, &current_desired_max_speed, number_of_people, sim_list[i]);
            }
        }

        //check_square_distance_vec(oracle_position, current_position, error_vector, 2 * number_of_people);

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
*   Function that runs all the testcases contained inside the list for each of the 
*   implementations added. 
*   Returns 1 if some test fails, 0 otherwise.
*/
int run_testcases()
{
    int error_check = 0;
    for (int k = 0; k < N_TESTS; k++)
    {
        if (testcases[k] == NULL)
            continue;
        char *id;
        switch (k)
        {
        case TDIR:
            id = "direction";
            break;
        case TACC:
            id = "acceleration";
            break;
        case TFRC:
            id = "social force";
            break;
        case TPOS:
            id = "position";
            break;
        default:
            break;
        }
        PRINT_RED;
        printf("Running %s tests\n", id);
        PRINT_DEF;
        float *result;
        size_t result_size;
        float *control;
        int np;
        char *testname;
        testcase_t *cur;
        for (int i = 0; i < counter[2 * k + 1]; i++)
        {
            PRINT_GREEN;
            printf("Implementation %d:\n", i);
            PRINT_DEF;
            for (int j = 0; j < counter[2 * k]; j++)
            {
                cur = testcases[k][j];
                np = cur->n;
                testname = cur->name;
                switch (k)
                {
                case TDIR:
                    result_size = np * 2;
                    void (*f0)(float *, float *, float *, int) = direction_ptr_v[i];
                    result = (float *)calloc(result_size, sizeof(float));
                    f0(cur->pos, cur->des, result, np);
                    control = cur->dir;
                    break;
                case TACC:
                    result_size = np * 2;
                    void (*f1)(float *, float *, float *, float *, int) = acceleration_ptr_v[i];
                    result = (float *)calloc(result_size, sizeof(float));
                    f1(cur->dir, result, cur->vel, cur->spe, np);
                    control = cur->acc;
                    break;
                case TFRC:
                    result_size = np * 2;
                    void (*f4)(float *, float *, float *, float *, int, int) = social_ptr_v[i];
                    result = (float *)calloc(result_size, sizeof(float));
                    f4(cur->acc, cur->pre, cur->bre, result, np, cur->n_borders);
                    control = cur->frc;
                    break;
                case TPOS:
                    result_size = np * 2;
                    void (*f5)(float *, float *, float *, float *, float *, float *, int) = pos_ptr_v[i];
                    f5(cur->pos, cur->dir, cur->spe, cur->frc, cur->vel, cur->dspe, np);
                    control = cur->des;
                    result = cur->pos; //TODO: right now it only check the final position, not the direction
                    break;
                default:
                    break;
                }
                if (check_square_distance(control, result, np, 0))
                {
                    printf("\tTest %s '%s': ERROR!\n", id, testname);
                    error_check = 1;
                }
                else
                {
                    printf("\tTest %s '%s': CORRECT!\n", id, testname);
                }
            }
        }
    }
    return error_check;
}

/*
*   Function used to check the square distance between the expected result of a testcase
*   and the result of the function testes.
*   This function returns 1 if the distance is greater than the threshold, 0 otherwise.
*/
int check_square_distance(float *expected, float *res, int n, int case_n)
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
            printf("x: %lf %lf %lf %d \n", expected[IndexX(i)], res[IndexX(i)], acc, IndexX(i));

            if (isnan(acc) || acc > EPS_NEW)
            {
                wrong = 1;
                printf("x: %lf %lf %lf %d \n", expected[IndexX(i)], res[IndexX(i)], acc, IndexX(i));
            }
        }

        for (int i = 0; i < n; i++)
        {
            acc = fabs(expected[IndexY(i, n)] - res[IndexY(i, n)]);
            printf("y: %lf %lf %lf %d \n", expected[IndexY(i, n)], res[IndexY(i, n)], acc, IndexY(i, n));

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

void add_direction_testcase(char *name, float *pos, float *des, float *output, int n)
{
    counter[2 * TDIR]++;
    testcases[TDIR] = realloc(testcases[TDIR], sizeof(testcase_t *) * counter[2 * TDIR]);
    testcases[TDIR][counter[2 * TDIR] - 1] = malloc(sizeof(testcase_t));
    testcases[TDIR][counter[2 * TDIR] - 1]->name = name;
    testcases[TDIR][counter[2 * TDIR] - 1]->pos = pos;
    testcases[TDIR][counter[2 * TDIR] - 1]->des = des;
    testcases[TDIR][counter[2 * TDIR] - 1]->dir = output;
    testcases[TDIR][counter[2 * TDIR] - 1]->n = n;
}

void add_acceleration_testcase(char *name, float *dir, float *vel, float *spe, float *output, int n)
{
    counter[2 * TACC]++;
    testcases[TACC] = realloc(testcases[TACC], sizeof(testcase_t *) * counter[2 * TACC]);
    testcases[TACC][counter[2 * TACC] - 1] = malloc(sizeof(testcase_t));
    testcases[TACC][counter[2 * TACC] - 1]->name = name;
    testcases[TACC][counter[2 * TACC] - 1]->dir = dir;
    testcases[TACC][counter[2 * TACC] - 1]->vel = vel;
    testcases[TACC][counter[2 * TACC] - 1]->spe = spe;
    testcases[TACC][counter[2 * TACC] - 1]->acc = output;
    testcases[TACC][counter[2 * TACC] - 1]->n = n;
}

void add_compute_social_force_testcase(char *name, float *acc, float *people_rep, float *border_rep, float *output, int n, int n_borders)
{
    counter[2 * TFRC]++;
    testcases[TFRC] = realloc(testcases[TFRC], sizeof(testcase_t *) * counter[2 * TFRC]);
    testcases[TFRC][counter[2 * TFRC] - 1] = malloc(sizeof(testcase_t));
    testcases[TFRC][counter[2 * TFRC] - 1]->name = name;
    testcases[TFRC][counter[2 * TFRC] - 1]->acc = acc;
    testcases[TFRC][counter[2 * TFRC] - 1]->pre = people_rep;
    testcases[TFRC][counter[2 * TFRC] - 1]->bre = border_rep;
    testcases[TFRC][counter[2 * TFRC] - 1]->frc = output;
    testcases[TFRC][counter[2 * TFRC] - 1]->n = n;
    testcases[TFRC][counter[2 * TFRC] - 1]->n_borders = n_borders;
}

void add_position_testcase(char *name, float *pos, float *dir, float *speed, float *force, float *actual_vel, float *desired_speed, float *output_pos, int n)
{
    counter[2 * TPOS]++;
    testcases[TPOS] = realloc(testcases[TPOS], sizeof(testcase_t *) * counter[2 * TPOS]);
    testcases[TPOS][counter[2 * TPOS] - 1] = malloc(sizeof(testcase_t));
    testcases[TPOS][counter[2 * TPOS] - 1]->name = name;
    testcases[TPOS][counter[2 * TPOS] - 1]->pos = pos;
    testcases[TPOS][counter[2 * TPOS] - 1]->spe = speed;
    testcases[TPOS][counter[2 * TPOS] - 1]->frc = force;
    testcases[TPOS][counter[2 * TPOS] - 1]->vel = actual_vel;
    testcases[TPOS][counter[2 * TPOS] - 1]->dir = dir;
    testcases[TPOS][counter[2 * TPOS] - 1]->dspe = desired_speed;
    testcases[TPOS][counter[2 * TPOS] - 1]->des = output_pos;
    testcases[TPOS][counter[2 * TPOS] - 1]->n = n;
}

void add_direction_implementation(void (*f)(float *, float *, float *, int))
{
    counter[2 * TDIR + 1]++;
    direction_ptr_v = realloc(direction_ptr_v, counter[2 * TDIR + 1] * sizeof(void (*)(float *, float *, float *, int)));
    direction_ptr_v[counter[2 * TDIR + 1] - 1] = f;
}

void add_acceleration_implementation(void (*f)(float *, float *, float *, float *, int))
{
    counter[2 * TACC + 1]++;
    acceleration_ptr_v = realloc(acceleration_ptr_v, counter[2 * TACC + 1] * sizeof(void (*)(float *, float *, float *, float *, int)));
    acceleration_ptr_v[counter[2 * TACC + 1] - 1] = f;
}

void add_social_implementation(void (*f)(float *, float *, float *, float *, int, int))
{
    counter[2 * TFRC + 1]++;
    social_ptr_v = realloc(social_ptr_v, counter[2 * TFRC + 1] * sizeof(void (*)(float *, float *, float *, float *, int, int)));
    social_ptr_v[counter[2 * TFRC + 1] - 1] = f;
}

void add_pos_implementation(void (*f)(float *, float *, float *, float *, float *, float *, int))
{
    counter[2 * TPOS + 1]++;
    pos_ptr_v = realloc(pos_ptr_v, counter[2 * TPOS + 1] * sizeof(void (*)(float *, float *, float *, float *, float *, float *, int)));
    pos_ptr_v[counter[2 * TPOS + 1] - 1] = f;
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

void copy_init_new(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
                   float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n)
{

    *pos = (float *)aligned_malloc(n * 2 * sizeof(float), 32);
    *dir = (float *)aligned_malloc(n * 2 * sizeof(float), 32);
    *fdes = (float *)aligned_malloc(n * 2 * sizeof(float), 32);
    *bor = (float *)aligned_malloc(N_BORDERS * sizeof(float), 32);
    *spe = (float *)aligned_malloc(n * sizeof(float), 32);
    *mspe = (float *)aligned_malloc(n * sizeof(float), 32);

    for (int i = 0; i < n; i++)
    {
        float *posk = *pos;
        float *dirk = *dir;
        float *fdesk = *fdes;

        posk[IndexX(i)] = s_pos[IndexX(i)];
        posk[IndexY(i, n)] = s_pos[IndexY(i, n)];

        dirk[IndexX(i)] = s_dir[IndexX(i)];
        dirk[IndexY(i, n)] = s_dir[IndexY(i, n)];

        fdesk[IndexX(i)] = s_fdes[IndexX(i)];
        fdesk[IndexY(i, n)] = s_fdes[IndexY(i, n)];
    }

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

void copy_state_new(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
                    float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n, sim_t *sim)
{

    for (int i = 0; i < n; i++)
    {
        float *posk = *pos;
        float *dirk = *dir;
        float *fdesk = *fdes;

        posk[IndexX(i)] = s_pos[IndexX(i)];
        posk[IndexY(i, n)] = s_pos[IndexY(i, n)];

        if (!strstr("vectorize_2_5_1", sim->name))
        {
            dirk[IndexX(i)] = s_dir[IndexX(i)];
            dirk[IndexY(i, n)] = s_dir[IndexY(i, n)];

            fdesk[IndexX(i)] = s_fdes[IndexX(i)];
            fdesk[IndexY(i, n)] = s_fdes[IndexY(i, n)];
        }
    }
    if (!strstr("vectorize_2_5_1", sim->name))
    {
        memcpy(*bor, s_bor, N_BORDERS * sizeof(float));
        memcpy(*spe, s_spe, n * sizeof(float));
        memcpy(*mspe, s_mspe, n * sizeof(float));
    }
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

/* finite-differences functions */

float compute_people_repulsion_fd(float *position, float *parameters, float *desired_direction, float *actual_speed, int i, int j)
{
    float rx_ab = parameters[0];
    float ry_ab = parameters[1];
    float ex_a = desired_direction[i * 2];
    float ey_a = desired_direction[i * 2 + 1];
    float ex_b = desired_direction[j * 2];
    float ey_b = desired_direction[j * 2 + 1];
    float vb = actual_speed[j];
    float delta_b = vb * TIMESTEP;

    float r_ab_norm = sqrt(rx_ab * rx_ab + ry_ab * ry_ab);

    float rx_ab_mex = rx_ab - delta_b * ex_b;
    float ry_ab_mey = ry_ab - delta_b * ey_b;

    float r_ab_me_norm = sqrt(rx_ab_mex * rx_ab_mex + ry_ab_mey * ry_ab_mey);

    float b = sqrt((r_ab_norm + r_ab_me_norm) * (r_ab_norm + r_ab_me_norm) - (delta_b * delta_b)) / 2;

    return V_ALPHA_BETA * exp(-b / SIGMA);
}

void add_estimated_gradient_people_repulsion(float *gradient, float *parameters, int n_parameters, float *position, float *desired_direction, float *actual_speed, int i, int j)
{

    float dp = 1.0e-7;
    float f_P, f_M;

    for (int p = 0; p < n_parameters; p++)
    {
        float tmpVal = parameters[p];
        parameters[p] = tmpVal + dp;
        f_P = compute_people_repulsion_fd(position, parameters, desired_direction, actual_speed, i, j);

        parameters[p] = tmpVal - dp;
        f_M = compute_people_repulsion_fd(position, parameters, desired_direction, actual_speed, i, j);

        gradient[p] = (f_P - f_M) / (2 * dp);
    }
    float ex_a = desired_direction[i * 2];
    float ey_a = desired_direction[i * 2 + 1];
    float check = ex_a * (gradient[0]) + ey_a * (gradient[1]);
    float threshold = sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]) * cos(PSI);
    float w = check >= threshold ? 1 : INFLUENCE;
    gradient[0] = w * -gradient[0];
    gradient[1] = w * -gradient[1];
}

void test_people_repulsion_with_FD(float *people_repulsion_term, int n, float *position, float *desired_direction, float *actual_speed)
{
    int num_errors = 0;
    float tol = 1e-4;
    float eps = 1e-10;

    float fd_gradient[2];
    float parameters[2];
    float *analytic_gradient = people_repulsion_term;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            float rx_ab = position[i * 2] - position[j * 2];
            float ry_ab = position[i * 2 + 1] - position[j * 2 + 1];
            parameters[0] = rx_ab;
            parameters[1] = ry_ab;
            add_estimated_gradient_people_repulsion(fd_gradient, parameters, 2, position, desired_direction, actual_speed, i, j);
            for (int p = 0; p < 2; p++)
            {
                float absErr = fabs(fd_gradient[p] - analytic_gradient[i * (2 * n) + 2 * j + p]);
                float relError = 2 * absErr / (eps + analytic_gradient[i * (2 * n) + 2 * j + p] + fd_gradient[p]);

                if (relError > tol && absErr > EPS)
                {
                    printf("Mismatch in people_repulsion element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
                    printf("Mismatch in people_repulsion element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
                    num_errors++;
                    if (num_errors > 51)
                    {
                        exit(1);
                    }
                }
            }
        }
    }
}
// default values

float compute_border_repulsion_fd(float *position, float *parameters, float *borders, int i, int j)
{
    float rx_a = parameters[0];
    float ry_a = parameters[1];

    float rx_aB = 0.0;
    float ry_aB = ry_a - borders[j];

    float r_aB_norm = fabs(ry_aB);

    return U_ALPHA_B * exp(-r_aB_norm / R);
}

void add_estimated_gradient_border_repulsion(float *gradient, float *parameters, int n_parameters, float *position, float *borders, int i, int j)
{

    float dp = 1.0e-7;
    float f_P, f_M;

    for (int p = 0; p < n_parameters; p++)
    {
        float tmpVal = parameters[p];
        parameters[p] = tmpVal + dp;
        f_P = compute_border_repulsion_fd(position, parameters, borders, i, j);

        parameters[p] = tmpVal - dp;
        f_M = compute_border_repulsion_fd(position, parameters, borders, i, j);

        gradient[p] = (f_P - f_M) / (2 * dp);
    }
    gradient[0] *= -1;
    gradient[1] *= -1;
}

void test_border_repulsion_with_FD(float *border_repulsion_term, float *position, float *borders, int n_borders, int n)
{
    int num_errors = 0;
    float tol = 1e-4;
    float eps = 1e-10;

    float fd_gradient[2];
    float parameters[2];
    float *analytic_gradient = border_repulsion_term;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n_borders; j++)
        {
            float rx_a = position[i * 2];
            float ry_a = position[i * 2 + 1];
            float rx_aB = 0.0;
            float ry_aB = ry_a - borders[j];

            parameters[0] = rx_aB;
            parameters[1] = ry_aB;
            add_estimated_gradient_border_repulsion(fd_gradient, parameters, 2, position, borders, i, j);
            for (int p = 0; p < 2; p++)
            {
                float absErr = fabs(fd_gradient[p] - analytic_gradient[i * (2 * n_borders) + 2 * j + p]);
                float relError = 2 * absErr / (eps + analytic_gradient[i * (2 * n_borders) + 2 * j + p] + fd_gradient[p]);

                if (relError > tol && absErr > EPS)
                {
                    printf("Mismatch in border_repulsion: element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n_borders) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
                    printf("Mismatch in border_repulsion: element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n_borders) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
                    num_errors++;
                    if (num_errors > 51)
                    {
                        exit(1);
                    }
                }
            }
        }
    }
}
