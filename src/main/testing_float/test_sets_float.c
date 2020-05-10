

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
#include "../social_force_model_basic_simplified.h"
#include "utility.h"

// basic direction test case for 4 persons
float direction_position0[] = {0, 0, 10, 10, 0, 0, 10, 10};
float direction_fdest0[] = {10, 10, 0, 0, 10, 10, 0, 0};
float direction_expected0[] = {0.7071, 0.7071, -0.7071, -0.7071, 0.7071, 0.7071, -0.7071, -0.7071};
int direction_n0 = 4;

float acceleration_desired_speed[] = {1.34, 1.34, 1.34, 1.34};
// acceleration from start straight and 45°
float acceleration_direction0[] = {1, 0, -1, 0, 1, 1, -1, 1};
float acceleration_expected0[] = {2.68, 0, -2.68, 0, 2.68, 2.68, -2.68, 2.68};
float acceleration_vel0[] = {0, 0, 0, 0, 0, 0, 0, 0};
int acceleration_n0 = 4;
// acceleration with velocity straight and 45°
float acceleration_direction1[] = {1, 0, -1, 0, 1, 1, -1, 1};
float acceleration_expected1[] = {1.68, 0, -1.68, 0, 1.68, 2.68, -1.68, 2.68};
float acceleration_vel1[] = {0.5, 0, -0.5, 0, 0.5, 0, -0.5, 0};
int acceleration_n1 = 4;
// acceleration to zero straight and 45°
float acceleration_direction2[] = {1, 0, -1, 0, 1, 1, -1, 1};
float acceleration_expected2[] = {0, 0, 0, 0, 0, 0, 0, 0};
float acceleration_vel2[] = {1.34, 0, -1.34, 0, 1.34, 1.34, -1.34, 1.34};
int acceleration_n2 = 4;
// deacceleration straight and 45°
float acceleration_direction3[] = {1, 0, -1, 0, 1, 1, -1, 1};
float acceleration_expected3[] = {-0.8, 0, 0.8, 0, -0.8, 2.68, 0.8, 2.68};
float acceleration_vel3[] = {1.74, 0, -1.74, 0, 1.74, 0, -1.74, 0};
int acceleration_n3 = 4;

float repulsion_position0[] = {5, 3, 8, 2, 9, 3};
float repulsion_speed0[] = {0.7, 0.5, 1.3};
float repulsion_direction0[] = {0, 1, 1, 0, 4, 3};
int repulsion_n0 = 3;
float repulsion_expected_result0[] = {0, 0, -0.000020266568837, 0.000005901638773, -0.000000000000016, -0.000000000000004,
                                      0.000033293428727, -0.000018112519707, 0, 0, -0.000000000027150, -0.000000000023879,
                                      0.000004591945972, -0.000000780381604, 0.034905007547125, 0.084268142615003, 0, 0};

float border_pos0[] = {0, 0, 10, 10};
float border_borders0[] = {0, 10};
float border_expected0[] = {1, 0, 1, 1, 0, 0, 1, 0};
int border_n0 = 2;
int border_nb0 = 2;

// basic social force test for 4 persons and 2 walls
float social_acc0[] = {1, 2, -3, 5, 2, 2, -4, -2};
float social_prep0[] = {0, 0, 1, 1, 2, 3, -1, 1,
                        -1, -1, 0, 0, 2, 2, -2, 1,
                        -2, -3, -2, -2, 0, 0, 3, -1,
                        1, -1, 2, -1, -3, 1, 0, 0};
float social_brep0[] = {0, 1, 1, 1,
                        -1, 1, 0, 1,
                        0, 0, 1, 2,
                        2, 1, -1, 2};
float social_expected0[] = {4, 9, -5, 9, 2, -2, -3, 0};
int social_n0 = 4;
int social_nb0 = 2;

float position_desired_speed[] = {1.34, 1.34, 1.34, 1.34};
float position_dir[] = {0, 0, 0, 0, 0, 0, 0, 0}; // not needed for computation, only written to
float position_speed[] = {0, 0, 0, 0};           // not needed for computation, only written to
// Test0, 4 unit directions, 0 starting speed, unit social force
float position_pos0[] = {0, 0, 0, 0, 0, 0, 0, 0};     // starting at origin
float position_force0[] = {1, 0, 0, 1, -1, 0, 0, -1}; // unit social force
float position_vel0[] = {0, 0, 0, 0, 0, 0, 0, 0};
float position_expected0[] = {0.04, 0, 0, 0.04, -0.04, 0, 0, -0.04};
int position_n0 = 4;
// Test1, 4 unit directions, 1.742 starting speed, unit social force
float position_pos1[] = {0, 0, 0, 0, 0, 0, 0, 0};     // starting at origin
float position_force1[] = {1, 0, 0, 1, -1, 0, 0, -1}; // unit social force
float position_vel1[] = {1.742, 0, 0, 1.742, -1.742, 0, 0, -1.742};
float position_expected1[] = {0.3484, 0, 0, 0.3484, -0.3484, 0, 0, -0.3484};
int position_n1 = 4;

/* global variables for test cases */

testcase_t **testcases[N_TESTS];

int counter[2 * N_TESTS];

void (**direction_ptr_v)(float *, float *, float *, int);
void (**acceleration_ptr_v)(float *, float *, float *, float *, int);
void (**social_ptr_v)(float *, float *, float *, float *, int, int);
void (**pos_ptr_v)(float *, float *, float *, float *, float *, float *, int);

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
    add_direction_implementation(update_desired_direction_simplified);

    //update acceleration implementations
    add_acceleration_implementation(update_acceleration_term_simplified);

    //compute social force implementations
    add_social_implementation(compute_social_force_simplified);

    //update position implementations
    add_pos_implementation(update_position_simplified);
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
                if (check_absolute_distance(control, result, np, 0))
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
