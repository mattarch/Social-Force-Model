

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// fix aligenment
#include "../aligned_free.h"
#include "../aligned_malloc.h"

#include "../parse_args.h"
#include "../social_force.h"
#include "../social_force_model_basic_simplified_double.h"
#include "../utility.h"

#include "compare_simulations_double.h"
#include "test_sets_double.h"

// basic direction test case for 4 persons
double direction_position0_double[] = {0, 0, 10, 10, 0, 0, 10, 10};
double direction_fdest0_double[] = {10, 10, 0, 0, 10, 10, 0, 0};
double direction_expected0_double[] = {0.7071, 0.7071, -0.7071, -0.7071, 0.7071, 0.7071, -0.7071, -0.7071};
int direction_n0_double = 4;

double acceleration_desired_speed_double[] = {1.34, 1.34, 1.34, 1.34};
// acceleration from start straight and 45°
double acceleration_direction0_double[] = {1, 0, -1, 0, 1, 1, -1, 1};
double acceleration_expected0_double[] = {2.68, 0, -2.68, 0, 2.68, 2.68, -2.68, 2.68};
double acceleration_vel0_double[] = {0, 0, 0, 0, 0, 0, 0, 0};
int acceleration_n0_double = 4;
// acceleration with velocity straight and 45°
double acceleration_direction1_double[] = {1, 0, -1, 0, 1, 1, -1, 1};
double acceleration_expected1_double[] = {1.68, 0, -1.68, 0, 1.68, 2.68, -1.68, 2.68};
double acceleration_vel1_double[] = {0.5, 0, -0.5, 0, 0.5, 0, -0.5, 0};
int acceleration_n1_double = 4;
// acceleration to zero straight and 45°
double acceleration_direction2_double[] = {1, 0, -1, 0, 1, 1, -1, 1};
double acceleration_expected2_doubel[] = {0, 0, 0, 0, 0, 0, 0, 0};
double acceleration_vel2_double[] = {1.34, 0, -1.34, 0, 1.34, 1.34, -1.34, 1.34};
int acceleration_n2_double = 4;
// deacceleration straight and 45°
double acceleration_direction3_double[] = {1, 0, -1, 0, 1, 1, -1, 1};
double acceleration_expected3_double[] = {-0.8, 0, 0.8, 0, -0.8, 2.68, 0.8, 2.68};
double acceleration_vel3_double[] = {1.74, 0, -1.74, 0, 1.74, 0, -1.74, 0};
int acceleration_n3_double = 4;

double repulsion_position0_double[] = {5, 3, 8, 2, 9, 3};
double repulsion_speed0_double[] = {0.7, 0.5, 1.3};
double repulsion_direction0_double[] = {0, 1, 1, 0, 4, 3};
int repulsion_n0_double = 3;
double repulsion_expected_result0_double[] = {0, 0, -0.000020266568837, 0.000005901638773, -0.000000000000016, -0.000000000000004,
                                              0.000033293428727, -0.000018112519707, 0, 0, -0.000000000027150, -0.000000000023879,
                                              0.000004591945972, -0.000000780381604, 0.034905007547125, 0.084268142615003, 0, 0};

double border_pos0_double[] = {0, 0, 10, 10};
double border_borders0_double[] = {0, 10};
double border_expected0_double[] = {1, 0, 1, 1, 0, 0, 1, 0};
int border_n0_double = 2;
int border_nb0_double = 2;

// basic social force test for 4 persons and 2 walls
double social_acc0_double[] = {1, 2, -3, 5, 2, 2, -4, -2};
double social_prep0_double[] = {0, 0, 1, 1, 2, 3, -1, 1,
                                -1, -1, 0, 0, 2, 2, -2, 1,
                                -2, -3, -2, -2, 0, 0, 3, -1,
                                1, -1, 2, -1, -3, 1, 0, 0};
double social_brep0_double[] = {0, 1, 1, 1,
                                -1, 1, 0, 1,
                                0, 0, 1, 2,
                                2, 1, -1, 2};
double social_expected0_double[] = {4, 9, -5, 9, 2, -2, -3, 0};
int social_n0_double = 4;
int social_nb0_double = 2;

double position_desired_speed_double[] = {1.34, 1.34, 1.34, 1.34};
double position_dir_double[] = {0, 0, 0, 0, 0, 0, 0, 0}; // not needed for computation, only written to
double position_speed_double[] = {0, 0, 0, 0};           // not needed for computation, only written to
// Test0, 4 unit directions, 0 starting speed, unit social force
double position_pos0_double[] = {0, 0, 0, 0, 0, 0, 0, 0};     // starting at origin
double position_force0_double[] = {1, 0, 0, 1, -1, 0, 0, -1}; // unit social force
double position_vel0_double[] = {0, 0, 0, 0, 0, 0, 0, 0};
double position_expected0_double[] = {0.04, 0, 0, 0.04, -0.04, 0, 0, -0.04};
int position_n0_double = 4;
// Test1, 4 unit directions, 1.742 starting speed, unit social force
double position_pos1_double[] = {0, 0, 0, 0, 0, 0, 0, 0};     // starting at origin
double position_force1_double[] = {1, 0, 0, 1, -1, 0, 0, -1}; // unit social force
double position_vel1_double[] = {1.742, 0, 0, 1.742, -1.742, 0, 0, -1.742};
double position_expected1_double[] = {0.3484, 0, 0, 0.3484, -0.3484, 0, 0, -0.3484};
int position_n1_double = 4;

/* global variables for test cases */

testcase_t_double **testcases[N_TESTS];

int counter[2 * N_TESTS];

void (**direction_ptr_v)(double *, double *, double *, int);
void (**acceleration_ptr_v)(double *, double *, double *, double *, int);
void (**social_ptr_v)(double *, double *, double *, double *, int, int);
void (**pos_ptr_v)(double *, double *, double *, double *, double *, double *, int);

/*
*   Function that adds all the testcases to try to the list
*/
void add_testcases_double()
{
    //direction testcases
    add_direction_testcase_double("basic direction test", direction_position0_double, direction_fdest0_double, direction_expected0_double, direction_n0_double);

    //acceleration testcases
    add_acceleration_testcase_double("acceleration_test_from_start_straight_and_with_angle", acceleration_direction0_double, acceleration_vel0_double, acceleration_desired_speed_double, acceleration_expected0_double, acceleration_n0_double);
    add_acceleration_testcase_double("acceleration_test_with_velocity_straight_and_with_angle", acceleration_direction1_double, acceleration_vel1_double, acceleration_desired_speed_double, acceleration_expected1_double, acceleration_n1_double);
    add_acceleration_testcase_double("acceleration_test_to_0_straight_and_with_angle", acceleration_direction2_double, acceleration_vel2_double, acceleration_desired_speed_double, acceleration_expected2_doubel, acceleration_n2_double);
    add_acceleration_testcase_double("deacceleration_test_straight_and_with_angle", acceleration_direction3_double, acceleration_vel3_double, acceleration_desired_speed_double, acceleration_expected3_double, acceleration_n3_double);

    //social force
    add_compute_social_force_testcase_double("basic test", social_acc0_double, social_prep0_double, social_brep0_double, social_expected0_double, social_n0_double, social_nb0_double);

    //position
    add_position_testcase_double("position_test_from_origin_0_speed_unit_social_force", position_pos0_double, position_dir_double, position_speed_double, position_force0_double, position_vel0_double, position_desired_speed_double, position_expected0_double, position_n0_double);
}

/*
*   function that adds different implementations to test for correctness
*/
void add_function_implementations_double()
{
    //update direction implementations
    add_direction_implementation_double(update_desired_direction_simplified_double);

    //update acceleration implementations
    add_acceleration_implementation_double(update_acceleration_term_simplified_double);

    //compute social force implementations
    add_social_implementation_double(compute_social_force_simplified_double);

    //update position implementations
    add_pos_implementation_double(update_position_simplified_double);
}

int run_tests_double(sim_t **sim_list, int sim_counter)
{
    //init auxilliary parameters
    for (int i = 0; i < N_TESTS * 2; i++)
    {
        counter[i] = 0;
    }

    //add tests
    add_testcases_double();
    add_function_implementations_double();

    //run tests
    int error = run_testcases_double();
    int errorsim = compare_simulations_double(sim_list, sim_counter);
    return error + errorsim;
}

/*
*   Function that runs all the testcases contained inside the list for each of the 
*   implementations added. 
*   Returns 1 if some test fails, 0 otherwise.
*/
int run_testcases_double()
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
        double *result;
        size_t result_size;
        double *control;
        int np;
        char *testname;
        testcase_t_double *cur;
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
                    void (*f0)(double *, double *, double *, int) = direction_ptr_v[i];
                    result = (double *)calloc(result_size, sizeof(double));
                    f0(cur->pos, cur->des, result, np);
                    control = cur->dir;
                    break;
                case TACC:
                    result_size = np * 2;
                    void (*f1)(double *, double *, double *, double *, int) = acceleration_ptr_v[i];
                    result = (double *)calloc(result_size, sizeof(double));
                    f1(cur->dir, result, cur->vel, cur->spe, np);
                    control = cur->acc;
                    break;
                case TFRC:
                    result_size = np * 2;
                    void (*f4)(double *, double *, double *, double *, int, int) = social_ptr_v[i];
                    result = (double *)calloc(result_size, sizeof(double));
                    f4(cur->acc, cur->pre, cur->bre, result, np, cur->n_borders);
                    control = cur->frc;
                    break;
                case TPOS:
                    result_size = np * 2;
                    void (*f5)(double *, double *, double *, double *, double *, double *, int) = pos_ptr_v[i];
                    f5(cur->pos, cur->dir, cur->spe, cur->frc, cur->vel, cur->dspe, np);
                    control = cur->des;
                    result = cur->pos; //TODO: right now it only check the final position, not the direction
                    break;
                default:
                    break;
                }
                if (check_absolute_distance_double(control, result, np, 0))
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

void add_direction_testcase_double(char *name, double *pos, double *des, double *output, int n)
{
    counter[2 * TDIR]++;
    testcases[TDIR] = realloc(testcases[TDIR], sizeof(testcase_t_double *) * counter[2 * TDIR]);
    testcases[TDIR][counter[2 * TDIR] - 1] = malloc(sizeof(testcase_t_double));
    testcases[TDIR][counter[2 * TDIR] - 1]->name = name;
    testcases[TDIR][counter[2 * TDIR] - 1]->pos = pos;
    testcases[TDIR][counter[2 * TDIR] - 1]->des = des;
    testcases[TDIR][counter[2 * TDIR] - 1]->dir = output;
    testcases[TDIR][counter[2 * TDIR] - 1]->n = n;
}

void add_acceleration_testcase_double(char *name, double *dir, double *vel, double *spe, double *output, int n)
{
    counter[2 * TACC]++;
    testcases[TACC] = realloc(testcases[TACC], sizeof(testcase_t_double *) * counter[2 * TACC]);
    testcases[TACC][counter[2 * TACC] - 1] = malloc(sizeof(testcase_t_double));
    testcases[TACC][counter[2 * TACC] - 1]->name = name;
    testcases[TACC][counter[2 * TACC] - 1]->dir = dir;
    testcases[TACC][counter[2 * TACC] - 1]->vel = vel;
    testcases[TACC][counter[2 * TACC] - 1]->spe = spe;
    testcases[TACC][counter[2 * TACC] - 1]->acc = output;
    testcases[TACC][counter[2 * TACC] - 1]->n = n;
}

void add_compute_social_force_testcase_double(char *name, double *acc, double *people_rep, double *border_rep, double *output, int n, int n_borders)
{
    counter[2 * TFRC]++;
    testcases[TFRC] = realloc(testcases[TFRC], sizeof(testcase_t_double *) * counter[2 * TFRC]);
    testcases[TFRC][counter[2 * TFRC] - 1] = malloc(sizeof(testcase_t_double));
    testcases[TFRC][counter[2 * TFRC] - 1]->name = name;
    testcases[TFRC][counter[2 * TFRC] - 1]->acc = acc;
    testcases[TFRC][counter[2 * TFRC] - 1]->pre = people_rep;
    testcases[TFRC][counter[2 * TFRC] - 1]->bre = border_rep;
    testcases[TFRC][counter[2 * TFRC] - 1]->frc = output;
    testcases[TFRC][counter[2 * TFRC] - 1]->n = n;
    testcases[TFRC][counter[2 * TFRC] - 1]->n_borders = n_borders;
}

void add_position_testcase_double(char *name, double *pos, double *dir, double *speed, double *force, double *actual_vel, double *desired_speed, double *output_pos, int n)
{
    counter[2 * TPOS]++;
    testcases[TPOS] = realloc(testcases[TPOS], sizeof(testcase_t_double *) * counter[2 * TPOS]);
    testcases[TPOS][counter[2 * TPOS] - 1] = malloc(sizeof(testcase_t_double));
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

void add_direction_implementation_double(void (*f)(double *, double *, double *, int))
{
    counter[2 * TDIR + 1]++;
    direction_ptr_v = realloc(direction_ptr_v, counter[2 * TDIR + 1] * sizeof(void (*)(double *, double *, double *, int)));
    direction_ptr_v[counter[2 * TDIR + 1] - 1] = f;
}

void add_acceleration_implementation_double(void (*f)(double *, double *, double *, double *, int))
{
    counter[2 * TACC + 1]++;
    acceleration_ptr_v = realloc(acceleration_ptr_v, counter[2 * TACC + 1] * sizeof(void (*)(double *, double *, double *, double *, int)));
    acceleration_ptr_v[counter[2 * TACC + 1] - 1] = f;
}

void add_social_implementation_double(void (*f)(double *, double *, double *, double *, int, int))
{
    counter[2 * TFRC + 1]++;
    social_ptr_v = realloc(social_ptr_v, counter[2 * TFRC + 1] * sizeof(void (*)(double *, double *, double *, double *, int, int)));
    social_ptr_v[counter[2 * TFRC + 1] - 1] = f;
}

void add_pos_implementation_double(void (*f)(double *, double *, double *, double *, double *, double *, int))
{
    counter[2 * TPOS + 1]++;
    pos_ptr_v = realloc(pos_ptr_v, counter[2 * TPOS + 1] * sizeof(void (*)(double *, double *, double *, double *, double *, double *, int)));
    pos_ptr_v[counter[2 * TPOS + 1] - 1] = f;
}
