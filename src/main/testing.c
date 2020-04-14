/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "testing.h"
#include "test_sets.h"
#include "social_force.h"
#include "social_force_model_basic.h"
#include "utility.h"

testcase_t **testcases[N_TESTS];

int counter[2 * N_TESTS];

void (**direction_ptr_v)(double *, double *, double *, int);
void (**acceleration_ptr_v)(double *, double *, double *, double *, int);
void (**social_ptr_v)(double *, double *, double *, double *, int, int);
void (**pos_ptr_v)(double *, double *, double *, double *, double *, double *, int);

simul_t *sim_list;
int sim_counter;

/*
*   Function that adds all the testcases to try to the list
*/
void add_tests()
{
    //direction testcases
    add_direction_testcase("basic direction test", direction_position0, direction_fdest0, direction_expected0, direction_n0);

    //acceleration testcases
    add_acceleration_testcase("acceleration_test_from_start_straight_and_with_angle", acceleration_direction0, acceleration_vel0, acceleration_desired_speed , acceleration_expected0, acceleration_n0);
    add_acceleration_testcase("acceleration_test_with_velocity_straight_and_with_angle", acceleration_direction1, acceleration_vel1, acceleration_desired_speed , acceleration_expected1, acceleration_n1);
    add_acceleration_testcase("acceleration_test_to_0_straight_and_with_angle", acceleration_direction2, acceleration_vel2, acceleration_desired_speed , acceleration_expected2, acceleration_n2);
    add_acceleration_testcase("deacceleration_test_straight_and_with_angle", acceleration_direction3, acceleration_vel3, acceleration_desired_speed , acceleration_expected3, acceleration_n3);

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
    add_simulation_implementation(simulation_basic);
    add_simulation_implementation(simulation_basic);
    add_simulation_implementation(simulation_basic);
    add_simulation_implementation(simulation_basic);

    //update direction implementations
    add_direction_implementation(update_desired_direction);

    //update acceleration implementations
    add_acceleration_implementation(update_acceleration_term);

    //compute social force implementations
    add_social_implementation(compute_social_force);

    //update position implementations
    add_pos_implementation(update_position);
}

int run_tests()
{
    //init auxilliary parameters
    for (int i = 0; i < N_TESTS * 2; i++)
    {
        counter[i] = 0;
    }

    //add tests
    add_tests();
    add_function_implementations();

    //run tests
    int error = run();
    int errorsim = compare_simulations();
    return error + errorsim;
}


/*
* Runs all the simulation and compare the result to the basic version
* Returns 0 if OK, returns 1 if there's an error.
*/
int compare_simulations()
{
    int error_check = 0;
    int number_of_people = 10;
    int n_timesteps = 50;
    double *oracle_position, *oracle_speed, *oracle_desired_direction, *oracle_final_destination,
        *oracle_borders, *oracle_actual_velocity, *oracle_acceleration_term, *oracle_people_repulsion_term,
        *oracle_border_repulsion_term, *oracle_social_force, *oracle_desired_speed;
    double *current_position, *current_speed, *current_desired_direction, *current_final_destination,
        *current_borders, *current_actual_velocity, *current_acceleration_term, *current_people_repulsion_term,
        *current_border_repulsion_term, *current_social_force, *current_desired_speed;
    // allocate memory
    double *starting_position = (double *)calloc(number_of_people * 2, sizeof(double));
    double *starting_desired_direction = (double *)calloc(number_of_people * 2, sizeof(double));
    double *starting_final_destination = (double *)calloc(number_of_people * 2, sizeof(double));
    double *starting_borders = (double *)calloc(N_BORDERS, sizeof(double));
    double *starting_desired_speed = (double *)calloc(number_of_people, sizeof(double));
    // // check if calloc worked correctly
    if (starting_position == NULL || starting_desired_direction == NULL || starting_final_destination == NULL || starting_borders == NULL || starting_desired_speed == NULL)
    {
        printf("Error: calloc failed\n");
        return 1;
    }

    // initialize starting position
    initialize_people(starting_position, starting_desired_direction, starting_final_destination, starting_desired_speed, number_of_people);
    initialize_borders(starting_borders, N_BORDERS);

    copy_init(starting_position, starting_desired_direction, starting_final_destination, starting_borders, starting_desired_speed,
              &oracle_position, &oracle_desired_direction, &oracle_final_destination, &oracle_borders, &oracle_desired_speed, number_of_people);
    allocate_arrays(&oracle_speed, &oracle_actual_velocity, &oracle_acceleration_term, &oracle_people_repulsion_term,
                    &oracle_border_repulsion_term, &oracle_social_force, number_of_people);

    simul_t f = sim_list[0];
    f(number_of_people, n_timesteps, oracle_position, oracle_speed, oracle_desired_direction, oracle_final_destination,
      oracle_borders, oracle_actual_velocity, oracle_acceleration_term,
      oracle_people_repulsion_term, oracle_border_repulsion_term, oracle_social_force, oracle_desired_speed);
    //test all implementations
    for (int i = 1; i < sim_counter; i++)
    {
        int check = 0;

        copy_init(starting_position, starting_desired_direction, starting_final_destination, starting_borders, starting_desired_speed,
                  &current_position, &current_desired_direction, &current_final_destination, &current_borders, &current_desired_speed, number_of_people);
        allocate_arrays(&current_speed, &current_actual_velocity, &current_acceleration_term, &current_people_repulsion_term,
                        &current_border_repulsion_term, &current_social_force, number_of_people);

        f = sim_list[i];
        f(number_of_people, n_timesteps, current_position, current_speed, current_desired_direction, current_final_destination,
          current_borders, current_actual_velocity, current_acceleration_term,
          current_people_repulsion_term, current_border_repulsion_term, current_social_force, current_desired_speed);

        check += check_square_distance(oracle_position, current_position, number_of_people * 2);
        check += check_square_distance(oracle_speed, current_speed, number_of_people);
        check += check_square_distance(oracle_desired_direction, current_desired_direction, number_of_people * 2);

        if (check)
        {
            printf("ERROR: implementation %d differs from the base\n", i);
            error_check = 1;
        }
        else
        {
            printf("CORRECT!\n");
        }

        free_all(11, &current_position, &current_speed, &current_desired_direction, &current_final_destination, &current_borders, &current_actual_velocity,
                 &current_acceleration_term, &current_people_repulsion_term, &current_border_repulsion_term, &current_social_force, &current_desired_speed);
    }

    free_all(11, &oracle_position, &oracle_speed, &oracle_desired_direction, &oracle_final_destination, &oracle_borders, &oracle_actual_velocity,
             &oracle_acceleration_term, &oracle_people_repulsion_term, &oracle_border_repulsion_term, &oracle_social_force, &oracle_desired_speed);

    return error_check;
}

/*
*   Function that runs all the testcases contained inside the list for each of the 
*   implementations added. 
*   Returns 1 if some test fails, 0 otherwise.
*/
int run()
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
                if (check_square_distance(control, result, np))
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
int check_square_distance(double *expected, double *res, int n)
{

    double acc = 0.0;
    for (int i = 0; i < n; i++)
    {
        acc += (expected[i] - res[i]) * (expected[i] - res[i]);
    }

    return isnan(acc) || acc > EPS;
}

void add_direction_testcase(char *name, double *pos, double *des, double *output, int n)
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

void add_acceleration_testcase(char *name, double *dir, double *vel, double *spe, double *output, int n)
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

void add_compute_social_force_testcase(char *name, double *acc, double *people_rep, double *border_rep, double *output, int n, int n_borders)
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

void add_position_testcase(char *name, double *pos, double *dir, double *speed, double *force, double *actual_vel, double *desired_speed, double *output_pos, int n)
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

void add_simulation_implementation(simul_t f)
{
    sim_counter++;
    sim_list = realloc(sim_list, sim_counter * sizeof(simul_t));
    sim_list[sim_counter - 1] = f;
}

void add_direction_implementation(void (*f)(double *, double *, double *, int))
{
    counter[2 * TDIR + 1]++;
    direction_ptr_v = realloc(direction_ptr_v, counter[2 * TDIR + 1] * sizeof(void (*)(double *, double *, double *, int)));
    direction_ptr_v[counter[2 * TDIR + 1] - 1] = f;
}

void add_acceleration_implementation(void (*f)(double *, double *, double *, double *, int))
{
    counter[2 * TACC + 1]++;
    acceleration_ptr_v = realloc(acceleration_ptr_v, counter[2 * TACC + 1] * sizeof(void (*)(double *, double *, double *, double *, int)));
    acceleration_ptr_v[counter[2 * TACC + 1] - 1] = f;
}

void add_social_implementation(void (*f)(double *, double *, double *, double *, int, int))
{
    counter[2 * TFRC + 1]++;
    social_ptr_v = realloc(social_ptr_v, counter[2 * TFRC + 1] * sizeof(void (*)(double *, double *, double *, double *, int, int)));
    social_ptr_v[counter[2 * TFRC + 1] - 1] = f;
}

void add_pos_implementation(void (*f)(double *, double *, double *, double *, double *, double *, int))
{
    counter[2 * TPOS + 1]++;
    pos_ptr_v = realloc(pos_ptr_v, counter[2 * TPOS + 1] * sizeof(void (*)(double *, double *, double *, double *, double *, double *, int)));
    pos_ptr_v[counter[2 * TPOS + 1] - 1] = f;
}

void copy_init(double *s_pos, double *s_dir, double *s_fdes, double *s_bor, double *s_spe,
               double **pos, double **dir, double **fdes, double **bor, double **spe, int n)
{
    *pos = malloc(n * 2 * sizeof(double));
    *dir = malloc(n * 2 * sizeof(double));
    *fdes = malloc(n * 2 * sizeof(double));
    *bor = malloc(N_BORDERS * sizeof(double));
    *spe = malloc(n * sizeof(double));

    memcpy(*pos, s_pos, n * 2 * sizeof(double));   //
    memcpy(*dir, s_dir, n * 2 * sizeof(double));   //
    memcpy(*fdes, s_fdes, n * 2 * sizeof(double)); //
    memcpy(*bor, s_bor, N_BORDERS * sizeof(double));
    memcpy(*spe, s_spe, n * sizeof(double));
}

void allocate_arrays(double **spe, double **vel, double **acc, double **prep, double **brep,
                     double **frc, int n)
{
    *spe = (double *)calloc(n, sizeof(double));
    *vel = (double *)calloc(n * 2, sizeof(double));
    *acc = (double *)calloc(n * 2, sizeof(double));
    *prep = (double *)calloc(n * n * 2, sizeof(double));
    *brep = (double *)calloc(n * N_BORDERS * 2, sizeof(double));
    *frc = (double *)calloc(n * 2, sizeof(double));
}


/* finite-differences functions */

double compute_people_repulsion_fd(double *position, double *parameters, double *desired_direction, double *actual_speed, int i, int j)
{
    double rx_ab = parameters[0];
    double ry_ab = parameters[1];
    double ex_a = desired_direction[i * 2];
    double ey_a = desired_direction[i * 2 + 1];
    double ex_b = desired_direction[j * 2];
    double ey_b = desired_direction[j * 2 + 1];
    double vb = actual_speed[j];
    double delta_b = vb * TIMESTEP;

    double r_ab_norm = sqrt(rx_ab * rx_ab + ry_ab * ry_ab);

    double rx_ab_mex = rx_ab - delta_b * ex_b;
    double ry_ab_mey = ry_ab - delta_b * ey_b;

    double r_ab_me_norm = sqrt(rx_ab_mex * rx_ab_mex + ry_ab_mey * ry_ab_mey);

    double b = sqrt((r_ab_norm + r_ab_me_norm) * (r_ab_norm + r_ab_me_norm) - (delta_b * delta_b)) / 2;

    return V_ALPHA_BETA * exp(-b / SIGMA);
}

void add_estimated_gradient_people_repulsion(double *gradient, double *parameters, int n_parameters, double *position, double *desired_direction, double *actual_speed, int i, int j)
{

    double dp = 1.0e-7;
    double f_P, f_M;

    for (int p = 0; p < n_parameters; p++)
    {
        double tmpVal = parameters[p];
        parameters[p] = tmpVal + dp;
        f_P = compute_people_repulsion_fd(position, parameters, desired_direction, actual_speed, i, j);

        parameters[p] = tmpVal - dp;
        f_M = compute_people_repulsion_fd(position, parameters, desired_direction, actual_speed, i, j);

        gradient[p] = (f_P - f_M) / (2 * dp);
    }
    double ex_a = desired_direction[i * 2];
    double ey_a = desired_direction[i * 2 + 1];
    double check = ex_a * (gradient[0]) + ey_a * (gradient[1]);
    double threshold = sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]) * cos(PSI);
    double w = check >= threshold ? 1 : INFLUENCE;
    gradient[0] = w * -gradient[0];
    gradient[1] = w * -gradient[1];
}

void test_people_repulsion_with_FD(double *people_repulsion_term, int n, double *position, double *desired_direction, double *actual_speed)
{
    double tol = 1e-4;
    double eps = 1e-10;

    double fd_gradient[2];
    double parameters[2];
    double *analytic_gradient = people_repulsion_term;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            double rx_ab = position[i * 2] - position[j * 2];
            double ry_ab = position[i * 2 + 1] - position[j * 2 + 1];
            parameters[0] = rx_ab;
            parameters[1] = ry_ab;
            add_estimated_gradient_people_repulsion(fd_gradient, parameters, 2, position, desired_direction, actual_speed, i, j);
            for (int p = 0; p < 2; p++)
            {
                double absErr = fabs(fd_gradient[p] - analytic_gradient[i * (2 * n) + 2 * j + p]);
                double relError = 2 * absErr / (eps + analytic_gradient[i * (2 * n) + 2 * j + p] + fd_gradient[p]);

                if (relError > tol && absErr > 1e-6)
                {
                    printf("Mismatch in people_repulsion element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
                    printf("Mismatch in people_repulsion element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
                }
            }
        }
    }
}
// default values

double compute_border_repulsion_fd(double *position, double *parameters, double *borders, int i, int j)
{
    double rx_a = parameters[0];
    double ry_a = parameters[1];

    double rx_aB = 0.0;
    double ry_aB = ry_a - borders[j];

    double r_aB_norm = fabs(ry_aB);

    return U_ALPHA_B * exp(-r_aB_norm / R);
}

void add_estimated_gradient_border_repulsion(double *gradient, double *parameters, int n_parameters, double *position, double *borders, int i, int j)
{

    double dp = 1.0e-7;
    double f_P, f_M;

    for (int p = 0; p < n_parameters; p++)
    {
        double tmpVal = parameters[p];
        parameters[p] = tmpVal + dp;
        f_P = compute_border_repulsion_fd(position, parameters, borders, i, j);

        parameters[p] = tmpVal - dp;
        f_M = compute_border_repulsion_fd(position, parameters, borders, i, j);

        gradient[p] = (f_P - f_M) / (2 * dp);
    }
    gradient[0] *= -1;
    gradient[1] *= -1;
}

void test_border_repulsion_with_FD(double *border_repulsion_term, double *position, double *borders, int n_borders, int n)
{
    double tol = 1e-4;
    double eps = 1e-10;

    double fd_gradient[2];
    double parameters[2];
    double *analytic_gradient = border_repulsion_term;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n_borders; j++)
        {
            double rx_a = position[i * 2];
            double ry_a = position[i * 2 + 1];
            double rx_aB = 0.0;
            double ry_aB = ry_a - borders[j];

            parameters[0] = rx_aB;
            parameters[1] = ry_aB;
            add_estimated_gradient_border_repulsion(fd_gradient, parameters, 2, position, borders, i, j);
            for (int p = 0; p < 2; p++)
            {
                double absErr = fabs(fd_gradient[p] - analytic_gradient[i * (2 * n_borders) + 2 * j + p]);
                double relError = 2 * absErr / (eps + analytic_gradient[i * (2 * n_borders) + 2 * j + p] + fd_gradient[p]);

                if (relError > tol && absErr > 1e-6)
                {
                    printf("Mismatch in border_repulsion: element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n_borders) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
                    printf("Mismatch in border_repulsion: element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n_borders) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
                }
            }
        }
    }
}