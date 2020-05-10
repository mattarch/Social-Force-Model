

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../social_force.h"
#include "../aligned_malloc.h"
#include "../aligned_free.h"
#include "../utility.h"

#include "test_gradients_double.h"
#include "compare_simulations_double.h"

/*
* Runs all the function in the test list. These functions are used to check that the
* gradient computed analytically are correct by comparing the result to a finite differences
* implementation of the gradients.
*/
void run_finite_differences_double(sim_func_double f, struct arguments arguments)
{
    int number_of_people = arguments.n_people;
    int n_timesteps = arguments.n_timesteps;
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
    double *people_repulsion_term = (double *)aligned_malloc(number_of_people * number_of_people * 2 * sizeof(double), 32);
    set_zero_double(people_repulsion_term, number_of_people * number_of_people * 2);
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

/* central-differences functions */

double compute_people_repulsion_fd_double(double *position, double *parameters, double *desired_direction, double *actual_speed, int i, int j, int n)
{
    double rx_ab = parameters[0];
    double ry_ab = parameters[1];
    double ex_a = desired_direction[IndexX(i)];
    double ey_a = desired_direction[IndexY(i, n)];
    double ex_b = desired_direction[IndexX(j)];
    double ey_b = desired_direction[IndexY(j, n)];
    double vb = actual_speed[j];
    double delta_b = vb * TIMESTEP;

    double r_ab_norm = sqrt(rx_ab * rx_ab + ry_ab * ry_ab);

    double rx_ab_mex = rx_ab - delta_b * ex_b;
    double ry_ab_mey = ry_ab - delta_b * ey_b;

    double r_ab_me_norm = sqrt(rx_ab_mex * rx_ab_mex + ry_ab_mey * ry_ab_mey);

    double b = sqrt((r_ab_norm + r_ab_me_norm) * (r_ab_norm + r_ab_me_norm) - (delta_b * delta_b)) / 2;

    return V_ALPHA_BETA * exp(-b / SIGMA);
}

void add_estimated_gradient_people_repulsion_double(double *gradient, double *parameters, int n_parameters, double *position, double *desired_direction, double *actual_speed, int i, int j, int n)
{
    double dp = 1.0e-7;
    double f_P, f_M;

    for (int p = 0; p < n_parameters; p++)
    {
        double tmpVal = parameters[p];
        parameters[p] = tmpVal + dp;
        f_P = compute_people_repulsion_fd_double(position, parameters, desired_direction, actual_speed, i, j, n);

        parameters[p] = tmpVal - dp;
        f_M = compute_people_repulsion_fd_double(position, parameters, desired_direction, actual_speed, i, j, n);

        gradient[p] = (f_P - f_M) / (2 * dp);
    }
    double ex_a = desired_direction[IndexX(i)];
    double ey_a = desired_direction[IndexY(i, n)];
    double check = ex_a * (gradient[0]) + ey_a * (gradient[1]);
    double threshold = sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]) * cos(PSI);
    double w = check >= threshold ? 1 : INFLUENCE;
    gradient[0] = w * -gradient[0];
    gradient[1] = w * -gradient[1];
}

void test_people_repulsion_with_FD_double(double *people_repulsion_term, int n, double *position, double *desired_direction, double *actual_speed)
{
    int num_errors = 0;
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
            double rx_ab = position[IndexX(i)] - position[IndexX(j)];
            double ry_ab = position[IndexY(i, n)] - position[IndexY(j, n)];
            parameters[0] = rx_ab;
            parameters[1] = ry_ab;
            add_estimated_gradient_people_repulsion_double(fd_gradient, parameters, 2, position, desired_direction, actual_speed, i, j, n);

            double absErr = fabs(fd_gradient[0] - analytic_gradient[IndexX_matrix(i, j, n) ]);
            double relError = 2 * absErr / (eps + analytic_gradient[IndexX_matrix(i, j, n) ] + fd_gradient[0]);

            if (relError > tol && absErr > EPS)
            {
                printf("Mismatch in people_repulsion element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[IndexX_matrix(i, j, n)], fd_gradient[0], absErr, relError * 100);
                printf("Mismatch in people_repulsion element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[IndexX_matrix(i, j, n)], fd_gradient[0], absErr, relError * 100);
                num_errors++;
                if (num_errors > 51)
                {
                    return;
                }
            }

            absErr = fabs(fd_gradient[1] - analytic_gradient[IndexY_matrix(i, j, n)]);
            relError = 2 * absErr / (eps + analytic_gradient[IndexY_matrix(i, j, n)] + fd_gradient[1]);

            if (relError > tol && absErr > EPS)
            {
                printf("Mismatch in people_repulsion element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[IndexY_matrix(i, j, n)], fd_gradient[1], absErr, relError * 100);
                printf("Mismatch in people_repulsion element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[IndexY_matrix(i, j, n)], fd_gradient[1], absErr, relError * 100);
                num_errors++;
                if (num_errors > 51)
                {
                    return;
                }
            }
        }
    }
}
// default values

double compute_border_repulsion_fd_double(double *position, double *parameters, double *borders, int i, int j, int n)
{
    double rx_a = parameters[0];
    double ry_a = parameters[1];

    double rx_aB = 0.0;
    double ry_aB = ry_a - borders[j];

    double r_aB_norm = fabs(ry_aB);

    return U_ALPHA_B * exp(-r_aB_norm / R);
}

void add_estimated_gradient_border_repulsion_double(double *gradient, double *parameters, int n_parameters, double *position, double *borders, int i, int j, int n)
{
    double dp = 1.0e-7;
    double f_P, f_M;

    for (int p = 0; p < n_parameters; p++)
    {
        double tmpVal = parameters[p];
        parameters[p] = tmpVal + dp;
        f_P = compute_border_repulsion_fd_double(position, parameters, borders, i, j, n);

        parameters[p] = tmpVal - dp;
        f_M = compute_border_repulsion_fd_double(position, parameters, borders, i, j, n);

        gradient[p] = (f_P - f_M) / (2 * dp);
    }
    gradient[0] *= -1;
    gradient[1] *= -1;
}

void test_border_repulsion_with_FD_double(double *border_repulsion_term, double *position, double *borders, int n_borders, int n)
{
    int num_errors = 0;
    double tol = 1e-4;
    double eps = 1e-10;

    double fd_gradient[2];
    double parameters[2];
    double *analytic_gradient = border_repulsion_term;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n_borders; j++)
        {
            double rx_a = position[IndexX(i)];
            double ry_a = position[IndexY(i, n)];
            double rx_aB = 0.0;
            double ry_aB = ry_a - borders[j];

            parameters[0] = rx_aB;
            parameters[1] = ry_aB;
            add_estimated_gradient_border_repulsion_double(fd_gradient, parameters, 2, position, borders, i, j, n);
            for (int p = 0; p < 1; p++)
            {
                double absErr = fabs(fd_gradient[p] - analytic_gradient[IndexX_border(i, j, n) + p]);
                double relError = 2 * absErr / (eps + analytic_gradient[IndexX_border(i, j, n) + p] + fd_gradient[p]);

                if (relError > tol && absErr > EPS)
                {
                    printf("Mismatch in border_repulsion: element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[IndexX_border(i, j, n) + p], fd_gradient[p], absErr, relError * 100);
                    printf("Mismatch in border_repulsion: element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[IndexX_border(i, j, n) + p], fd_gradient[p], absErr, relError * 100);
                    num_errors++;
                    if (num_errors > 51)
                    {
                        return;
                    }
                }
            }
        }
    }
}