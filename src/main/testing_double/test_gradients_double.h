#ifndef TEST_GRADIENTS_DOUBLE_H
#define TEST_GRADIENTS_DOUBLE_H

#include "../social_force.h"
#include "../parse_args.h"

/* finite-differences functions */
void run_finite_differences_double(sim_func_double f, struct arguments arguments);

double compute_people_repulsion_fd_double(double *position, double *parameters, double *desired_direction, double *actual_speed, int i, int j, int n);
void add_estimated_gradient_people_repulsion_double(double *gradient, double *parameters, int n_parameters, double *position, double *desired_direction, double *actual_speed, int i, int j, int n);
void test_people_repulsion_with_FD_double(double *people_repulsion_term, int n, double *position, double *desired_direction, double *actual_speed);
double compute_border_repulsion_fd_double(double *position, double *parameters, double *borders, int i, int j, int n);
void add_estimated_gradient_border_repulsion_double(double *gradient, double *parameters, int n_parameters, double *position, double *borders, int i, int j, int n);
void test_border_repulsion_with_FD_double(double *border_repulsion_term, double *position, double *borders, int n_borders, int n);

#endif