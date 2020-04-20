/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef TESTING_H_ /* Include guard */
#define TESTING_H_

#include "parse_args.h"

#define PRINT_RED printf("\033[0;31m")
#define PRINT_GREEN printf("\033[0;32m");
#define PRINT_DEF printf("\033[0m");

#define EPS 1e-03
#define N_TESTS 4
#define TDIR 0
#define TACC 1
#define TFRC 2
#define TPOS 3

void add_function_implementations();
void add_testcases();
int run_tests(sim_t **sim_list, int sim_counter);
void run_finite_differences(sim_func f, struct arguments arguments);
int run_testcases();
int compare_simulations(sim_t **sim_list, int sim_counter);
int check_square_distance(double *expected, double *res, int n);
void copy_init(double *s_pos, double *s_dir, double *s_fdes, double *s_bor, double *s_spe,
			   double **pos, double **dir, double **fdes, double **bor, double **spe, int n);
void allocate_arrays(double **spe, double **vel, double **acc, double **prep, double **brep,
					 double **frc, int n);
void add_direction_testcase(char *name, double *pos, double *speed, double *output, int n);
void add_acceleration_testcase(char *name, double *dir, double *vel, double *spe, double *output, int n);
void add_compute_social_force_testcase(char *name, double *acc, double *people_rep, double *border_rep, double *output, int n, int n_borders);
void add_position_testcase(char *name, double *pos, double *dir, double *speed, double *force, double *actual_vel, double *desired_speed, double *output_pos, int n);
void add_direction_implementation(void (*f)(double *, double *, double *, int));
void add_acceleration_implementation(void (*f)(double *, double *, double *, double *, int));
void add_social_implementation(void (*f)(double *, double *, double *, double *, int, int));
void add_pos_implementation(void (*f)(double *, double *, double *, double *, double *, double *, int));

/* finite-differences functions */

double compute_people_repulsion_fd(double *position, double *parameters, double *desired_direction, double *actual_speed, int i, int j);
void add_estimated_gradient_people_repulsion(double *gradient, double *parameters, int n_parameters, double *position, double *desired_direction, double *actual_speed, int i, int j);
void test_people_repulsion_with_FD(double *people_repulsion_term, int n, double *position, double *desired_direction, double *actual_speed);
double compute_border_repulsion_fd(double *position, double *parameters, double *borders, int i, int j);
void add_estimated_gradient_border_repulsion(double *gradient, double *parameters, int n_parameters, double *position, double *borders, int i, int j);
void test_border_repulsion_with_FD(double *border_repulsion_term, double *position, double *borders, int n_borders, int n);

typedef struct T
{
	char *name;
	int n;
	int n_borders;
	double *dir;
	double *pos;
	double *spe;
	double *dspe;
	double *vel;
	double *acc;
	double *pre;
	double *bre;
	double *bor;
	double *frc;
	double *des;
} testcase_t;

#endif
