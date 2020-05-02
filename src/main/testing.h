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
#define EPS_NEW 1e-01
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
int check_square_distance(float *expected, float *res, int n, int case_n);
void copy_init(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
			   float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n);
void copy_init_new(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
				   float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n);
void copy_state(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
				float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n);
void copy_state_new(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
				float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n);
void allocate_arrays(float **spe, float **vel, float **acc, float **prep, float **brep,
					 float **frc, int n);
void add_direction_testcase(char *name, float *pos, float *speed, float *output, int n);
void add_acceleration_testcase(char *name, float *dir, float *vel, float *spe, float *output, int n);
void add_compute_social_force_testcase(char *name, float *acc, float *people_rep, float *border_rep, float *output, int n, int n_borders);
void add_position_testcase(char *name, float *pos, float *dir, float *speed, float *force, float *actual_vel, float *desired_speed, float *output_pos, int n);
void add_direction_implementation(void (*f)(float *, float *, float *, int));
void add_acceleration_implementation(void (*f)(float *, float *, float *, float *, int));
void add_social_implementation(void (*f)(float *, float *, float *, float *, int, int));
void add_pos_implementation(void (*f)(float *, float *, float *, float *, float *, float *, int));

/* finite-differences functions */

float compute_people_repulsion_fd(float *position, float *parameters, float *desired_direction, float *actual_speed, int i, int j);
void add_estimated_gradient_people_repulsion(float *gradient, float *parameters, int n_parameters, float *position, float *desired_direction, float *actual_speed, int i, int j);
void test_people_repulsion_with_FD(float *people_repulsion_term, int n, float *position, float *desired_direction, float *actual_speed);
float compute_border_repulsion_fd(float *position, float *parameters, float *borders, int i, int j);
void add_estimated_gradient_border_repulsion(float *gradient, float *parameters, int n_parameters, float *position, float *borders, int i, int j);
void test_border_repulsion_with_FD(float *border_repulsion_term, float *position, float *borders, int n_borders, int n);

typedef struct T
{
	char *name;
	int n;
	int n_borders;
	float *dir;
	float *pos;
	float *spe;
	float *dspe;
	float *vel;
	float *acc;
	float *pre;
	float *bre;
	float *bor;
	float *frc;
	float *des;
} testcase_t;

#endif
