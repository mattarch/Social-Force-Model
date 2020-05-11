#ifndef TEST_SET_DOUBLE_H
#define TEST_SET_DOUBLE_H

#define PRINT_RED printf("\033[0;31m")
#define PRINT_GREEN printf("\033[0;32m");
#define PRINT_DEF printf("\033[0m");

#define N_TESTS 4
#define TDIR 0
#define TACC 1
#define TFRC 2
#define TPOS 3

typedef struct T_double
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
} testcase_t_double;

void add_function_implementations_double();
void add_testcases_double();
int run_tests_double(sim_t **sim_list, int sim_counter);
int run_testcases_double();

void add_direction_testcase_double(char *name, double *pos, double *speed, double *output, int n);
void add_acceleration_testcase_double(char *name, double *dir, double *vel, double *spe, double *output, int n);
void add_compute_social_force_testcase_double(char *name, double *acc, double *people_rep, double *border_rep, double *output, int n, int n_borders);
void add_position_testcase_double(char *name, double *pos, double *dir, double *speed, double *force, double *actual_vel, double *desired_speed, double *output_pos, int n);
void add_direction_implementation_double(void (*f)(double *, double *, double *, int));
void add_acceleration_implementation_double(void (*f)(double *, double *, double *, double *, int));
void add_social_implementation_double(void (*f)(double *, double *, double *, double *, int, int));
void add_pos_implementation_double(void (*f)(double *, double *, double *, double *, double *, double *, int));

#endif