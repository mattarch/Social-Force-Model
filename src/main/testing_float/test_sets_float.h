#ifndef TEST_SET_FLOAT_H
#define TEST_SET_FLOAT_H

#define PRINT_RED printf("\033[0;31m")
#define PRINT_GREEN printf("\033[0;32m");
#define PRINT_DEF printf("\033[0m");

#define N_TESTS 4
#define TDIR 0
#define TACC 1
#define TFRC 2
#define TPOS 3

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

void add_function_implementations();
void add_testcases();
int run_tests(sim_t **sim_list, int sim_counter);
int run_testcases();

void add_direction_testcase(char *name, float *pos, float *speed, float *output, int n);
void add_acceleration_testcase(char *name, float *dir, float *vel, float *spe, float *output, int n);
void add_compute_social_force_testcase(char *name, float *acc, float *people_rep, float *border_rep, float *output, int n, int n_borders);
void add_position_testcase(char *name, float *pos, float *dir, float *speed, float *force, float *actual_vel, float *desired_speed, float *output_pos, int n);
void add_direction_implementation(void (*f)(float *, float *, float *, int));
void add_acceleration_implementation(void (*f)(float *, float *, float *, float *, int));
void add_social_implementation(void (*f)(float *, float *, float *, float *, int, int));
void add_pos_implementation(void (*f)(float *, float *, float *, float *, float *, float *, int));

#endif