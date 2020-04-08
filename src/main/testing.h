#ifndef TESTING_H_ /* Include guard */
#define TESTING_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "social_force_model_basic.h"
#include "test_sets.h"

#define PRINT_RED printf("\033[0;31m")
#define PRINT_GREEN printf("\033[0;32m");
#define PRINT_DEF printf("\033[0m");

#define EPS 1e-06
#define N_TESTS 6
#define TDIR 0
#define TACC 1
#define TPREP 2
#define TBREP 3
#define TFRC 4
#define TPOS 5

void add_implementations();
void add_tests();
int run_tests();
int run();
int check_square_distance(double *expected, double *res, int n);
void add_direction_testcase(char *name, double *pos, double *speed, double *output, int n);
void add_acceleration_testcase(char *name, double *dir, double *vel, double *output, int n);
void add_people_repulsion_testcase(char *name, double *pos, double *speed, double *directions, double *output, int n);
void add_border_repulsion_testcase(char *name, double *pos, double *borders, double *output, int n, int n_borders);
void add_compute_social_force_testcase(char *name, double *acc, double *people_rep, double *border_rep, double *output, int n, int n_borders);
void add_position_testcase(char *name, double *pos, double *dir, double *speed, double *force, double *actual_vel, double *output_pos, int n);
void add_direction_implementation(void (*f)(double *, double *, double *, int));
void add_acceleration_implementation(void (*f)(double *, double *, double *, int));
void add_people_repulsion_implementation(void (*f)(double *, double *, double *, double *, int));
void add_border_repulsion_implementation(void (*f)(double *, double *, double *, int, int));
void add_social_implementation(void (*f)(double *, double *, double *, double *, int, int));
void add_pos_implementation(void (*f)(double *, double *, double *, double *, double *, int));

typedef struct T
{
	char *name;
	int n;
	int n_borders;
	double *dir;
	double *pos;
	double *spe;
	double *vel;
	double *acc;
	double *pre;
	double *bre;
	double *bor;
	double *frc;
	double *des;
} testcase_t;

testcase_t **testcases[N_TESTS];

int counter[2 * N_TESTS];

void (**direction_ptr_v)(double *, double *, double *, int);
void (**acceleration_ptr_v)(double *, double *, double *, int);
void (**repulsion_ptr_v)(double *, double *, double *, double *, int);
void (**border_repulsion_ptr_v)(double *, double *, double *, int, int);
void (**social_ptr_v)(double *, double *, double *, double *, int, int);
void (**pos_ptr_v)(double *, double *, double *, double *, double *, int);

void add_tests()
{
	//direction testcases
	//add_direction_testcase("test testing", direction_position0, direction_fdest0, direction_expected0, direction_n0);
	//add_direction_testcase("error testing", direction_position0, direction_fdest0, direction_expected1, direction_n0);

	//acceleration testcases
	add_acceleration_testcase("acceleration_test_from_start_straight", acceleration_direction0, acceleration_vel0, acceleration_expected0, acceleration_n0);
	add_acceleration_testcase("acceleration_test_from_start_45", acceleration_direction1, acceleration_vel1, acceleration_expected1, acceleration_n1);
	add_acceleration_testcase("acceleration_test_with_velocity_straight", acceleration_direction2, acceleration_vel2, acceleration_expected2, acceleration_n2);
	add_acceleration_testcase("acceleration_test_with_velocity_45", acceleration_direction3, acceleration_vel3, acceleration_expected3, acceleration_n3);
	add_acceleration_testcase("acceleration_test_to_0", acceleration_direction4, acceleration_vel4, acceleration_expected4, acceleration_n4);
	add_acceleration_testcase("deacceleration_test_straight", acceleration_direction5, acceleration_vel5, acceleration_expected5, acceleration_n5);
	add_acceleration_testcase("deacceleration_test_45", acceleration_direction6, acceleration_vel6, acceleration_expected6, acceleration_n6);

	//people repulsion testcases
	add_people_repulsion_testcase("test1", repulsion_position0, repulsion_speed0, repulsion_direction0, repulsion_expected_result0, repulsion_n0);

	//border repulsion
	//add_border_repulsion_testcase("segfault test", border_pos0, border_borders0, border_expected0, border_n0, border_nb0);

	//social force
	add_compute_social_force_testcase("basic test", social_acc0, social_prep0, social_brep0, social_expected0, social_n0, social_nb0);

	//position
	//add_position_testcase("segfault test", position_pos0, position_dir0, position_speed0, position_force0, position_vel0, position_expected0, position_n0);
}

void add_implementations()
{
	//update direction implementations
	//add_direction_implementation(update_desired_direction);

	//update acceleration implementations
	add_acceleration_implementation(update_acceleration_term);

	//update people repulsion implementations
	add_people_repulsion_implementation(update_people_repulsion_term);

	//update border repulsion implementations
	//add_border_repulsion_implementation(update_border_repulsion_term);

	//compute social force implementations
	//add_social_implementation(compute_social_force);

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
	add_implementations();

	//run tests
	int error = run();
	return error;
}

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
		case TPREP:
			id = "people repulsion";
			break;
		case TBREP:
			id = "border repulsion";
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
					void (*f1)(double *, double *, double *, int) = acceleration_ptr_v[i];
					result = (double *)calloc(result_size, sizeof(double));
					f1(cur->dir, result, cur->vel, np);
					control = cur->acc;
					break;
				case TPREP:
					result_size = np * np * 2;
					void (*f2)(double *, double *, double *, double *, int) = repulsion_ptr_v[i];
					result = (double *)calloc(result_size, sizeof(double));
					f2(cur->pos, cur->dir, cur->spe, result, np);
					control = cur->pre;
					break;
				case TBREP:
					result_size = np * cur->n_borders * 2;
					void (*f3)(double *, double *, double *, int, int) = border_repulsion_ptr_v[i];
					result = (double *)calloc(result_size, sizeof(double));
					f3(cur->pos, cur->bor, result, np, cur->n_borders);
					control = cur->bre;
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
					void (*f5)(double *, double *, double *, double *, double *, int) = pos_ptr_v[i];
					f5(cur->pos, cur->dir, cur->spe, cur->frc, cur->vel, np);
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

//returns 1 if error, returns 0 if OK
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

void add_acceleration_testcase(char *name, double *dir, double *vel, double *output, int n)
{
	counter[2 * TACC]++;
	testcases[TACC] = realloc(testcases[TACC], sizeof(testcase_t *) * counter[2 * TACC]);
	testcases[TACC][counter[2 * TACC] - 1] = malloc(sizeof(testcase_t));
	testcases[TACC][counter[2 * TACC] - 1]->name = name;
	testcases[TACC][counter[2 * TACC] - 1]->dir = dir;
	testcases[TACC][counter[2 * TACC] - 1]->vel = vel;
	testcases[TACC][counter[2 * TACC] - 1]->acc = output;
	testcases[TACC][counter[2 * TACC] - 1]->n = n;
}

void add_people_repulsion_testcase(char *name, double *pos, double *speed, double *directions, double *output, int n)
{
	counter[2 * TPREP]++;
	testcases[TPREP] = realloc(testcases[TPREP], sizeof(testcase_t *) * counter[2 * TPREP]);
	testcases[TPREP][counter[2 * TPREP] - 1] = malloc(sizeof(testcase_t));
	testcases[TPREP][counter[2 * TPREP] - 1]->name = name;
	testcases[TPREP][counter[2 * TPREP] - 1]->pos = pos;
	testcases[TPREP][counter[2 * TPREP] - 1]->spe = speed;
	testcases[TPREP][counter[2 * TPREP] - 1]->dir = directions;
	testcases[TPREP][counter[2 * TPREP] - 1]->pre = output;
	testcases[TPREP][counter[2 * TPREP] - 1]->n = n;
}

void add_border_repulsion_testcase(char *name, double *pos, double *borders, double *output, int n, int n_borders)
{
	counter[2 * TBREP]++;
	testcases[TBREP] = realloc(testcases[TBREP], counter[2 * TBREP] * sizeof(testcase_t *));
	testcases[TBREP][counter[2 * TBREP] - 1] = malloc(sizeof(testcase_t));
	testcases[TBREP][counter[2 * TBREP] - 1]->name = name;
	testcases[TBREP][counter[2 * TBREP] - 1]->pos = pos;
	testcases[TBREP][counter[2 * TBREP] - 1]->bor = borders;
	testcases[TBREP][counter[2 * TBREP] - 1]->bre = output;
	testcases[TBREP][counter[2 * TBREP] - 1]->n = n;
	testcases[TBREP][counter[2 * TBREP] - 1]->n_borders = n_borders;
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

void add_position_testcase(char *name, double *pos, double *dir, double *speed, double *force, double *actual_vel, double *output_pos, int n)
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
	testcases[TPOS][counter[2 * TPOS] - 1]->des = output_pos;
	testcases[TPOS][counter[2 * TPOS] - 1]->n = n;
}

void add_direction_implementation(void (*f)(double *, double *, double *, int))
{
	counter[2 * TDIR + 1]++;
	direction_ptr_v = realloc(direction_ptr_v, counter[2 * TDIR + 1] * sizeof(void (*)(double *, double *, double *, int)));
	direction_ptr_v[counter[2 * TDIR + 1] - 1] = f;
}

void add_acceleration_implementation(void (*f)(double *, double *, double *, int))
{
	counter[2 * TACC + 1]++;
	acceleration_ptr_v = realloc(acceleration_ptr_v, counter[2 * TACC + 1] * sizeof(void (*)(double *, double *, double *, int)));
	acceleration_ptr_v[counter[2 * TACC + 1] - 1] = f;
}

void add_people_repulsion_implementation(void (*f)(double *, double *, double *, double *, int))
{
	counter[2 * TPREP + 1]++;
	repulsion_ptr_v = realloc(repulsion_ptr_v, counter[2 * TPREP + 1] * sizeof(void (*)(double *, double *, double *, double *, int)));
	repulsion_ptr_v[counter[2 * TPREP + 1] - 1] = f;
}

void add_border_repulsion_implementation(void (*f)(double *, double *, double *, int, int))
{
	counter[2 * TBREP + 1]++;
	border_repulsion_ptr_v = realloc(border_repulsion_ptr_v, counter[2 * TBREP + 1] * sizeof(void (*)(double *, double *, double *, double *, int, int)));
	border_repulsion_ptr_v[counter[2 * TBREP + 1] - 1] = f;
}

void add_social_implementation(void (*f)(double *, double *, double *, double *, int, int))
{
	counter[2 * TFRC + 1]++;
	social_ptr_v = realloc(social_ptr_v, counter[2 * TFRC + 1] * sizeof(void (*)(double *, double *, double *, double *, int, int)));
	social_ptr_v[counter[2 * TFRC + 1] - 1] = f;
}

void add_pos_implementation(void (*f)(double *, double *, double *, double *, double *, int))
{
	counter[2 * TPOS + 1]++;
	pos_ptr_v = realloc(pos_ptr_v, counter[2 * TPOS + 1] * sizeof(void (*)(double *, double *, double *, double *, double *, int)));
	pos_ptr_v[counter[2 * TPOS + 1] - 1] = f;
}

double compute_people_repulsion_fd(double *position, double *parameters, double *desired_direction, double *actual_speed, int i, int j)
{
	double rx_ab = parameters[0];
	double ry_ab = parameters[1];
	double ex_a = desired_direction[i * 2];
	double ey_a = desired_direction[i * 2 + 1];
	double ex_b = desired_direction[j * 2];
	double ey_b = desired_direction[j * 2 + 1];
	double vb = actual_speed[j];
	double delta_b = vb * TIMESTEP; // 1 mul => 1 flops

	double r_ab_norm = sqrt(rx_ab * rx_ab + ry_ab * ry_ab); //(1) 1 sqrt, 1 add, 2 mul => 4 flops

	//me stands for "minus e"
	double rx_ab_mex = rx_ab - delta_b * ex_b; //1 add, 1 mul => 2 flops
	double ry_ab_mey = ry_ab - delta_b * ey_b; //1 add, 1 mul => 2 flops

	double r_ab_me_norm = sqrt(rx_ab_mex * rx_ab_mex + ry_ab_mey * ry_ab_mey); //(2)  1 sqrt, 2 mul, 1 add => 4 flops

	double b = sqrt((r_ab_norm + r_ab_me_norm) * (r_ab_norm + r_ab_me_norm) - (delta_b * delta_b)) / 2; //1 sqrt, 3 add, 2 mul, 1 div => 7 flops

	return V_ALPHA_BETA * exp(-b / SIGMA);
}

// only for 2 people
void add_estimated_gradient(double *gradient, double *parameters, int n_parameters, double *position, double *desired_direction, double *actual_speed, int i, int j)
{

	double dp = 1.0e-6;
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
	double check = ex_a * (-gradient[0]) + ey_a * (-gradient[1]);							   //2 mult, 1 add => 3 flops
	double threshold = sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]) * cos(PSI); //1 sqrt, 2 mults, 1 add => 4 flops
	double w = check >= threshold ? 1 : INFLUENCE;
	gradient[0] =  w * -gradient[0];
	gradient[1] =  w * -gradient[1];
}


void test_people_repulsion_with_FD(double *people_repulsion_term, int n, double *position, double *desired_direction, double *actual_speed)
{

	double tol = 1e-4;
	double eps = 1e-10;

	double fd_gradient[2]; //compute it on the fly
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
			add_estimated_gradient(fd_gradient, parameters, 2, position, desired_direction, actual_speed, i, j);
			for (int p = 0; p < 2; p++)
			{
				double absErr = fabs(fd_gradient[p] - analytic_gradient[i * (2 * n) + 2 * j + p]);
				double relError = 2 * absErr / (eps + analytic_gradient[i * (2 * n) + 2 * j + p] + fd_gradient[p]);
				//printf("%f, %f \n",fd_gradient[p],analytic_gradient[i * (2 * n) + 2 * j + p]);
				if (relError > tol && absErr > 1e-6)
				{
					printf("Mismatch element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
					printf("Mismatch element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf(%lf%%)\n", i, j, analytic_gradient[i * (2 * n) + 2 * j + p], fd_gradient[p], absErr, relError * 100);
				}
			}
		}
	}
}


#endif