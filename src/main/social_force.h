/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef SOCIAL_FORCE_H_ /* Include guard */
#define SOCIAL_FORCE_H_

// default values
#define AVG_SPEED 1.34 // in basic scenario this is the desired speed of every person
#define RELAX_TIME 0.5
#define INV_RELAX_TIME 2 // == 1 / REALX_TIME in [seconds]; used for smooth acceleration
#define N_BORDERS 2      // number of biarders used in the scenario
#define TIMESTEP 0.2     // in [seconds]; timestep for simulation

// parameters model PAGE 8
#define V_ALPHA_BETA 2.1 // in m^{2}s^{-2}
#define SIGMA 0.3        // in m
#define INV_SIGMA 3.333333333
#define U_ALPHA_B 10.0 // in m^{2}s^{-2}
#define UTIMESR 50.0   // U_ALPHA_B * INV_R
#define R 0.2          // in m
#define INV_R 5
#define PSI 1.75                   // in radians
#define PROJECTION_FACTOR -0.17825 // == cos(PSI), pure
#define INFLUENCE 0.5              // pure
#define DIV_FACTOR 1.75            // == V_ALPHA_BETA / 4 / SIGMA , pure

// benchmark parameters
#define REP 15
#define CYCLES_REQUIRED 1e9

#define NTESTS_FINITE_DIFFERENCES 10

#define EPS 1e-02
#define EPS_NEW 1e-02

#ifdef DEBUG
#define CONSOLE_PRINT(x) printf x
#else
#define CONSOLE_PRINT(x) \
  do                     \
  {                      \
  } while (0)
#endif

#define IndexX(i) ((i))
#define IndexY(i, n) ((n + i))

#define IndexX_old(i) ((2 * i))
#define IndexY_old(i, n) ((2 * i + 1))

#define IndexX_matrix(i, j, n) ((i * n + j))
#define IndexY_matrix(i, j, n) ((n * n + i * n + j))

#define IndexX_matrix_old(i, j, n) ((i * (2 * n) + 2 * j))
#define IndexY_matrix_old(i, j, n) ((i * (2 * n) + 2 * j + 1))

#define IndexX_border(i, j, n) ((j * 2 * n + i))
#define IndexY_border(i, j, n) ((n + j * 2 * n + i))

#define IndexX_border_old(i, j, n) ((i * (2 * 2) + 2 * j))
#define IndexY_border_old(i, j, n) ((i * (2 * 2) + 2 * j + 1))

#define IS_FLOAT 0
#define IS_DOUBLE 1

// typedefs
typedef void (*sim_func)(int, int, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *);
typedef void (*sim_func_double)(int, int, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

typedef long long unsigned int (*flops_func)(int);
typedef double (*op_int_func)(int, int);

typedef struct
{
  sim_func f;
  sim_func_double f_double;
  flops_func flops_f;
  op_int_func op_f;
  int is_double;
  char *name;
} sim_t;

void add_implementations_float(sim_t **sim_list, int *sim_counter, sim_t **test_functions_list, int *test_func_counter);
void add_implementations_double(sim_t **sim_list, int *sim_counter, sim_t **test_functions_list, int *test_func_counter);
void add_implementations_double_restructured(sim_t **sim_list, int *sim_counter, sim_t **test_functions_list, int *test_func_counter);

void add_function(sim_t **sim_list, int *sim_counter, sim_func f, sim_func_double f_double, flops_func flops_f, int is_double, op_int_func op_f, char *name);
void add_test_function(sim_t **test_functions_list, sim_func f, sim_func_double f_double, int is_double, int *test_func_counter);

/* initialization */
void initialize_people_float(float *position, float *desired_direction, float *final_destination, float *desired_speed, int n);
void initialize_max_speed_float(float *desired_speed, float *desired_max_speed, int n);
void initialize_borders_float(float *borders, int n_borders);

void initialize_people_double(double *position, double *desired_direction, double *final_destination, double *desired_speed, int n);
void initialize_max_speed_double(double *desired_speed, double *desired_max_speed, int n);
void initialize_borders_double(double *borders, int n_borders);

/* benchmarking */
void run_bench_float(sim_t sim);
void run_bench_double(sim_t sim);
int compare(const void *a, const void *b);
long long unsigned int compute_simplified_flops(int number_of_people);
long long unsigned compute_251_flops(int number_of_people);
double compute_operational_intensity_251(int number_of_people, int number_of_timesteps);
// double compute_operational_intensity_0(int n, int number_of_timesteps);
double compute_operational_intensity_float(int n, int number_of_timesteps);
double compute_operational_intensity_optimized_float(int n, int number_of_timesteps);
double compute_operational_intensity_double(int n, int number_of_timesteps);
double compute_operational_intensity_optimized_double(int n, int number_of_timesteps);
long long unsigned lower_memory_bound_float(int number_of_people);
long long unsigned lower_memory_bound_optimized_float(int number_of_people);
long long unsigned lower_memory_bound_double(int number_of_people);
long long unsigned lower_memory_bound_optimized_double(int number_of_people);
unsigned int flush_cache(int n);

#endif