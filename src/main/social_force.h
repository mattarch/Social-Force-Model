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
#define V_ALPHA_BETA 2.1           // in m^{2}s^{-2}
#define SIGMA 0.3                  // in m
#define U_ALPHA_B 10.0             // in m^{2}s^{-2}
#define R 0.2                      // in m
#define PSI 1.75                   // in radians
#define PROJECTION_FACTOR -0.17825 // == cos(PSI), pure
#define INFLUENCE 0.5              // pure
#define DIV_FACTOR 1.75            // == V_ALPHA_BETA / 4 / SIGMA , pure

// benchmark parameters
#define REP 15
#define CYCLES_REQUIRED 1e8

#define NTESTS_FINITE_DIFFERENCES 10

#ifdef DEBUG
#define CONSOLE_PRINT(x) printf x
#else
#define CONSOLE_PRINT(x) \
  do                     \
  {                      \
  } while (0)
#endif

#define IndexX(i) ((2 * i))        //(i)
#define IndexY(i, n) ((2 * i + 1)) //((n + i))

#define IndexX_matrix(i, j, n) ((i * (2 * n) + 2 * j))     // ((i * n + j))
#define IndexY_matrix(i, j, n) ((i * (2 * n) + 2 * j + 1)) //((n * n + i * n + j))

#define IndexX_border(i, j, n) ((i * (2 * 2) + 2 * j))     //((j * 2 * n + i))
#define IndexY_border(i, j, n) ((i * (2 * 2) + 2 * j + 1)) //((n + j * 2 * n + i))

#define IndexX_border_old(i, j, n) ((i * (2 * 2) + 2 * j))
#define IndexY_border_old(i, j, n) ((i * (2 * 2) + 2 * j + 1))

#define IndexX_old(i) ((2 * i))
#define IndexY_old(i, n) ((2 * i + 1))

// typedefs
typedef void (*sim_func)(int, int, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
typedef long long unsigned int (*flops_func)(int);

typedef struct
{
  sim_func f;
  flops_func flops_f;
  char *name;
} sim_t;

void add_implementations(sim_t **sim_list, int *sim_counter, sim_func *test_functions_list, int *test_func_counter);
void initialize_people(double *position, double *desired_direction, double *final_destination, double *desired_speed, int n);
void compute_max_speed(double *desired_speed, double *desired_max_speed, int n);
void initialize_borders(double *borders, int n_borders);
void run_bench(sim_t sim);
int compare(const void *a, const void *b);
long long unsigned int compute_basic_flops(int number_of_people);
long long unsigned int compute_simplified_flops(int number_of_people);
void add_function(sim_t **sim_list, int *sim_counter, sim_func f, flops_func flops_f, char *name);
void add_test_function(sim_func *test_functions_list, sim_func f, int *test_func_counter);

#endif