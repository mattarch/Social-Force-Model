#ifndef COMPARE_SIMULATIONS_FLOAT_H
#define COMPARE_SIMULATIONS_FLOAT_H

#include "../parse_args.h"
#include "../social_force.h"

int compare_simulations(sim_t **sim_list, int sim_counter);
int check_absolute_distance(float *expected, float *res, int n, int case_n);
void copy_init(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
               float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n);
void copy_state(float *s_pos, float *s_dir, float *s_fdes, float *s_bor, float *s_spe, float *s_mspe,
                float **pos, float **dir, float **fdes, float **bor, float **spe, float **mspe, int n,sim_t sim);
void allocate_arrays(float **spe, float **vel, float **acc, float **prep, float **brep,
                     float **frc, int n);
#endif