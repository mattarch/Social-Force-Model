/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef COMPARE_SIMULATIONS_DOUBLE_H
#define COMPARE_SIMULATIONS_DOUBLE_H

#include "../parse_args.h"
#include "../social_force.h"

int compare_simulations_double(sim_t **sim_list, int sim_counter);
int check_absolute_distance_double(double *expected, double *res, int n, int case_n);
void copy_init_double(double *s_pos, double *s_dir, double *s_fdes, double *s_bor, double *s_spe, double *s_mspe,
               double **pos, double **dir, double **fdes, double **bor, double **spe, double **mspe, int n);
void copy_state_double(double *s_pos, double *s_dir, double *s_fdes, double *s_bor, double *s_spe, double *s_mspe,
                double **pos, double **dir, double **fdes, double **bor, double **spe, double **mspe, int n, sim_t sim);
void allocate_arrays_double(double **spe, double **vel, double **acc, double **prep, double **brep,
                     double **frc, int n);

#endif