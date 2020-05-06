/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#include <stdio.h>
#include <math.h>
#include <immintrin.h>
#include <assert.h>

#include "social_force_model_vectorize_2_5_1.h"
//#include "../tsc_x86.h"
#include "../social_force.h"
#include "../utility.h"

__m256 exp_fast_vec_2_5_1(__m256 x, __m256 one, __m256 exp_constant)
{
  x = _mm256_fmadd_ps(x, exp_constant, one);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  x = _mm256_mul_ps(x, x);
  return x;
}

void simulation_basic_vectorize_2_5_1(int number_of_people, int n_timesteps, float *position, float *speed, float *desired_direction, float *final_destination,
                                      float *borders, float *actual_velocity, float *acceleration_term, float *people_repulsion_term, float *border_repulsion_term,
                                      float *social_force, float *desired_speed, float *desired_max_speed)
{
  const float inv_sigma = 1 / SIGMA; // 1 div -> 1 flop
  const float cospsi = PROJECTION_FACTOR;
  int n = number_of_people;
  const float fb = borders[0]; // first border, this is positive and equal to walkway width
  const float sb = borders[1]; // second border, this is 0
  //myInt64 fma_cycles_1 = 0, sqrt_cycles_1 = 0;
  //myInt64 fma_cycles_2 = 0, sqrt_cycles_2 = 0;
  for (int step = 0; step < n_timesteps; step++)
  {
    //precomputation of the space walked in TIMESTEP [s], speed can be overwritten because it
    //is only used in the people repulsion term
    for (int j = 0; j < n - 3; j += 4)
    {
      speed[j] *= TIMESTEP; // 1 mult -> 1 flop
      speed[j + 1] *= TIMESTEP;
      speed[j + 2] *= TIMESTEP;
      speed[j + 3] *= TIMESTEP;
    } // n* (1 mult) -> n flops
    // iterate over every person
    for (int i = 0; i < n - 7; i += 8)
    {
      /************************************************/
      // LOADS
      /************************************************/
      float sfx0 = social_force[IndexX(i)]; //social force x
      float sfy0 = social_force[IndexY(i, n)];
      float sfx1 = social_force[IndexX(i + 1)];
      float sfy1 = social_force[IndexY(i + 1, n)];
      float sfx2 = social_force[IndexX(i + 2)];
      float sfy2 = social_force[IndexY(i + 2, n)];
      float sfx3 = social_force[IndexX(i + 3)];
      float sfy3 = social_force[IndexY(i + 3, n)];
      float sfx4 = social_force[IndexX(i + 4)];
      float sfy4 = social_force[IndexY(i + 4, n)];
      float sfx5 = social_force[IndexX(i + 5)];
      float sfy5 = social_force[IndexY(i + 5, n)];
      float sfx6 = social_force[IndexX(i + 6)];
      float sfy6 = social_force[IndexY(i + 6, n)];
      float sfx7 = social_force[IndexX(i + 7)];
      float sfy7 = social_force[IndexY(i + 7, n)];

      float rxa0 = position[IndexX(i)];
      float rya0 = position[IndexY(i, n)];
      float rxa1 = position[IndexX(i + 1)];
      float rya1 = position[IndexY(i + 1, n)];
      float rxa2 = position[IndexX(i + 2)];
      float rya2 = position[IndexY(i + 2, n)];
      float rxa3 = position[IndexX(i + 3)];
      float rya3 = position[IndexY(i + 3, n)];
      float rxa4 = position[IndexX(i + 4)];
      float rya4 = position[IndexY(i + 4, n)];
      float rxa5 = position[IndexX(i + 5)];
      float rya5 = position[IndexY(i + 5, n)];
      float rxa6 = position[IndexX(i + 6)];
      float rya6 = position[IndexY(i + 6, n)];
      float rxa7 = position[IndexX(i + 7)];
      float rya7 = position[IndexY(i + 7, n)];

      float exa0 = desired_direction[IndexX(i)];
      float eya0 = desired_direction[IndexY(i, n)];
      float exa1 = desired_direction[IndexX(i + 1)];
      float eya1 = desired_direction[IndexY(i + 1, n)];
      float exa2 = desired_direction[IndexX(i + 2)];
      float eya2 = desired_direction[IndexY(i + 2, n)];
      float exa3 = desired_direction[IndexX(i + 3)];
      float eya3 = desired_direction[IndexY(i + 3, n)];
      float exa4 = desired_direction[IndexX(i + 4)];
      float eya4 = desired_direction[IndexY(i + 4, n)];
      float exa5 = desired_direction[IndexX(i + 5)];
      float eya5 = desired_direction[IndexY(i + 5, n)];
      float exa6 = desired_direction[IndexX(i + 6)];
      float eya6 = desired_direction[IndexY(i + 6, n)];
      float exa7 = desired_direction[IndexX(i + 7)];
      float eya7 = desired_direction[IndexY(i + 7, n)];

      float avx0 = actual_velocity[IndexX(i)];
      float avy0 = actual_velocity[IndexY(i, n)];
      float avx1 = actual_velocity[IndexX(i + 1)];
      float avy1 = actual_velocity[IndexY(i + 1, n)];
      float avx2 = actual_velocity[IndexX(i + 2)];
      float avy2 = actual_velocity[IndexY(i + 2, n)];
      float avx3 = actual_velocity[IndexX(i + 3)];
      float avy3 = actual_velocity[IndexY(i + 3, n)];
      float avx4 = actual_velocity[IndexX(i + 4)];
      float avy4 = actual_velocity[IndexY(i + 4, n)];
      float avx5 = actual_velocity[IndexX(i + 5)];
      float avy5 = actual_velocity[IndexY(i + 5, n)];
      float avx6 = actual_velocity[IndexX(i + 6)];
      float avy6 = actual_velocity[IndexY(i + 6, n)];
      float avx7 = actual_velocity[IndexX(i + 7)];
      float avy7 = actual_velocity[IndexY(i + 7, n)];

      float dsv0 = desired_speed[i];     //desired speed value
      float dsv1 = desired_speed[i + 1]; //desired speed value
      float dsv2 = desired_speed[i + 2]; //desired speed value
      float dsv3 = desired_speed[i + 3]; //desired speed value
      float dsv4 = desired_speed[i + 4]; //desired speed value
      float dsv5 = desired_speed[i + 5]; //desired speed value
      float dsv6 = desired_speed[i + 6]; //desired speed value
      float dsv7 = desired_speed[i + 7]; //desired speed value

      float da0 = speed[i];
      float da1 = speed[i + 1];
      float da2 = speed[i + 2];
      float da3 = speed[i + 3];
      float da4 = speed[i + 4];
      float da5 = speed[i + 5];
      float da6 = speed[i + 6];
      float da7 = speed[i + 7];

      /************************************************/
      // UPDATE BORDER REPULSION TERM
      /************************************************/

      float ryaB00 = rya0 - fb; //1 add
      float ryaB10 = rya1 - fb; //1 add
      float ryaB20 = rya2 - fb; //1 add
      float ryaB30 = rya3 - fb; //1 add
      float ryaB40 = rya4 - fb; //1 add
      float ryaB50 = rya5 - fb; //1 add
      float ryaB60 = rya6 - fb; //1 add
      float ryaB70 = rya7 - fb; //1 add

      float se00 = exp_fast(ryaB00 * INV_R) * UTIMESR / (-ryaB00); //1 exp, 2 mult, 1 div
      float se10 = exp_fast(ryaB10 * INV_R) * UTIMESR / (-ryaB10); //1 exp, 2 mult, 1 div
      float se20 = exp_fast(ryaB20 * INV_R) * UTIMESR / (-ryaB20); //1 exp, 2 mult, 1 div
      float se30 = exp_fast(ryaB30 * INV_R) * UTIMESR / (-ryaB30); //1 exp, 2 mult, 1 div
      float se40 = exp_fast(ryaB40 * INV_R) * UTIMESR / (-ryaB40); //1 exp, 2 mult, 1 div
      float se50 = exp_fast(ryaB50 * INV_R) * UTIMESR / (-ryaB50); //1 exp, 2 mult, 1 div
      float se60 = exp_fast(ryaB60 * INV_R) * UTIMESR / (-ryaB60); //1 exp, 2 mult, 1 div
      float se70 = exp_fast(ryaB70 * INV_R) * UTIMESR / (-ryaB70); //1 exp, 2 mult, 1 div

      float rb00 = se00 * ryaB00; // repulsion border 0 , 1 mult
      float rb10 = se10 * ryaB10; // repulsion border 0 , 1 mult
      float rb20 = se20 * ryaB20; // repulsion border 0 , 1 mult
      float rb30 = se30 * ryaB30; // repulsion border 0 , 1 mult
      float rb40 = se40 * ryaB40; // repulsion border 0 , 1 mult
      float rb50 = se50 * ryaB50; // repulsion border 0 , 1 mult
      float rb60 = se60 * ryaB60; // repulsion border 0 , 1 mult
      float rb70 = se70 * ryaB70; // repulsion border 0 , 1 mult

      float ryaB01 = rya0 - sb; //1 add
      float ryaB11 = rya1 - sb; //1 add
      float ryaB21 = rya2 - sb; //1 add
      float ryaB31 = rya3 - sb; //1 add
      float ryaB41 = rya4 - sb; //1 add
      float ryaB51 = rya5 - sb; //1 add
      float ryaB61 = rya6 - sb; //1 add
      float ryaB71 = rya7 - sb; //1 add

      float se01 = exp_fast(-ryaB01 * INV_R) * UTIMESR / (ryaB01); //1 exp, 2 mult, 1 div
      float se11 = exp_fast(-ryaB11 * INV_R) * UTIMESR / (ryaB11); //1 exp, 2 mult, 1 div
      float se21 = exp_fast(-ryaB21 * INV_R) * UTIMESR / (ryaB21); //1 exp, 2 mult, 1 div
      float se31 = exp_fast(-ryaB31 * INV_R) * UTIMESR / (ryaB31); //1 exp, 2 mult, 1 div
      float se41 = exp_fast(-ryaB41 * INV_R) * UTIMESR / (ryaB41); //1 exp, 2 mult, 1 div
      float se51 = exp_fast(-ryaB51 * INV_R) * UTIMESR / (ryaB51); //1 exp, 2 mult, 1 div
      float se61 = exp_fast(-ryaB61 * INV_R) * UTIMESR / (ryaB61); //1 exp, 2 mult, 1 div
      float se71 = exp_fast(-ryaB71 * INV_R) * UTIMESR / (ryaB71); //1 exp, 2 mult, 1 div

      float rb01 = rb00 + se01 * ryaB01; //1 mult
      float rb11 = rb10 + se11 * ryaB11; //1 mult
      float rb21 = rb20 + se21 * ryaB21; //1 mult
      float rb31 = rb30 + se31 * ryaB31; //1 mult
      float rb41 = rb40 + se41 * ryaB41; //1 mult
      float rb51 = rb50 + se51 * ryaB51; //1 mult
      float rb61 = rb60 + se61 * ryaB61; //1 mult
      float rb71 = rb70 + se71 * ryaB71; //1 mult

      sfy0 += rb01; //1 add
      sfy1 += rb11; //1 add
      sfy2 += rb21; //1 add
      sfy3 += rb31; //1 add
      sfy4 += rb41; //1 add
      sfy5 += rb51; //1 add
      sfy6 += rb61; //1 add
      sfy7 += rb71; //1 add

      // n_borders * (3 addds, 4 mults, 1 div, 1 exp)

      /************************************************/
      // UPDATE PEOPLE REPULSION TERM
      /************************************************/

      float rxa0a1 = rxa0 - rxa1;
      float rya0a1 = rya0 - rya1;
      float rxa0a2 = rxa0 - rxa2;
      float rya0a2 = rya0 - rya2;
      float rxa0a3 = rxa0 - rxa3;
      float rya0a3 = rya0 - rya3;
      float rxa0a4 = rxa0 - rxa4;
      float rya0a4 = rya0 - rya4;
      float rxa0a5 = rxa0 - rxa5;
      float rya0a5 = rya0 - rya5;
      float rxa0a6 = rxa0 - rxa6;
      float rya0a6 = rya0 - rya6;
      float rxa0a7 = rxa0 - rxa7;
      float rya0a7 = rya0 - rya7;

      float rxa1a2 = rxa1 - rxa2;
      float rya1a2 = rya1 - rya2;
      float rxa1a3 = rxa1 - rxa3;
      float rya1a3 = rya1 - rya3;
      float rxa1a4 = rxa1 - rxa4;
      float rya1a4 = rya1 - rya4;
      float rxa1a5 = rxa1 - rxa5;
      float rya1a5 = rya1 - rya5;
      float rxa1a6 = rxa1 - rxa6;
      float rya1a6 = rya1 - rya6;
      float rxa1a7 = rxa1 - rxa7;
      float rya1a7 = rya1 - rya7;

      float rxa2a3 = rxa2 - rxa3;
      float rya2a3 = rya2 - rya3;
      float rxa2a4 = rxa2 - rxa4;
      float rya2a4 = rya2 - rya4;
      float rxa2a5 = rxa2 - rxa5;
      float rya2a5 = rya2 - rya5;
      float rxa2a6 = rxa2 - rxa6;
      float rya2a6 = rya2 - rya6;
      float rxa2a7 = rxa2 - rxa7;
      float rya2a7 = rya2 - rya7;

      float rxa3a4 = rxa3 - rxa4;
      float rya3a4 = rya3 - rya4;
      float rxa3a5 = rxa3 - rxa5;
      float rya3a5 = rya3 - rya5;
      float rxa3a6 = rxa3 - rxa6;
      float rya3a6 = rya3 - rya6;
      float rxa3a7 = rxa3 - rxa7;
      float rya3a7 = rya3 - rya7;

      float rxa4a5 = rxa4 - rxa5;
      float rya4a5 = rya4 - rya5;
      float rxa4a6 = rxa4 - rxa6;
      float rya4a6 = rya4 - rya6;
      float rxa4a7 = rxa4 - rxa7;
      float rya4a7 = rya4 - rya7;

      float rxa5a6 = rxa5 - rxa6;
      float rya5a6 = rya5 - rya6;
      float rxa5a7 = rxa5 - rxa7;
      float rya5a7 = rya5 - rya7;

      float rxa6a7 = rxa6 - rxa7;
      float rya6a7 = rya6 - rya7;

      float ra0a1 = rxa0a1 * rxa0a1 + rya0a1 * rya0a1;
      float ra0a2 = rxa0a2 * rxa0a2 + rya0a2 * rya0a2;
      float ra0a3 = rxa0a3 * rxa0a3 + rya0a3 * rya0a3;
      float ra0a4 = rxa0a4 * rxa0a4 + rya0a4 * rya0a4;
      float ra0a5 = rxa0a5 * rxa0a5 + rya0a5 * rya0a5;
      float ra0a6 = rxa0a6 * rxa0a6 + rya0a6 * rya0a6;
      float ra0a7 = rxa0a7 * rxa0a7 + rya0a7 * rya0a7;

      float ra1a2 = rxa1a2 * rxa1a2 + rya1a2 * rya1a2;
      float ra1a3 = rxa1a3 * rxa1a3 + rya1a3 * rya1a3;
      float ra1a4 = rxa1a4 * rxa1a4 + rya1a4 * rya1a4;
      float ra1a5 = rxa1a5 * rxa1a5 + rya1a5 * rya1a5;
      float ra1a6 = rxa1a6 * rxa1a6 + rya1a6 * rya1a6;
      float ra1a7 = rxa1a7 * rxa1a7 + rya1a7 * rya1a7;

      float ra2a3 = rxa2a3 * rxa2a3 + rya2a3 * rya2a3;
      float ra2a4 = rxa2a4 * rxa2a4 + rya2a4 * rya2a4;
      float ra2a5 = rxa2a5 * rxa2a5 + rya2a5 * rya2a5;
      float ra2a6 = rxa2a6 * rxa2a6 + rya2a6 * rya2a6;
      float ra2a7 = rxa2a7 * rxa2a7 + rya2a7 * rya2a7;

      float ra3a4 = rxa3a4 * rxa3a4 + rya3a4 * rya3a4;
      float ra3a5 = rxa3a5 * rxa3a5 + rya3a5 * rya3a5;
      float ra3a6 = rxa3a6 * rxa3a6 + rya3a6 * rya3a6;
      float ra3a7 = rxa3a7 * rxa3a7 + rya3a7 * rya3a7;

      float ra4a5 = rxa4a5 * rxa4a5 + rya4a5 * rya4a5;
      float ra4a6 = rxa4a6 * rxa4a6 + rya4a6 * rya4a6;
      float ra4a7 = rxa4a7 * rxa4a7 + rya4a7 * rya4a7;

      float ra5a6 = rxa5a6 * rxa5a6 + rya5a6 * rya5a6;
      float ra5a7 = rxa5a7 * rxa5a7 + rya5a7 * rya5a7;

      float ra6a7 = rxa6a7 * rxa6a7 + rya6a7 * rya6a7;

      float ra0a1norm = sqrt(ra0a1);
      float ra0a2norm = sqrt(ra0a2);
      float ra0a3norm = sqrt(ra0a3);
      float ra0a4norm = sqrt(ra0a4);
      float ra0a5norm = sqrt(ra0a5);
      float ra0a6norm = sqrt(ra0a6);
      float ra0a7norm = sqrt(ra0a7);

      float ra1a2norm = sqrt(ra1a2);
      float ra1a3norm = sqrt(ra1a3);
      float ra1a4norm = sqrt(ra1a4);
      float ra1a5norm = sqrt(ra1a5);
      float ra1a6norm = sqrt(ra1a6);
      float ra1a7norm = sqrt(ra1a7);

      float ra2a3norm = sqrt(ra2a3);
      float ra2a4norm = sqrt(ra2a4);
      float ra2a5norm = sqrt(ra2a5);
      float ra2a6norm = sqrt(ra2a6);
      float ra2a7norm = sqrt(ra2a7);

      float ra3a4norm = sqrt(ra3a4);
      float ra3a5norm = sqrt(ra3a5);
      float ra3a6norm = sqrt(ra3a6);
      float ra3a7norm = sqrt(ra3a7);

      float ra4a5norm = sqrt(ra4a5);
      float ra4a6norm = sqrt(ra4a6);
      float ra4a7norm = sqrt(ra4a7);

      float ra5a6norm = sqrt(ra5a6);
      float ra5a7norm = sqrt(ra5a7);

      float ra6a7norm = sqrt(ra6a7);

      //everything twice

      float rxa0a1mex = rxa0a1 - da1 * exa1; //1 add, 1 mult
      float rya0a1mey = rya0a1 - da1 * eya1; //1 add, 1 mult
      float rxa0a2mex = rxa0a2 - da2 * exa2; //1 add, 1 mult
      float rya0a2mey = rya0a2 - da2 * eya2; //1 add, 1 mult
      float rxa0a3mex = rxa0a3 - da3 * exa3; //1 add, 1 mult
      float rya0a3mey = rya0a3 - da3 * eya3; //1 add, 1 mult
      float rxa0a4mex = rxa0a4 - da4 * exa4; //1 add, 1 mult
      float rya0a4mey = rya0a4 - da4 * eya4; //1 add, 1 mult
      float rxa0a5mex = rxa0a5 - da5 * exa5; //1 add, 1 mult
      float rya0a5mey = rya0a5 - da5 * eya5; //1 add, 1 mult
      float rxa0a6mex = rxa0a6 - da6 * exa6; //1 add, 1 mult
      float rya0a6mey = rya0a6 - da6 * eya6; //1 add, 1 mult
      float rxa0a7mex = rxa0a7 - da7 * exa7; //1 add, 1 mult
      float rya0a7mey = rya0a7 - da7 * eya7; //1 add, 1 mult

      float rxa1a0mex = -rxa0a1 - da0 * exa0; //1 add, 1 mult
      float rya1a0mey = -rya0a1 - da0 * eya0; //1 add, 1 mult
      float rxa2a0mex = -rxa0a2 - da0 * exa0; //1 add, 1 mult
      float rya2a0mey = -rya0a2 - da0 * eya0; //1 add, 1 mult
      float rxa3a0mex = -rxa0a3 - da0 * exa0; //1 add, 1 mult
      float rya3a0mey = -rya0a3 - da0 * eya0; //1 add, 1 mult
      float rxa4a0mex = -rxa0a4 - da0 * exa0; //1 add, 1 mult
      float rya4a0mey = -rya0a4 - da0 * eya0; //1 add, 1 mult
      float rxa5a0mex = -rxa0a5 - da0 * exa0; //1 add, 1 mult
      float rya5a0mey = -rya0a5 - da0 * eya0; //1 add, 1 mult
      float rxa6a0mex = -rxa0a6 - da0 * exa0; //1 add, 1 mult
      float rya6a0mey = -rya0a6 - da0 * eya0; //1 add, 1 mult
      float rxa7a0mex = -rxa0a7 - da0 * exa0; //1 add, 1 mult
      float rya7a0mey = -rya0a7 - da0 * eya0; //1 add, 1 mult

      float rxa1a2mex = rxa1a2 - da2 * exa2; //1 add, 1 mult
      float rya1a2mey = rya1a2 - da2 * eya2; //1 add, 1 mult
      float rxa1a3mex = rxa1a3 - da3 * exa3; //1 add, 1 mult
      float rya1a3mey = rya1a3 - da3 * eya3; //1 add, 1 mult
      float rxa1a4mex = rxa1a4 - da4 * exa4; //1 add, 1 mult
      float rya1a4mey = rya1a4 - da4 * eya4; //1 add, 1 mult
      float rxa1a5mex = rxa1a5 - da5 * exa5; //1 add, 1 mult
      float rya1a5mey = rya1a5 - da5 * eya5; //1 add, 1 mult
      float rxa1a6mex = rxa1a6 - da6 * exa6; //1 add, 1 mult
      float rya1a6mey = rya1a6 - da6 * eya6; //1 add, 1 mult
      float rxa1a7mex = rxa1a7 - da7 * exa7; //1 add, 1 mult
      float rya1a7mey = rya1a7 - da7 * eya7; //1 add, 1 mult

      float rxa2a1mex = -rxa1a2 - da1 * exa1; //1 add, 1 mult
      float rya2a1mey = -rya1a2 - da1 * eya1; //1 add, 1 mult
      float rxa3a1mex = -rxa1a3 - da1 * exa1; //1 add, 1 mult
      float rya3a1mey = -rya1a3 - da1 * eya1; //1 add, 1 mult
      float rxa4a1mex = -rxa1a4 - da1 * exa1; //1 add, 1 mult
      float rya4a1mey = -rya1a4 - da1 * eya1; //1 add, 1 mult
      float rxa5a1mex = -rxa1a5 - da1 * exa1; //1 add, 1 mult
      float rya5a1mey = -rya1a5 - da1 * eya1; //1 add, 1 mult
      float rxa6a1mex = -rxa1a6 - da1 * exa1; //1 add, 1 mult
      float rya6a1mey = -rya1a6 - da1 * eya1; //1 add, 1 mult
      float rxa7a1mex = -rxa1a7 - da1 * exa1; //1 add, 1 mult
      float rya7a1mey = -rya1a7 - da1 * eya1; //1 add, 1 mult

      float rxa2a3mex = rxa2a3 - da3 * exa3; //1 add, 1 mult
      float rya2a3mey = rya2a3 - da3 * eya3; //1 add, 1 mult
      float rxa2a4mex = rxa2a4 - da4 * exa4; //1 add, 1 mult
      float rya2a4mey = rya2a4 - da4 * eya4; //1 add, 1 mult
      float rxa2a5mex = rxa2a5 - da5 * exa5; //1 add, 1 mult
      float rya2a5mey = rya2a5 - da5 * eya5; //1 add, 1 mult
      float rxa2a6mex = rxa2a6 - da6 * exa6; //1 add, 1 mult
      float rya2a6mey = rya2a6 - da6 * eya6; //1 add, 1 mult
      float rxa2a7mex = rxa2a7 - da7 * exa7; //1 add, 1 mult
      float rya2a7mey = rya2a7 - da7 * eya7; //1 add, 1 mult

      float rxa3a2mex = -rxa2a3 - da2 * exa2; //1 add, 1 mult
      float rya3a2mey = -rya2a3 - da2 * eya2; //1 add, 1 mult
      float rxa4a2mex = -rxa2a4 - da2 * exa2; //1 add, 1 mult
      float rya4a2mey = -rya2a4 - da2 * eya2; //1 add, 1 mult
      float rxa5a2mex = -rxa2a5 - da2 * exa2; //1 add, 1 mult
      float rya5a2mey = -rya2a5 - da2 * eya2; //1 add, 1 mult
      float rxa6a2mex = -rxa2a6 - da2 * exa2; //1 add, 1 mult
      float rya6a2mey = -rya2a6 - da2 * eya2; //1 add, 1 mult
      float rxa7a2mex = -rxa2a7 - da2 * exa2; //1 add, 1 mult
      float rya7a2mey = -rya2a7 - da2 * eya2; //1 add, 1 mult

      float rxa3a4mex = rxa3a4 - da4 * exa4; //1 add, 1 mult
      float rya3a4mey = rya3a4 - da4 * eya4; //1 add, 1 mult
      float rxa3a5mex = rxa3a5 - da5 * exa5; //1 add, 1 mult
      float rya3a5mey = rya3a5 - da5 * eya5; //1 add, 1 mult
      float rxa3a6mex = rxa3a6 - da6 * exa6; //1 add, 1 mult
      float rya3a6mey = rya3a6 - da6 * eya6; //1 add, 1 mult
      float rxa3a7mex = rxa3a7 - da7 * exa7; //1 add, 1 mult
      float rya3a7mey = rya3a7 - da7 * eya7; //1 add, 1 mult

      float rxa4a3mex = -rxa3a4 - da3 * exa3; //1 add, 1 mult
      float rya4a3mey = -rya3a4 - da3 * eya3; //1 add, 1 mult
      float rxa5a3mex = -rxa3a5 - da3 * exa3; //1 add, 1 mult
      float rya5a3mey = -rya3a5 - da3 * eya3; //1 add, 1 mult
      float rxa6a3mex = -rxa3a6 - da3 * exa3; //1 add, 1 mult
      float rya6a3mey = -rya3a6 - da3 * eya3; //1 add, 1 mult
      float rxa7a3mex = -rxa3a7 - da3 * exa3; //1 add, 1 mult
      float rya7a3mey = -rya3a7 - da3 * eya3; //1 add, 1 mult

      float rxa4a5mex = rxa4a5 - da5 * exa5; //1 add, 1 mult
      float rya4a5mey = rya4a5 - da5 * eya5; //1 add, 1 mult
      float rxa4a6mex = rxa4a6 - da6 * exa6; //1 add, 1 mult
      float rya4a6mey = rya4a6 - da6 * eya6; //1 add, 1 mult
      float rxa4a7mex = rxa4a7 - da7 * exa7; //1 add, 1 mult
      float rya4a7mey = rya4a7 - da7 * eya7; //1 add, 1 mult

      float rxa5a4mex = -rxa4a5 - da4 * exa4; //1 add, 1 mult
      float rya5a4mey = -rya4a5 - da4 * eya4; //1 add, 1 mult
      float rxa6a4mex = -rxa4a6 - da4 * exa4; //1 add, 1 mult
      float rya6a4mey = -rya4a6 - da4 * eya4; //1 add, 1 mult
      float rxa7a4mex = -rxa4a7 - da4 * exa4; //1 add, 1 mult
      float rya7a4mey = -rya4a7 - da4 * eya4; //1 add, 1 mult

      float rxa5a6mex = rxa5a6 - da6 * exa6; //1 add, 1 mult
      float rya5a6mey = rya5a6 - da6 * eya6; //1 add, 1 mult
      float rxa5a7mex = rxa5a7 - da7 * exa7; //1 add, 1 mult
      float rya5a7mey = rya5a7 - da7 * eya7; //1 add, 1 mult

      float rxa6a5mex = -rxa5a6 - da5 * exa5; //1 add, 1 mult
      float rya6a5mey = -rya5a6 - da5 * eya5; //1 add, 1 mult
      float rxa7a5mex = -rxa5a7 - da5 * exa5; //1 add, 1 mult
      float rya7a5mey = -rya5a7 - da5 * eya5; //1 add, 1 mult

      float rxa6a7mex = rxa6a7 - da7 * exa7; //1 add, 1 mult
      float rya6a7mey = rya6a7 - da7 * eya7; //1 add, 1 mult

      float rxa7a6mex = -rxa6a7 - da6 * exa6; //1 add, 1 mult
      float rya7a6mey = -rya6a7 - da6 * eya6; //1 add, 1 mult

      float ra0a1me = rxa0a1mex * rxa0a1mex + rya0a1mey * rya0a1mey;
      float ra0a2me = rxa0a2mex * rxa0a2mex + rya0a2mey * rya0a2mey;
      float ra0a3me = rxa0a3mex * rxa0a3mex + rya0a3mey * rya0a3mey;
      float ra0a4me = rxa0a4mex * rxa0a4mex + rya0a4mey * rya0a4mey;
      float ra0a5me = rxa0a5mex * rxa0a5mex + rya0a5mey * rya0a5mey;
      float ra0a6me = rxa0a6mex * rxa0a6mex + rya0a6mey * rya0a6mey;
      float ra0a7me = rxa0a7mex * rxa0a7mex + rya0a7mey * rya0a7mey;

      float ra1a0me = rxa1a0mex * rxa1a0mex + rya1a0mey * rya1a0mey;
      float ra2a0me = rxa2a0mex * rxa2a0mex + rya2a0mey * rya2a0mey;
      float ra3a0me = rxa3a0mex * rxa3a0mex + rya3a0mey * rya3a0mey;
      float ra4a0me = rxa4a0mex * rxa4a0mex + rya4a0mey * rya4a0mey;
      float ra5a0me = rxa5a0mex * rxa5a0mex + rya5a0mey * rya5a0mey;
      float ra6a0me = rxa6a0mex * rxa6a0mex + rya6a0mey * rya6a0mey;
      float ra7a0me = rxa7a0mex * rxa7a0mex + rya7a0mey * rya7a0mey;

      float ra1a2me = rxa1a2mex * rxa1a2mex + rya1a2mey * rya1a2mey;
      float ra1a3me = rxa1a3mex * rxa1a3mex + rya1a3mey * rya1a3mey;
      float ra1a4me = rxa1a4mex * rxa1a4mex + rya1a4mey * rya1a4mey;
      float ra1a5me = rxa1a5mex * rxa1a5mex + rya1a5mey * rya1a5mey;
      float ra1a6me = rxa1a6mex * rxa1a6mex + rya1a6mey * rya1a6mey;
      float ra1a7me = rxa1a7mex * rxa1a7mex + rya1a7mey * rya1a7mey;

      float ra2a1me = rxa2a1mex * rxa2a1mex + rya2a1mey * rya2a1mey;
      float ra3a1me = rxa3a1mex * rxa3a1mex + rya3a1mey * rya3a1mey;
      float ra4a1me = rxa4a1mex * rxa4a1mex + rya4a1mey * rya4a1mey;
      float ra5a1me = rxa5a1mex * rxa5a1mex + rya5a1mey * rya5a1mey;
      float ra6a1me = rxa6a1mex * rxa6a1mex + rya6a1mey * rya6a1mey;
      float ra7a1me = rxa7a1mex * rxa7a1mex + rya7a1mey * rya7a1mey;

      float ra2a3me = rxa2a3mex * rxa2a3mex + rya2a3mey * rya2a3mey;
      float ra2a4me = rxa2a4mex * rxa2a4mex + rya2a4mey * rya2a4mey;
      float ra2a5me = rxa2a5mex * rxa2a5mex + rya2a5mey * rya2a5mey;
      float ra2a6me = rxa2a6mex * rxa2a6mex + rya2a6mey * rya2a6mey;
      float ra2a7me = rxa2a7mex * rxa2a7mex + rya2a7mey * rya2a7mey;

      float ra3a2me = rxa3a2mex * rxa3a2mex + rya3a2mey * rya3a2mey;
      float ra4a2me = rxa4a2mex * rxa4a2mex + rya4a2mey * rya4a2mey;
      float ra5a2me = rxa5a2mex * rxa5a2mex + rya5a2mey * rya5a2mey;
      float ra6a2me = rxa6a2mex * rxa6a2mex + rya6a2mey * rya6a2mey;
      float ra7a2me = rxa7a2mex * rxa7a2mex + rya7a2mey * rya7a2mey;

      float ra3a4me = rxa3a4mex * rxa3a4mex + rya3a4mey * rya3a4mey;
      float ra3a5me = rxa3a5mex * rxa3a5mex + rya3a5mey * rya3a5mey;
      float ra3a6me = rxa3a6mex * rxa3a6mex + rya3a6mey * rya3a6mey;
      float ra3a7me = rxa3a7mex * rxa3a7mex + rya3a7mey * rya3a7mey;

      float ra4a3me = rxa4a3mex * rxa4a3mex + rya4a3mey * rya4a3mey;
      float ra5a3me = rxa5a3mex * rxa5a3mex + rya5a3mey * rya5a3mey;
      float ra6a3me = rxa6a3mex * rxa6a3mex + rya6a3mey * rya6a3mey;
      float ra7a3me = rxa7a3mex * rxa7a3mex + rya7a3mey * rya7a3mey;

      float ra4a5me = rxa4a5mex * rxa4a5mex + rya4a5mey * rya4a5mey;
      float ra4a6me = rxa4a6mex * rxa4a6mex + rya4a6mey * rya4a6mey;
      float ra4a7me = rxa4a7mex * rxa4a7mex + rya4a7mey * rya4a7mey;

      float ra5a4me = rxa5a4mex * rxa5a4mex + rya5a4mey * rya5a4mey;
      float ra6a4me = rxa6a4mex * rxa6a4mex + rya6a4mey * rya6a4mey;
      float ra7a4me = rxa7a4mex * rxa7a4mex + rya7a4mey * rya7a4mey;

      float ra5a6me = rxa5a6mex * rxa5a6mex + rya5a6mey * rya5a6mey;
      float ra5a7me = rxa5a7mex * rxa5a7mex + rya5a7mey * rya5a7mey;

      float ra6a5me = rxa6a5mex * rxa6a5mex + rya6a5mey * rya6a5mey;
      float ra7a5me = rxa7a5mex * rxa7a5mex + rya7a5mey * rya7a5mey;

      float ra6a7me = rxa6a7mex * rxa6a7mex + rya6a7mey * rya6a7mey;

      float ra7a6me = rxa7a6mex * rxa7a6mex + rya7a6mey * rya7a6mey;

      float ra0a1menorm = sqrt(ra0a1me);
      float ra0a2menorm = sqrt(ra0a2me);
      float ra0a3menorm = sqrt(ra0a3me);
      float ra0a4menorm = sqrt(ra0a4me);
      float ra0a5menorm = sqrt(ra0a5me);
      float ra0a6menorm = sqrt(ra0a6me);
      float ra0a7menorm = sqrt(ra0a7me);

      float ra1a0menorm = sqrt(ra1a0me);
      float ra2a0menorm = sqrt(ra2a0me);
      float ra3a0menorm = sqrt(ra3a0me);
      float ra4a0menorm = sqrt(ra4a0me);
      float ra5a0menorm = sqrt(ra5a0me);
      float ra6a0menorm = sqrt(ra6a0me);
      float ra7a0menorm = sqrt(ra7a0me);

      float ra1a2menorm = sqrt(ra1a2me);
      float ra1a3menorm = sqrt(ra1a3me);
      float ra1a4menorm = sqrt(ra1a4me);
      float ra1a5menorm = sqrt(ra1a5me);
      float ra1a6menorm = sqrt(ra1a6me);
      float ra1a7menorm = sqrt(ra1a7me);

      float ra2a1menorm = sqrt(ra2a1me);
      float ra3a1menorm = sqrt(ra3a1me);
      float ra4a1menorm = sqrt(ra4a1me);
      float ra5a1menorm = sqrt(ra5a1me);
      float ra6a1menorm = sqrt(ra6a1me);
      float ra7a1menorm = sqrt(ra7a1me);

      float ra2a3menorm = sqrt(ra2a3me);
      float ra2a4menorm = sqrt(ra2a4me);
      float ra2a5menorm = sqrt(ra2a5me);
      float ra2a6menorm = sqrt(ra2a6me);
      float ra2a7menorm = sqrt(ra2a7me);

      float ra3a2menorm = sqrt(ra3a2me);
      float ra4a2menorm = sqrt(ra4a2me);
      float ra5a2menorm = sqrt(ra5a2me);
      float ra6a2menorm = sqrt(ra6a2me);
      float ra7a2menorm = sqrt(ra7a2me);

      float ra3a4menorm = sqrt(ra3a4me);
      float ra3a5menorm = sqrt(ra3a5me);
      float ra3a6menorm = sqrt(ra3a6me);
      float ra3a7menorm = sqrt(ra3a7me);

      float ra4a3menorm = sqrt(ra4a3me);
      float ra5a3menorm = sqrt(ra5a3me);
      float ra6a3menorm = sqrt(ra6a3me);
      float ra7a3menorm = sqrt(ra7a3me);

      float ra4a5menorm = sqrt(ra4a5me);
      float ra4a6menorm = sqrt(ra4a6me);
      float ra4a7menorm = sqrt(ra4a7me);

      float ra5a4menorm = sqrt(ra5a4me);
      float ra6a4menorm = sqrt(ra6a4me);
      float ra7a4menorm = sqrt(ra7a4me);

      float ra5a6menorm = sqrt(ra5a6me);
      float ra5a7menorm = sqrt(ra5a7me);

      float ra6a5menorm = sqrt(ra6a5me);
      float ra7a5menorm = sqrt(ra7a5me);

      float ra6a7menorm = sqrt(ra6a7me);

      float ra7a6menorm = sqrt(ra7a6me);

      float normsuma0a1 = ra0a1norm + ra0a1menorm;
      float normsuma0a2 = ra0a2norm + ra0a2menorm;
      float normsuma0a3 = ra0a3norm + ra0a3menorm;
      float normsuma0a4 = ra0a4norm + ra0a4menorm;
      float normsuma0a5 = ra0a5norm + ra0a5menorm;
      float normsuma0a6 = ra0a6norm + ra0a6menorm;
      float normsuma0a7 = ra0a7norm + ra0a7menorm;

      float normsuma1a0 = ra0a1norm + ra1a0menorm;
      float normsuma2a0 = ra0a2norm + ra2a0menorm;
      float normsuma3a0 = ra0a3norm + ra3a0menorm;
      float normsuma4a0 = ra0a4norm + ra4a0menorm;
      float normsuma5a0 = ra0a5norm + ra5a0menorm;
      float normsuma6a0 = ra0a6norm + ra6a0menorm;
      float normsuma7a0 = ra0a7norm + ra7a0menorm;

      float normsuma1a2 = ra1a2norm + ra1a2menorm;
      float normsuma1a3 = ra1a3norm + ra1a3menorm;
      float normsuma1a4 = ra1a4norm + ra1a4menorm;
      float normsuma1a5 = ra1a5norm + ra1a5menorm;
      float normsuma1a6 = ra1a6norm + ra1a6menorm;
      float normsuma1a7 = ra1a7norm + ra1a7menorm;

      float normsuma2a1 = ra1a2norm + ra2a1menorm;
      float normsuma3a1 = ra1a3norm + ra3a1menorm;
      float normsuma4a1 = ra1a4norm + ra4a1menorm;
      float normsuma5a1 = ra1a5norm + ra5a1menorm;
      float normsuma6a1 = ra1a6norm + ra6a1menorm;
      float normsuma7a1 = ra1a7norm + ra7a1menorm;

      float normsuma2a3 = ra2a3norm + ra2a3menorm;
      float normsuma2a4 = ra2a4norm + ra2a4menorm;
      float normsuma2a5 = ra2a5norm + ra2a5menorm;
      float normsuma2a6 = ra2a6norm + ra2a6menorm;
      float normsuma2a7 = ra2a7norm + ra2a7menorm;

      float normsuma3a2 = ra2a3norm + ra3a2menorm;
      float normsuma4a2 = ra2a4norm + ra4a2menorm;
      float normsuma5a2 = ra2a5norm + ra5a2menorm;
      float normsuma6a2 = ra2a6norm + ra6a2menorm;
      float normsuma7a2 = ra2a7norm + ra7a2menorm;

      float normsuma3a4 = ra3a4norm + ra3a4menorm;
      float normsuma3a5 = ra3a5norm + ra3a5menorm;
      float normsuma3a6 = ra3a6norm + ra3a6menorm;
      float normsuma3a7 = ra3a7norm + ra3a7menorm;

      float normsuma4a3 = ra3a4norm + ra4a3menorm;
      float normsuma5a3 = ra3a5norm + ra5a3menorm;
      float normsuma6a3 = ra3a6norm + ra6a3menorm;
      float normsuma7a3 = ra3a7norm + ra7a3menorm;

      float normsuma4a5 = ra4a5norm + ra4a5menorm;
      float normsuma4a6 = ra4a6norm + ra4a6menorm;
      float normsuma4a7 = ra4a7norm + ra4a7menorm;

      float normsuma5a4 = ra4a5norm + ra5a4menorm;
      float normsuma6a4 = ra4a6norm + ra6a4menorm;
      float normsuma7a4 = ra4a7norm + ra7a4menorm;

      float normsuma5a6 = ra5a6norm + ra5a6menorm;
      float normsuma5a7 = ra5a7norm + ra5a7menorm;

      float normsuma6a5 = ra5a6norm + ra6a5menorm;
      float normsuma7a5 = ra5a7norm + ra7a5menorm;

      float normsuma6a7 = ra6a7norm + ra6a7menorm;

      float normsuma7a6 = ra6a7norm + ra7a6menorm;

      float repxa0a1 = rxa0a1 / ra0a1norm + rxa0a1mex / ra0a1menorm;
      float repya0a1 = rya0a1 / ra0a1norm + rya0a1mey / ra0a1menorm;
      float repxa0a2 = rxa0a2 / ra0a2norm + rxa0a2mex / ra0a2menorm;
      float repya0a2 = rya0a2 / ra0a2norm + rya0a2mey / ra0a2menorm;
      float repxa0a3 = rxa0a3 / ra0a3norm + rxa0a3mex / ra0a3menorm;
      float repya0a3 = rya0a3 / ra0a3norm + rya0a3mey / ra0a3menorm;
      float repxa0a4 = rxa0a4 / ra0a4norm + rxa0a4mex / ra0a4menorm;
      float repya0a4 = rya0a4 / ra0a4norm + rya0a4mey / ra0a4menorm;
      float repxa0a5 = rxa0a5 / ra0a5norm + rxa0a5mex / ra0a5menorm;
      float repya0a5 = rya0a5 / ra0a5norm + rya0a5mey / ra0a5menorm;
      float repxa0a6 = rxa0a6 / ra0a6norm + rxa0a6mex / ra0a6menorm;
      float repya0a6 = rya0a6 / ra0a6norm + rya0a6mey / ra0a6menorm;
      float repxa0a7 = rxa0a7 / ra0a7norm + rxa0a7mex / ra0a7menorm;
      float repya0a7 = rya0a7 / ra0a7norm + rya0a7mey / ra0a7menorm;

      float repxa1a0 = -rxa0a1 / ra0a1norm + rxa1a0mex / ra1a0menorm;
      float repya1a0 = -rya0a1 / ra0a1norm + rya1a0mey / ra1a0menorm;
      float repxa2a0 = -rxa0a2 / ra0a2norm + rxa2a0mex / ra2a0menorm;
      float repya2a0 = -rya0a2 / ra0a2norm + rya2a0mey / ra2a0menorm;
      float repxa3a0 = -rxa0a3 / ra0a3norm + rxa3a0mex / ra3a0menorm;
      float repya3a0 = -rya0a3 / ra0a3norm + rya3a0mey / ra3a0menorm;
      float repxa4a0 = -rxa0a4 / ra0a4norm + rxa4a0mex / ra4a0menorm;
      float repya4a0 = -rya0a4 / ra0a4norm + rya4a0mey / ra4a0menorm;
      float repxa5a0 = -rxa0a5 / ra0a5norm + rxa5a0mex / ra5a0menorm;
      float repya5a0 = -rya0a5 / ra0a5norm + rya5a0mey / ra5a0menorm;
      float repxa6a0 = -rxa0a6 / ra0a6norm + rxa6a0mex / ra6a0menorm;
      float repya6a0 = -rya0a6 / ra0a6norm + rya6a0mey / ra6a0menorm;
      float repxa7a0 = -rxa0a7 / ra0a7norm + rxa7a0mex / ra7a0menorm;
      float repya7a0 = -rya0a7 / ra0a7norm + rya7a0mey / ra7a0menorm;

      float repxa1a2 = rxa1a2 / ra1a2norm + rxa1a2mex / ra1a2menorm;
      float repya1a2 = rya1a2 / ra1a2norm + rya1a2mey / ra1a2menorm;
      float repxa1a3 = rxa1a3 / ra1a3norm + rxa1a3mex / ra1a3menorm;
      float repya1a3 = rya1a3 / ra1a3norm + rya1a3mey / ra1a3menorm;
      float repxa1a4 = rxa1a4 / ra1a4norm + rxa1a4mex / ra1a4menorm;
      float repya1a4 = rya1a4 / ra1a4norm + rya1a4mey / ra1a4menorm;
      float repxa1a5 = rxa1a5 / ra1a5norm + rxa1a5mex / ra1a5menorm;
      float repya1a5 = rya1a5 / ra1a5norm + rya1a5mey / ra1a5menorm;
      float repxa1a6 = rxa1a6 / ra1a6norm + rxa1a6mex / ra1a6menorm;
      float repya1a6 = rya1a6 / ra1a6norm + rya1a6mey / ra1a6menorm;
      float repxa1a7 = rxa1a7 / ra1a7norm + rxa1a7mex / ra1a7menorm;
      float repya1a7 = rya1a7 / ra1a7norm + rya1a7mey / ra1a7menorm;

      float repxa2a1 = -rxa1a2 / ra1a2norm + rxa2a1mex / ra2a1menorm;
      float repya2a1 = -rya1a2 / ra1a2norm + rya2a1mey / ra2a1menorm;
      float repxa3a1 = -rxa1a3 / ra1a3norm + rxa3a1mex / ra3a1menorm;
      float repya3a1 = -rya1a3 / ra1a3norm + rya3a1mey / ra3a1menorm;
      float repxa4a1 = -rxa1a4 / ra1a4norm + rxa4a1mex / ra4a1menorm;
      float repya4a1 = -rya1a4 / ra1a4norm + rya4a1mey / ra4a1menorm;
      float repxa5a1 = -rxa1a5 / ra1a5norm + rxa5a1mex / ra5a1menorm;
      float repya5a1 = -rya1a5 / ra1a5norm + rya5a1mey / ra5a1menorm;
      float repxa6a1 = -rxa1a6 / ra1a6norm + rxa6a1mex / ra6a1menorm;
      float repya6a1 = -rya1a6 / ra1a6norm + rya6a1mey / ra6a1menorm;
      float repxa7a1 = -rxa1a7 / ra1a7norm + rxa7a1mex / ra7a1menorm;
      float repya7a1 = -rya1a7 / ra1a7norm + rya7a1mey / ra7a1menorm;

      float repxa2a3 = rxa2a3 / ra2a3norm + rxa2a3mex / ra2a3menorm;
      float repya2a3 = rya2a3 / ra2a3norm + rya2a3mey / ra2a3menorm;
      float repxa2a4 = rxa2a4 / ra2a4norm + rxa2a4mex / ra2a4menorm;
      float repya2a4 = rya2a4 / ra2a4norm + rya2a4mey / ra2a4menorm;
      float repxa2a5 = rxa2a5 / ra2a5norm + rxa2a5mex / ra2a5menorm;
      float repya2a5 = rya2a5 / ra2a5norm + rya2a5mey / ra2a5menorm;
      float repxa2a6 = rxa2a6 / ra2a6norm + rxa2a6mex / ra2a6menorm;
      float repya2a6 = rya2a6 / ra2a6norm + rya2a6mey / ra2a6menorm;
      float repxa2a7 = rxa2a7 / ra2a7norm + rxa2a7mex / ra2a7menorm;
      float repya2a7 = rya2a7 / ra2a7norm + rya2a7mey / ra2a7menorm;

      float repxa3a2 = -rxa2a3 / ra2a3norm + rxa3a2mex / ra3a2menorm;
      float repya3a2 = -rya2a3 / ra2a3norm + rya3a2mey / ra3a2menorm;
      float repxa4a2 = -rxa2a4 / ra2a4norm + rxa4a2mex / ra4a2menorm;
      float repya4a2 = -rya2a4 / ra2a4norm + rya4a2mey / ra4a2menorm;
      float repxa5a2 = -rxa2a5 / ra2a5norm + rxa5a2mex / ra5a2menorm;
      float repya5a2 = -rya2a5 / ra2a5norm + rya5a2mey / ra5a2menorm;
      float repxa6a2 = -rxa2a6 / ra2a6norm + rxa6a2mex / ra6a2menorm;
      float repya6a2 = -rya2a6 / ra2a6norm + rya6a2mey / ra6a2menorm;
      float repxa7a2 = -rxa2a7 / ra2a7norm + rxa7a2mex / ra7a2menorm;
      float repya7a2 = -rya2a7 / ra2a7norm + rya7a2mey / ra7a2menorm;

      float repxa3a4 = rxa3a4 / ra3a4norm + rxa3a4mex / ra3a4menorm;
      float repya3a4 = rya3a4 / ra3a4norm + rya3a4mey / ra3a4menorm;
      float repxa3a5 = rxa3a5 / ra3a5norm + rxa3a5mex / ra3a5menorm;
      float repya3a5 = rya3a5 / ra3a5norm + rya3a5mey / ra3a5menorm;
      float repxa3a6 = rxa3a6 / ra3a6norm + rxa3a6mex / ra3a6menorm;
      float repya3a6 = rya3a6 / ra3a6norm + rya3a6mey / ra3a6menorm;
      float repxa3a7 = rxa3a7 / ra3a7norm + rxa3a7mex / ra3a7menorm;
      float repya3a7 = rya3a7 / ra3a7norm + rya3a7mey / ra3a7menorm;

      float repxa4a3 = -rxa3a4 / ra3a4norm + rxa4a3mex / ra4a3menorm;
      float repya4a3 = -rya3a4 / ra3a4norm + rya4a3mey / ra4a3menorm;
      float repxa5a3 = -rxa3a5 / ra3a5norm + rxa5a3mex / ra5a3menorm;
      float repya5a3 = -rya3a5 / ra3a5norm + rya5a3mey / ra5a3menorm;
      float repxa6a3 = -rxa3a6 / ra3a6norm + rxa6a3mex / ra6a3menorm;
      float repya6a3 = -rya3a6 / ra3a6norm + rya6a3mey / ra6a3menorm;
      float repxa7a3 = -rxa3a7 / ra3a7norm + rxa7a3mex / ra7a3menorm;
      float repya7a3 = -rya3a7 / ra3a7norm + rya7a3mey / ra7a3menorm;

      float repxa4a5 = rxa4a5 / ra4a5norm + rxa4a5mex / ra4a5menorm;
      float repya4a5 = rya4a5 / ra4a5norm + rya4a5mey / ra4a5menorm;
      float repxa4a6 = rxa4a6 / ra4a6norm + rxa4a6mex / ra4a6menorm;
      float repya4a6 = rya4a6 / ra4a6norm + rya4a6mey / ra4a6menorm;
      float repxa4a7 = rxa4a7 / ra4a7norm + rxa4a7mex / ra4a7menorm;
      float repya4a7 = rya4a7 / ra4a7norm + rya4a7mey / ra4a7menorm;

      float repxa5a4 = -rxa4a5 / ra4a5norm + rxa5a4mex / ra5a4menorm;
      float repya5a4 = -rya4a5 / ra4a5norm + rya5a4mey / ra5a4menorm;
      float repxa6a4 = -rxa4a6 / ra4a6norm + rxa6a4mex / ra6a4menorm;
      float repya6a4 = -rya4a6 / ra4a6norm + rya6a4mey / ra6a4menorm;
      float repxa7a4 = -rxa4a7 / ra4a7norm + rxa7a4mex / ra7a4menorm;
      float repya7a4 = -rya4a7 / ra4a7norm + rya7a4mey / ra7a4menorm;

      float repxa5a6 = rxa5a6 / ra5a6norm + rxa5a6mex / ra5a6menorm;
      float repya5a6 = rya5a6 / ra5a6norm + rya5a6mey / ra5a6menorm;
      float repxa5a7 = rxa5a7 / ra5a7norm + rxa5a7mex / ra5a7menorm;
      float repya5a7 = rya5a7 / ra5a7norm + rya5a7mey / ra5a7menorm;

      float repxa6a5 = -rxa5a6 / ra5a6norm + rxa6a5mex / ra6a5menorm;
      float repya6a5 = -rya5a6 / ra5a6norm + rya6a5mey / ra6a5menorm;
      float repxa7a5 = -rxa5a7 / ra5a7norm + rxa7a5mex / ra7a5menorm;
      float repya7a5 = -rya5a7 / ra5a7norm + rya7a5mey / ra7a5menorm;

      float repxa6a7 = rxa6a7 / ra6a7norm + rxa6a7mex / ra6a7menorm;
      float repya6a7 = rya6a7 / ra6a7norm + rya6a7mey / ra6a7menorm;

      float repxa7a6 = -rxa6a7 / ra6a7norm + rxa7a6mex / ra7a6menorm;
      float repya7a6 = -rya6a7 / ra6a7norm + rya7a6mey / ra7a6menorm;

      float ba0a1 = normsuma0a1 * normsuma0a1 - da1 * da1;
      float ba0a2 = normsuma0a2 * normsuma0a2 - da2 * da2;
      float ba0a3 = normsuma0a3 * normsuma0a3 - da3 * da3;
      float ba0a4 = normsuma0a4 * normsuma0a4 - da4 * da4;
      float ba0a5 = normsuma0a5 * normsuma0a5 - da5 * da5;
      float ba0a6 = normsuma0a6 * normsuma0a6 - da6 * da6;
      float ba0a7 = normsuma0a7 * normsuma0a7 - da7 * da7;

      float ba1a0 = normsuma1a0 * normsuma1a0 - da0 * da0;
      float ba2a0 = normsuma2a0 * normsuma2a0 - da0 * da0;
      float ba3a0 = normsuma3a0 * normsuma3a0 - da0 * da0;
      float ba4a0 = normsuma4a0 * normsuma4a0 - da0 * da0;
      float ba5a0 = normsuma5a0 * normsuma5a0 - da0 * da0;
      float ba6a0 = normsuma6a0 * normsuma6a0 - da0 * da0;
      float ba7a0 = normsuma7a0 * normsuma7a0 - da0 * da0;

      float ba1a2 = normsuma1a2 * normsuma1a2 - da2 * da2;
      float ba1a3 = normsuma1a3 * normsuma1a3 - da3 * da3;
      float ba1a4 = normsuma1a4 * normsuma1a4 - da4 * da4;
      float ba1a5 = normsuma1a5 * normsuma1a5 - da5 * da5;
      float ba1a6 = normsuma1a6 * normsuma1a6 - da6 * da6;
      float ba1a7 = normsuma1a7 * normsuma1a7 - da7 * da7;

      float ba2a1 = normsuma2a1 * normsuma2a1 - da1 * da1;
      float ba3a1 = normsuma3a1 * normsuma3a1 - da1 * da1;
      float ba4a1 = normsuma4a1 * normsuma4a1 - da1 * da1;
      float ba5a1 = normsuma5a1 * normsuma5a1 - da1 * da1;
      float ba6a1 = normsuma6a1 * normsuma6a1 - da1 * da1;
      float ba7a1 = normsuma7a1 * normsuma7a1 - da1 * da1;

      float ba2a3 = normsuma2a3 * normsuma2a3 - da3 * da3;
      float ba2a4 = normsuma2a4 * normsuma2a4 - da4 * da4;
      float ba2a5 = normsuma2a5 * normsuma2a5 - da5 * da5;
      float ba2a6 = normsuma2a6 * normsuma2a6 - da6 * da6;
      float ba2a7 = normsuma2a7 * normsuma2a7 - da7 * da7;

      float ba3a2 = normsuma3a2 * normsuma3a2 - da2 * da2;
      float ba4a2 = normsuma4a2 * normsuma4a2 - da2 * da2;
      float ba5a2 = normsuma5a2 * normsuma5a2 - da2 * da2;
      float ba6a2 = normsuma6a2 * normsuma6a2 - da2 * da2;
      float ba7a2 = normsuma7a2 * normsuma7a2 - da2 * da2;

      float ba3a4 = normsuma3a4 * normsuma3a4 - da4 * da4;
      float ba3a5 = normsuma3a5 * normsuma3a5 - da5 * da5;
      float ba3a6 = normsuma3a6 * normsuma3a6 - da6 * da6;
      float ba3a7 = normsuma3a7 * normsuma3a7 - da7 * da7;

      float ba4a3 = normsuma4a3 * normsuma4a3 - da3 * da3;
      float ba5a3 = normsuma5a3 * normsuma5a3 - da3 * da3;
      float ba6a3 = normsuma6a3 * normsuma6a3 - da3 * da3;
      float ba7a3 = normsuma7a3 * normsuma7a3 - da3 * da3;

      float ba4a5 = normsuma4a5 * normsuma4a5 - da5 * da5;
      float ba4a6 = normsuma4a6 * normsuma4a6 - da6 * da6;
      float ba4a7 = normsuma4a7 * normsuma4a7 - da7 * da7;

      float ba5a4 = normsuma5a4 * normsuma5a4 - da4 * da4;
      float ba6a4 = normsuma6a4 * normsuma6a4 - da4 * da4;
      float ba7a4 = normsuma7a4 * normsuma7a4 - da4 * da4;

      float ba5a6 = normsuma5a6 * normsuma5a6 - da6 * da6;
      float ba5a7 = normsuma5a7 * normsuma5a7 - da7 * da7;

      float ba6a5 = normsuma6a5 * normsuma6a5 - da5 * da5;
      float ba7a5 = normsuma7a5 * normsuma7a5 - da5 * da5;

      float ba6a7 = normsuma6a7 * normsuma6a7 - da7 * da7;

      float ba7a6 = normsuma7a6 * normsuma7a6 - da6 * da6;

      float ba0a1norm = sqrt(ba0a1) * 0.5;
      float ba0a2norm = sqrt(ba0a2) * 0.5;
      float ba0a3norm = sqrt(ba0a3) * 0.5;
      float ba0a4norm = sqrt(ba0a4) * 0.5;
      float ba0a5norm = sqrt(ba0a5) * 0.5;
      float ba0a6norm = sqrt(ba0a6) * 0.5;
      float ba0a7norm = sqrt(ba0a7) * 0.5;

      float ba1a0norm = sqrt(ba1a0) * 0.5;
      float ba2a0norm = sqrt(ba2a0) * 0.5;
      float ba3a0norm = sqrt(ba3a0) * 0.5;
      float ba4a0norm = sqrt(ba4a0) * 0.5;
      float ba5a0norm = sqrt(ba5a0) * 0.5;
      float ba6a0norm = sqrt(ba6a0) * 0.5;
      float ba7a0norm = sqrt(ba7a0) * 0.5;

      float ba1a2norm = sqrt(ba1a2) * 0.5;
      float ba1a3norm = sqrt(ba1a3) * 0.5;
      float ba1a4norm = sqrt(ba1a4) * 0.5;
      float ba1a5norm = sqrt(ba1a5) * 0.5;
      float ba1a6norm = sqrt(ba1a6) * 0.5;
      float ba1a7norm = sqrt(ba1a7) * 0.5;

      float ba2a1norm = sqrt(ba2a1) * 0.5;
      float ba3a1norm = sqrt(ba3a1) * 0.5;
      float ba4a1norm = sqrt(ba4a1) * 0.5;
      float ba5a1norm = sqrt(ba5a1) * 0.5;
      float ba6a1norm = sqrt(ba6a1) * 0.5;
      float ba7a1norm = sqrt(ba7a1) * 0.5;

      float ba2a3norm = sqrt(ba2a3) * 0.5;
      float ba2a4norm = sqrt(ba2a4) * 0.5;
      float ba2a5norm = sqrt(ba2a5) * 0.5;
      float ba2a6norm = sqrt(ba2a6) * 0.5;
      float ba2a7norm = sqrt(ba2a7) * 0.5;

      float ba3a2norm = sqrt(ba3a2) * 0.5;
      float ba4a2norm = sqrt(ba4a2) * 0.5;
      float ba5a2norm = sqrt(ba5a2) * 0.5;
      float ba6a2norm = sqrt(ba6a2) * 0.5;
      float ba7a2norm = sqrt(ba7a2) * 0.5;

      float ba3a4norm = sqrt(ba3a4) * 0.5;
      float ba3a5norm = sqrt(ba3a5) * 0.5;
      float ba3a6norm = sqrt(ba3a6) * 0.5;
      float ba3a7norm = sqrt(ba3a7) * 0.5;

      float ba4a3norm = sqrt(ba4a3) * 0.5;
      float ba5a3norm = sqrt(ba5a3) * 0.5;
      float ba6a3norm = sqrt(ba6a3) * 0.5;
      float ba7a3norm = sqrt(ba7a3) * 0.5;

      float ba4a5norm = sqrt(ba4a5) * 0.5;
      float ba4a6norm = sqrt(ba4a6) * 0.5;
      float ba4a7norm = sqrt(ba4a7) * 0.5;

      float ba5a4norm = sqrt(ba5a4) * 0.5;
      float ba6a4norm = sqrt(ba6a4) * 0.5;
      float ba7a4norm = sqrt(ba7a4) * 0.5;

      float ba5a6norm = sqrt(ba5a6) * 0.5;
      float ba5a7norm = sqrt(ba5a7) * 0.5;

      float ba6a5norm = sqrt(ba6a5) * 0.5;
      float ba7a5norm = sqrt(ba7a5) * 0.5;

      float ba6a7norm = sqrt(ba6a7) * 0.5;

      float ba7a6norm = sqrt(ba7a6) * 0.5;

      float cfa0a1 = exp_fast(-ba0a1norm * inv_sigma) * normsuma0a1 * DIV_FACTOR / ba0a1norm;
      float cfa0a2 = exp_fast(-ba0a2norm * inv_sigma) * normsuma0a2 * DIV_FACTOR / ba0a2norm;
      float cfa0a3 = exp_fast(-ba0a3norm * inv_sigma) * normsuma0a3 * DIV_FACTOR / ba0a3norm;
      float cfa0a4 = exp_fast(-ba0a4norm * inv_sigma) * normsuma0a4 * DIV_FACTOR / ba0a4norm;
      float cfa0a5 = exp_fast(-ba0a5norm * inv_sigma) * normsuma0a5 * DIV_FACTOR / ba0a5norm;
      float cfa0a6 = exp_fast(-ba0a6norm * inv_sigma) * normsuma0a6 * DIV_FACTOR / ba0a6norm;
      float cfa0a7 = exp_fast(-ba0a7norm * inv_sigma) * normsuma0a7 * DIV_FACTOR / ba0a7norm;

      float cfa1a0 = exp_fast(-ba1a0norm * inv_sigma) * normsuma1a0 * DIV_FACTOR / ba1a0norm;
      float cfa2a0 = exp_fast(-ba2a0norm * inv_sigma) * normsuma2a0 * DIV_FACTOR / ba2a0norm;
      float cfa3a0 = exp_fast(-ba3a0norm * inv_sigma) * normsuma3a0 * DIV_FACTOR / ba3a0norm;
      float cfa4a0 = exp_fast(-ba4a0norm * inv_sigma) * normsuma4a0 * DIV_FACTOR / ba4a0norm;
      float cfa5a0 = exp_fast(-ba5a0norm * inv_sigma) * normsuma5a0 * DIV_FACTOR / ba5a0norm;
      float cfa6a0 = exp_fast(-ba6a0norm * inv_sigma) * normsuma6a0 * DIV_FACTOR / ba6a0norm;
      float cfa7a0 = exp_fast(-ba7a0norm * inv_sigma) * normsuma7a0 * DIV_FACTOR / ba7a0norm;

      float cfa1a2 = exp_fast(-ba1a2norm * inv_sigma) * normsuma1a2 * DIV_FACTOR / ba1a2norm;
      float cfa1a3 = exp_fast(-ba1a3norm * inv_sigma) * normsuma1a3 * DIV_FACTOR / ba1a3norm;
      float cfa1a4 = exp_fast(-ba1a4norm * inv_sigma) * normsuma1a4 * DIV_FACTOR / ba1a4norm;
      float cfa1a5 = exp_fast(-ba1a5norm * inv_sigma) * normsuma1a5 * DIV_FACTOR / ba1a5norm;
      float cfa1a6 = exp_fast(-ba1a6norm * inv_sigma) * normsuma1a6 * DIV_FACTOR / ba1a6norm;
      float cfa1a7 = exp_fast(-ba1a7norm * inv_sigma) * normsuma1a7 * DIV_FACTOR / ba1a7norm;

      float cfa2a1 = exp_fast(-ba2a1norm * inv_sigma) * normsuma2a1 * DIV_FACTOR / ba2a1norm;
      float cfa3a1 = exp_fast(-ba3a1norm * inv_sigma) * normsuma3a1 * DIV_FACTOR / ba3a1norm;
      float cfa4a1 = exp_fast(-ba4a1norm * inv_sigma) * normsuma4a1 * DIV_FACTOR / ba4a1norm;
      float cfa5a1 = exp_fast(-ba5a1norm * inv_sigma) * normsuma5a1 * DIV_FACTOR / ba5a1norm;
      float cfa6a1 = exp_fast(-ba6a1norm * inv_sigma) * normsuma6a1 * DIV_FACTOR / ba6a1norm;
      float cfa7a1 = exp_fast(-ba7a1norm * inv_sigma) * normsuma7a1 * DIV_FACTOR / ba7a1norm;

      float cfa2a3 = exp_fast(-ba2a3norm * inv_sigma) * normsuma2a3 * DIV_FACTOR / ba2a3norm;
      float cfa2a4 = exp_fast(-ba2a4norm * inv_sigma) * normsuma2a4 * DIV_FACTOR / ba2a4norm;
      float cfa2a5 = exp_fast(-ba2a5norm * inv_sigma) * normsuma2a5 * DIV_FACTOR / ba2a5norm;
      float cfa2a6 = exp_fast(-ba2a6norm * inv_sigma) * normsuma2a6 * DIV_FACTOR / ba2a6norm;
      float cfa2a7 = exp_fast(-ba2a7norm * inv_sigma) * normsuma2a7 * DIV_FACTOR / ba2a7norm;

      float cfa3a2 = exp_fast(-ba3a2norm * inv_sigma) * normsuma3a2 * DIV_FACTOR / ba3a2norm;
      float cfa4a2 = exp_fast(-ba4a2norm * inv_sigma) * normsuma4a2 * DIV_FACTOR / ba4a2norm;
      float cfa5a2 = exp_fast(-ba5a2norm * inv_sigma) * normsuma5a2 * DIV_FACTOR / ba5a2norm;
      float cfa6a2 = exp_fast(-ba6a2norm * inv_sigma) * normsuma6a2 * DIV_FACTOR / ba6a2norm;
      float cfa7a2 = exp_fast(-ba7a2norm * inv_sigma) * normsuma7a2 * DIV_FACTOR / ba7a2norm;

      float cfa3a4 = exp_fast(-ba3a4norm * inv_sigma) * normsuma3a4 * DIV_FACTOR / ba3a4norm;
      float cfa3a5 = exp_fast(-ba3a5norm * inv_sigma) * normsuma3a5 * DIV_FACTOR / ba3a5norm;
      float cfa3a6 = exp_fast(-ba3a6norm * inv_sigma) * normsuma3a6 * DIV_FACTOR / ba3a6norm;
      float cfa3a7 = exp_fast(-ba3a7norm * inv_sigma) * normsuma3a7 * DIV_FACTOR / ba3a7norm;

      float cfa4a3 = exp_fast(-ba4a3norm * inv_sigma) * normsuma4a3 * DIV_FACTOR / ba4a3norm;
      float cfa5a3 = exp_fast(-ba5a3norm * inv_sigma) * normsuma5a3 * DIV_FACTOR / ba5a3norm;
      float cfa6a3 = exp_fast(-ba6a3norm * inv_sigma) * normsuma6a3 * DIV_FACTOR / ba6a3norm;
      float cfa7a3 = exp_fast(-ba7a3norm * inv_sigma) * normsuma7a3 * DIV_FACTOR / ba7a3norm;

      float cfa4a5 = exp_fast(-ba4a5norm * inv_sigma) * normsuma4a5 * DIV_FACTOR / ba4a5norm;
      float cfa4a6 = exp_fast(-ba4a6norm * inv_sigma) * normsuma4a6 * DIV_FACTOR / ba4a6norm;
      float cfa4a7 = exp_fast(-ba4a7norm * inv_sigma) * normsuma4a7 * DIV_FACTOR / ba4a7norm;

      float cfa5a4 = exp_fast(-ba5a4norm * inv_sigma) * normsuma5a4 * DIV_FACTOR / ba5a4norm;
      float cfa6a4 = exp_fast(-ba6a4norm * inv_sigma) * normsuma6a4 * DIV_FACTOR / ba6a4norm;
      float cfa7a4 = exp_fast(-ba7a4norm * inv_sigma) * normsuma7a4 * DIV_FACTOR / ba7a4norm;

      float cfa5a6 = exp_fast(-ba5a6norm * inv_sigma) * normsuma5a6 * DIV_FACTOR / ba5a6norm;
      float cfa5a7 = exp_fast(-ba5a7norm * inv_sigma) * normsuma5a7 * DIV_FACTOR / ba5a7norm;

      float cfa6a5 = exp_fast(-ba6a5norm * inv_sigma) * normsuma6a5 * DIV_FACTOR / ba6a5norm;
      float cfa7a5 = exp_fast(-ba7a5norm * inv_sigma) * normsuma7a5 * DIV_FACTOR / ba7a5norm;

      float cfa6a7 = exp_fast(-ba6a7norm * inv_sigma) * normsuma6a7 * DIV_FACTOR / ba6a7norm;

      float cfa7a6 = exp_fast(-ba7a6norm * inv_sigma) * normsuma7a6 * DIV_FACTOR / ba7a6norm;

      float repxa0a11 = repxa0a1 * cfa0a1;
      float repya0a11 = repya0a1 * cfa0a1;
      float repxa0a21 = repxa0a2 * cfa0a2;
      float repya0a21 = repya0a2 * cfa0a2;
      float repxa0a31 = repxa0a3 * cfa0a3;
      float repya0a31 = repya0a3 * cfa0a3;
      float repxa0a41 = repxa0a4 * cfa0a4;
      float repya0a41 = repya0a4 * cfa0a4;
      float repxa0a51 = repxa0a5 * cfa0a5;
      float repya0a51 = repya0a5 * cfa0a5;
      float repxa0a61 = repxa0a6 * cfa0a6;
      float repya0a61 = repya0a6 * cfa0a6;
      float repxa0a71 = repxa0a7 * cfa0a7;
      float repya0a71 = repya0a7 * cfa0a7;

      float repxa1a01 = repxa1a0 * cfa1a0;
      float repya1a01 = repya1a0 * cfa1a0;
      float repxa2a01 = repxa2a0 * cfa2a0;
      float repya2a01 = repya2a0 * cfa2a0;
      float repxa3a01 = repxa3a0 * cfa3a0;
      float repya3a01 = repya3a0 * cfa3a0;
      float repxa4a01 = repxa4a0 * cfa4a0;
      float repya4a01 = repya4a0 * cfa4a0;
      float repxa5a01 = repxa5a0 * cfa5a0;
      float repya5a01 = repya5a0 * cfa5a0;
      float repxa6a01 = repxa6a0 * cfa6a0;
      float repya6a01 = repya6a0 * cfa6a0;
      float repxa7a01 = repxa7a0 * cfa7a0;
      float repya7a01 = repya7a0 * cfa7a0;

      float repxa1a21 = repxa1a2 * cfa1a2;
      float repya1a21 = repya1a2 * cfa1a2;
      float repxa1a31 = repxa1a3 * cfa1a3;
      float repya1a31 = repya1a3 * cfa1a3;
      float repxa1a41 = repxa1a4 * cfa1a4;
      float repya1a41 = repya1a4 * cfa1a4;
      float repxa1a51 = repxa1a5 * cfa1a5;
      float repya1a51 = repya1a5 * cfa1a5;
      float repxa1a61 = repxa1a6 * cfa1a6;
      float repya1a61 = repya1a6 * cfa1a6;
      float repxa1a71 = repxa1a7 * cfa1a7;
      float repya1a71 = repya1a7 * cfa1a7;

      float repxa2a11 = repxa2a1 * cfa2a1;
      float repya2a11 = repya2a1 * cfa2a1;
      float repxa3a11 = repxa3a1 * cfa3a1;
      float repya3a11 = repya3a1 * cfa3a1;
      float repxa4a11 = repxa4a1 * cfa4a1;
      float repya4a11 = repya4a1 * cfa4a1;
      float repxa5a11 = repxa5a1 * cfa5a1;
      float repya5a11 = repya5a1 * cfa5a1;
      float repxa6a11 = repxa6a1 * cfa6a1;
      float repya6a11 = repya6a1 * cfa6a1;
      float repxa7a11 = repxa7a1 * cfa7a1;
      float repya7a11 = repya7a1 * cfa7a1;

      float repxa2a31 = repxa2a3 * cfa2a3;
      float repya2a31 = repya2a3 * cfa2a3;
      float repxa2a41 = repxa2a4 * cfa2a4;
      float repya2a41 = repya2a4 * cfa2a4;
      float repxa2a51 = repxa2a5 * cfa2a5;
      float repya2a51 = repya2a5 * cfa2a5;
      float repxa2a61 = repxa2a6 * cfa2a6;
      float repya2a61 = repya2a6 * cfa2a6;
      float repxa2a71 = repxa2a7 * cfa2a7;
      float repya2a71 = repya2a7 * cfa2a7;

      float repxa3a21 = repxa3a2 * cfa3a2;
      float repya3a21 = repya3a2 * cfa3a2;
      float repxa4a21 = repxa4a2 * cfa4a2;
      float repya4a21 = repya4a2 * cfa4a2;
      float repxa5a21 = repxa5a2 * cfa5a2;
      float repya5a21 = repya5a2 * cfa5a2;
      float repxa6a21 = repxa6a2 * cfa6a2;
      float repya6a21 = repya6a2 * cfa6a2;
      float repxa7a21 = repxa7a2 * cfa7a2;
      float repya7a21 = repya7a2 * cfa7a2;

      float repxa3a41 = repxa3a4 * cfa3a4;
      float repya3a41 = repya3a4 * cfa3a4;
      float repxa3a51 = repxa3a5 * cfa3a5;
      float repya3a51 = repya3a5 * cfa3a5;
      float repxa3a61 = repxa3a6 * cfa3a6;
      float repya3a61 = repya3a6 * cfa3a6;
      float repxa3a71 = repxa3a7 * cfa3a7;
      float repya3a71 = repya3a7 * cfa3a7;

      float repxa4a31 = repxa4a3 * cfa4a3;
      float repya4a31 = repya4a3 * cfa4a3;
      float repxa5a31 = repxa5a3 * cfa5a3;
      float repya5a31 = repya5a3 * cfa5a3;
      float repxa6a31 = repxa6a3 * cfa6a3;
      float repya6a31 = repya6a3 * cfa6a3;
      float repxa7a31 = repxa7a3 * cfa7a3;
      float repya7a31 = repya7a3 * cfa7a3;

      float repxa4a51 = repxa4a5 * cfa4a5;
      float repya4a51 = repya4a5 * cfa4a5;
      float repxa4a61 = repxa4a6 * cfa4a6;
      float repya4a61 = repya4a6 * cfa4a6;
      float repxa4a71 = repxa4a7 * cfa4a7;
      float repya4a71 = repya4a7 * cfa4a7;

      float repxa5a41 = repxa5a4 * cfa5a4;
      float repya5a41 = repya5a4 * cfa5a4;
      float repxa6a41 = repxa6a4 * cfa6a4;
      float repya6a41 = repya6a4 * cfa6a4;
      float repxa7a41 = repxa7a4 * cfa7a4;
      float repya7a41 = repya7a4 * cfa7a4;

      float repxa5a61 = repxa5a6 * cfa5a6;
      float repya5a61 = repya5a6 * cfa5a6;
      float repxa5a71 = repxa5a7 * cfa5a7;
      float repya5a71 = repya5a7 * cfa5a7;

      float repxa6a51 = repxa6a5 * cfa6a5;
      float repya6a51 = repya6a5 * cfa6a5;
      float repxa7a51 = repxa7a5 * cfa7a5;
      float repya7a51 = repya7a5 * cfa7a5;

      float repxa6a71 = repxa6a7 * cfa6a7;
      float repya6a71 = repya6a7 * cfa6a7;

      float repxa7a61 = repxa7a6 * cfa7a6;
      float repya7a61 = repya7a6 * cfa7a6;

      float ca0a1 = exa0 * repxa0a11 + eya0 * repya0a11;
      float ca0a2 = exa0 * repxa0a21 + eya0 * repya0a21;
      float ca0a3 = exa0 * repxa0a31 + eya0 * repya0a31;
      float ca0a4 = exa0 * repxa0a41 + eya0 * repya0a41;
      float ca0a5 = exa0 * repxa0a51 + eya0 * repya0a51;
      float ca0a6 = exa0 * repxa0a61 + eya0 * repya0a61;
      float ca0a7 = exa0 * repxa0a71 + eya0 * repya0a71;

      float ca1a0 = exa1 * repxa1a01 + eya1 * repya1a01;
      float ca2a0 = exa2 * repxa2a01 + eya2 * repya2a01;
      float ca3a0 = exa3 * repxa3a01 + eya3 * repya3a01;
      float ca4a0 = exa4 * repxa4a01 + eya4 * repya4a01;
      float ca5a0 = exa5 * repxa5a01 + eya5 * repya5a01;
      float ca6a0 = exa6 * repxa6a01 + eya6 * repya6a01;
      float ca7a0 = exa7 * repxa7a01 + eya7 * repya7a01;

      float ca1a2 = exa1 * repxa1a21 + eya1 * repya1a21;
      float ca1a3 = exa1 * repxa1a31 + eya1 * repya1a31;
      float ca1a4 = exa1 * repxa1a41 + eya1 * repya1a41;
      float ca1a5 = exa1 * repxa1a51 + eya1 * repya1a51;
      float ca1a6 = exa1 * repxa1a61 + eya1 * repya1a61;
      float ca1a7 = exa1 * repxa1a71 + eya1 * repya1a71;

      float ca2a1 = exa2 * repxa2a11 + eya2 * repya2a11;
      float ca3a1 = exa3 * repxa3a11 + eya3 * repya3a11;
      float ca4a1 = exa4 * repxa4a11 + eya4 * repya4a11;
      float ca5a1 = exa5 * repxa5a11 + eya5 * repya5a11;
      float ca6a1 = exa6 * repxa6a11 + eya6 * repya6a11;
      float ca7a1 = exa7 * repxa7a11 + eya7 * repya7a11;

      float ca2a3 = exa2 * repxa2a31 + eya2 * repya2a31;
      float ca2a4 = exa2 * repxa2a41 + eya2 * repya2a41;
      float ca2a5 = exa2 * repxa2a51 + eya2 * repya2a51;
      float ca2a6 = exa2 * repxa2a61 + eya2 * repya2a61;
      float ca2a7 = exa2 * repxa2a71 + eya2 * repya2a71;

      float ca3a2 = exa3 * repxa3a21 + eya3 * repya3a21;
      float ca4a2 = exa4 * repxa4a21 + eya4 * repya4a21;
      float ca5a2 = exa5 * repxa5a21 + eya5 * repya5a21;
      float ca6a2 = exa6 * repxa6a21 + eya6 * repya6a21;
      float ca7a2 = exa7 * repxa7a21 + eya7 * repya7a21;

      float ca3a4 = exa3 * repxa3a41 + eya3 * repya3a41;
      float ca3a5 = exa3 * repxa3a51 + eya3 * repya3a51;
      float ca3a6 = exa3 * repxa3a61 + eya3 * repya3a61;
      float ca3a7 = exa3 * repxa3a71 + eya3 * repya3a71;

      float ca4a3 = exa4 * repxa4a31 + eya4 * repya4a31;
      float ca5a3 = exa5 * repxa5a31 + eya5 * repya5a31;
      float ca6a3 = exa6 * repxa6a31 + eya6 * repya6a31;
      float ca7a3 = exa7 * repxa7a31 + eya7 * repya7a31;

      float ca4a5 = exa4 * repxa4a51 + eya4 * repya4a51;
      float ca4a6 = exa4 * repxa4a61 + eya4 * repya4a61;
      float ca4a7 = exa4 * repxa4a71 + eya4 * repya4a71;

      float ca5a4 = exa5 * repxa5a41 + eya5 * repya5a41;
      float ca6a4 = exa6 * repxa6a41 + eya6 * repya6a41;
      float ca7a4 = exa7 * repxa7a41 + eya7 * repya7a41;

      float ca5a6 = exa5 * repxa5a61 + eya5 * repya5a61;
      float ca5a7 = exa5 * repxa5a71 + eya5 * repya5a71;

      float ca6a5 = exa6 * repxa6a51 + eya6 * repya6a51;
      float ca7a5 = exa7 * repxa7a51 + eya7 * repya7a51;

      float ca6a7 = exa6 * repxa6a71 + eya6 * repya6a71;

      float ca7a6 = exa7 * repxa7a61 + eya7 * repya7a61;

      float tha0a1 = sqrt(repxa0a11 * repxa0a11 + repya0a11 * repya0a11) * cospsi;
      float tha0a2 = sqrt(repxa0a21 * repxa0a21 + repya0a21 * repya0a21) * cospsi;
      float tha0a3 = sqrt(repxa0a31 * repxa0a31 + repya0a31 * repya0a31) * cospsi;
      float tha0a4 = sqrt(repxa0a41 * repxa0a41 + repya0a41 * repya0a41) * cospsi;
      float tha0a5 = sqrt(repxa0a51 * repxa0a51 + repya0a51 * repya0a51) * cospsi;
      float tha0a6 = sqrt(repxa0a61 * repxa0a61 + repya0a61 * repya0a61) * cospsi;
      float tha0a7 = sqrt(repxa0a71 * repxa0a71 + repya0a71 * repya0a71) * cospsi;

      float tha1a0 = sqrt(repxa1a01 * repxa1a01 + repya1a01 * repya1a01) * cospsi;
      float tha2a0 = sqrt(repxa2a01 * repxa2a01 + repya2a01 * repya2a01) * cospsi;
      float tha3a0 = sqrt(repxa3a01 * repxa3a01 + repya3a01 * repya3a01) * cospsi;
      float tha4a0 = sqrt(repxa4a01 * repxa4a01 + repya4a01 * repya4a01) * cospsi;
      float tha5a0 = sqrt(repxa5a01 * repxa5a01 + repya5a01 * repya5a01) * cospsi;
      float tha6a0 = sqrt(repxa6a01 * repxa6a01 + repya6a01 * repya6a01) * cospsi;
      float tha7a0 = sqrt(repxa7a01 * repxa7a01 + repya7a01 * repya7a01) * cospsi;

      float tha1a2 = sqrt(repxa1a21 * repxa1a21 + repya1a21 * repya1a21) * cospsi;
      float tha1a3 = sqrt(repxa1a31 * repxa1a31 + repya1a31 * repya1a31) * cospsi;
      float tha1a4 = sqrt(repxa1a41 * repxa1a41 + repya1a41 * repya1a41) * cospsi;
      float tha1a5 = sqrt(repxa1a51 * repxa1a51 + repya1a51 * repya1a51) * cospsi;
      float tha1a6 = sqrt(repxa1a61 * repxa1a61 + repya1a61 * repya1a61) * cospsi;
      float tha1a7 = sqrt(repxa1a71 * repxa1a71 + repya1a71 * repya1a71) * cospsi;

      float tha2a1 = sqrt(repxa2a11 * repxa2a11 + repya2a11 * repya2a11) * cospsi;
      float tha3a1 = sqrt(repxa3a11 * repxa3a11 + repya3a11 * repya3a11) * cospsi;
      float tha4a1 = sqrt(repxa4a11 * repxa4a11 + repya4a11 * repya4a11) * cospsi;
      float tha5a1 = sqrt(repxa5a11 * repxa5a11 + repya5a11 * repya5a11) * cospsi;
      float tha6a1 = sqrt(repxa6a11 * repxa6a11 + repya6a11 * repya6a11) * cospsi;
      float tha7a1 = sqrt(repxa7a11 * repxa7a11 + repya7a11 * repya7a11) * cospsi;

      float tha2a3 = sqrt(repxa2a31 * repxa2a31 + repya2a31 * repya2a31) * cospsi;
      float tha2a4 = sqrt(repxa2a41 * repxa2a41 + repya2a41 * repya2a41) * cospsi;
      float tha2a5 = sqrt(repxa2a51 * repxa2a51 + repya2a51 * repya2a51) * cospsi;
      float tha2a6 = sqrt(repxa2a61 * repxa2a61 + repya2a61 * repya2a61) * cospsi;
      float tha2a7 = sqrt(repxa2a71 * repxa2a71 + repya2a71 * repya2a71) * cospsi;

      float tha3a2 = sqrt(repxa3a21 * repxa3a21 + repya3a21 * repya3a21) * cospsi;
      float tha4a2 = sqrt(repxa4a21 * repxa4a21 + repya4a21 * repya4a21) * cospsi;
      float tha5a2 = sqrt(repxa5a21 * repxa5a21 + repya5a21 * repya5a21) * cospsi;
      float tha6a2 = sqrt(repxa6a21 * repxa6a21 + repya6a21 * repya6a21) * cospsi;
      float tha7a2 = sqrt(repxa7a21 * repxa7a21 + repya7a21 * repya7a21) * cospsi;

      float tha3a4 = sqrt(repxa3a41 * repxa3a41 + repya3a41 * repya3a41) * cospsi;
      float tha3a5 = sqrt(repxa3a51 * repxa3a51 + repya3a51 * repya3a51) * cospsi;
      float tha3a6 = sqrt(repxa3a61 * repxa3a61 + repya3a61 * repya3a61) * cospsi;
      float tha3a7 = sqrt(repxa3a71 * repxa3a71 + repya3a71 * repya3a71) * cospsi;

      float tha4a3 = sqrt(repxa4a31 * repxa4a31 + repya4a31 * repya4a31) * cospsi;
      float tha5a3 = sqrt(repxa5a31 * repxa5a31 + repya5a31 * repya5a31) * cospsi;
      float tha6a3 = sqrt(repxa6a31 * repxa6a31 + repya6a31 * repya6a31) * cospsi;
      float tha7a3 = sqrt(repxa7a31 * repxa7a31 + repya7a31 * repya7a31) * cospsi;

      float tha4a5 = sqrt(repxa4a51 * repxa4a51 + repya4a51 * repya4a51) * cospsi;
      float tha4a6 = sqrt(repxa4a61 * repxa4a61 + repya4a61 * repya4a61) * cospsi;
      float tha4a7 = sqrt(repxa4a71 * repxa4a71 + repya4a71 * repya4a71) * cospsi;

      float tha5a4 = sqrt(repxa5a41 * repxa5a41 + repya5a41 * repya5a41) * cospsi;
      float tha6a4 = sqrt(repxa6a41 * repxa6a41 + repya6a41 * repya6a41) * cospsi;
      float tha7a4 = sqrt(repxa7a41 * repxa7a41 + repya7a41 * repya7a41) * cospsi;

      float tha5a6 = sqrt(repxa5a61 * repxa5a61 + repya5a61 * repya5a61) * cospsi;
      float tha5a7 = sqrt(repxa5a71 * repxa5a71 + repya5a71 * repya5a71) * cospsi;

      float tha6a5 = sqrt(repxa6a51 * repxa6a51 + repya6a51 * repya6a51) * cospsi;
      float tha7a5 = sqrt(repxa7a51 * repxa7a51 + repya7a51 * repya7a51) * cospsi;

      float tha6a7 = sqrt(repxa6a71 * repxa6a71 + repya6a71 * repya6a71) * cospsi;

      float tha7a6 = sqrt(repxa7a61 * repxa7a61 + repya7a61 * repya7a61) * cospsi;

      float wa0a1 = -ca0a1 >= tha0a1 ? 1 : INFLUENCE;
      float wa0a2 = -ca0a2 >= tha0a2 ? 1 : INFLUENCE;
      float wa0a3 = -ca0a3 >= tha0a3 ? 1 : INFLUENCE;
      float wa0a4 = -ca0a4 >= tha0a4 ? 1 : INFLUENCE;
      float wa0a5 = -ca0a5 >= tha0a5 ? 1 : INFLUENCE;
      float wa0a6 = -ca0a6 >= tha0a6 ? 1 : INFLUENCE;
      float wa0a7 = -ca0a7 >= tha0a7 ? 1 : INFLUENCE;

      float wa1a0 = -ca1a0 >= tha1a0 ? 1 : INFLUENCE;
      float wa2a0 = -ca2a0 >= tha2a0 ? 1 : INFLUENCE;
      float wa3a0 = -ca3a0 >= tha3a0 ? 1 : INFLUENCE;
      float wa4a0 = -ca4a0 >= tha4a0 ? 1 : INFLUENCE;
      float wa5a0 = -ca5a0 >= tha5a0 ? 1 : INFLUENCE;
      float wa6a0 = -ca6a0 >= tha6a0 ? 1 : INFLUENCE;
      float wa7a0 = -ca7a0 >= tha7a0 ? 1 : INFLUENCE;

      float wa1a2 = -ca1a2 >= tha1a2 ? 1 : INFLUENCE;
      float wa1a3 = -ca1a3 >= tha1a3 ? 1 : INFLUENCE;
      float wa1a4 = -ca1a4 >= tha1a4 ? 1 : INFLUENCE;
      float wa1a5 = -ca1a5 >= tha1a5 ? 1 : INFLUENCE;
      float wa1a6 = -ca1a6 >= tha1a6 ? 1 : INFLUENCE;
      float wa1a7 = -ca1a7 >= tha1a7 ? 1 : INFLUENCE;

      float wa2a1 = -ca2a1 >= tha2a1 ? 1 : INFLUENCE;
      float wa3a1 = -ca3a1 >= tha3a1 ? 1 : INFLUENCE;
      float wa4a1 = -ca4a1 >= tha4a1 ? 1 : INFLUENCE;
      float wa5a1 = -ca5a1 >= tha5a1 ? 1 : INFLUENCE;
      float wa6a1 = -ca6a1 >= tha6a1 ? 1 : INFLUENCE;
      float wa7a1 = -ca7a1 >= tha7a1 ? 1 : INFLUENCE;

      float wa2a3 = -ca2a3 >= tha2a3 ? 1 : INFLUENCE;
      float wa2a4 = -ca2a4 >= tha2a4 ? 1 : INFLUENCE;
      float wa2a5 = -ca2a5 >= tha2a5 ? 1 : INFLUENCE;
      float wa2a6 = -ca2a6 >= tha2a6 ? 1 : INFLUENCE;
      float wa2a7 = -ca2a7 >= tha2a7 ? 1 : INFLUENCE;

      float wa3a2 = -ca3a2 >= tha3a2 ? 1 : INFLUENCE;
      float wa4a2 = -ca4a2 >= tha4a2 ? 1 : INFLUENCE;
      float wa5a2 = -ca5a2 >= tha5a2 ? 1 : INFLUENCE;
      float wa6a2 = -ca6a2 >= tha6a2 ? 1 : INFLUENCE;
      float wa7a2 = -ca7a2 >= tha7a2 ? 1 : INFLUENCE;

      float wa3a4 = -ca3a4 >= tha3a4 ? 1 : INFLUENCE;
      float wa3a5 = -ca3a5 >= tha3a5 ? 1 : INFLUENCE;
      float wa3a6 = -ca3a6 >= tha3a6 ? 1 : INFLUENCE;
      float wa3a7 = -ca3a7 >= tha3a7 ? 1 : INFLUENCE;

      float wa4a3 = -ca4a3 >= tha4a3 ? 1 : INFLUENCE;
      float wa5a3 = -ca5a3 >= tha5a3 ? 1 : INFLUENCE;
      float wa6a3 = -ca6a3 >= tha6a3 ? 1 : INFLUENCE;
      float wa7a3 = -ca7a3 >= tha7a3 ? 1 : INFLUENCE;

      float wa4a5 = -ca4a5 >= tha4a5 ? 1 : INFLUENCE;
      float wa4a6 = -ca4a6 >= tha4a6 ? 1 : INFLUENCE;
      float wa4a7 = -ca4a7 >= tha4a7 ? 1 : INFLUENCE;

      float wa5a4 = -ca5a4 >= tha5a4 ? 1 : INFLUENCE;
      float wa6a4 = -ca6a4 >= tha6a4 ? 1 : INFLUENCE;
      float wa7a4 = -ca7a4 >= tha7a4 ? 1 : INFLUENCE;

      float wa5a6 = -ca5a6 >= tha5a6 ? 1 : INFLUENCE;
      float wa5a7 = -ca5a7 >= tha5a7 ? 1 : INFLUENCE;

      float wa6a5 = -ca6a5 >= tha6a5 ? 1 : INFLUENCE;
      float wa7a5 = -ca7a5 >= tha7a5 ? 1 : INFLUENCE;

      float wa6a7 = -ca6a7 >= tha6a7 ? 1 : INFLUENCE;

      float wa7a6 = -ca7a6 >= tha7a6 ? 1 : INFLUENCE;

      sfx0 += (repxa0a11 * wa0a1) + (repxa0a21 * wa0a2) + (repxa0a31 * wa0a3) + (repxa0a41 * wa0a4) + (repxa0a51 * wa0a5) + (repxa0a61 * wa0a6) + (repxa0a71 * wa0a7);
      sfy0 += (repya0a11 * wa0a1) + (repya0a21 * wa0a2) + (repya0a31 * wa0a3) + (repya0a41 * wa0a4) + (repya0a51 * wa0a5) + (repya0a61 * wa0a6) + (repya0a71 * wa0a7);
      sfx1 += (repxa1a01 * wa1a0) + (repxa1a21 * wa1a2) + (repxa1a31 * wa1a3) + (repxa1a41 * wa1a4) + (repxa1a51 * wa1a5) + (repxa1a61 * wa1a6) + (repxa1a71 * wa1a7);
      sfy1 += (repya1a01 * wa1a0) + (repya1a21 * wa1a2) + (repya1a31 * wa1a3) + (repya1a41 * wa1a4) + (repya1a51 * wa1a5) + (repya1a61 * wa1a6) + (repya1a71 * wa1a7);
      sfx2 += (repxa2a11 * wa2a1) + (repxa2a01 * wa2a0) + (repxa2a31 * wa2a3) + (repxa2a41 * wa2a4) + (repxa2a51 * wa2a5) + (repxa2a61 * wa2a6) + (repxa2a71 * wa2a7);
      sfy2 += (repya2a11 * wa2a1) + (repya2a01 * wa2a0) + (repya2a31 * wa2a3) + (repya2a41 * wa2a4) + (repya2a51 * wa2a5) + (repya2a61 * wa2a6) + (repya2a71 * wa2a7);
      sfx3 += (repxa3a11 * wa3a1) + (repxa3a21 * wa3a2) + (repxa3a01 * wa3a0) + (repxa3a41 * wa3a4) + (repxa3a51 * wa3a5) + (repxa3a61 * wa3a6) + (repxa3a71 * wa3a7);
      sfy3 += (repya3a11 * wa3a1) + (repya3a21 * wa3a2) + (repya3a01 * wa3a0) + (repya3a41 * wa3a4) + (repya3a51 * wa3a5) + (repya3a61 * wa3a6) + (repya3a71 * wa3a7);
      sfx4 += (repxa4a11 * wa4a1) + (repxa4a21 * wa4a2) + (repxa4a31 * wa4a3) + (repxa4a01 * wa4a0) + (repxa4a51 * wa4a5) + (repxa4a61 * wa4a6) + (repxa4a71 * wa4a7);
      sfy4 += (repya4a11 * wa4a1) + (repya4a21 * wa4a2) + (repya4a31 * wa4a3) + (repya4a01 * wa4a0) + (repya4a51 * wa4a5) + (repya4a61 * wa4a6) + (repya4a71 * wa4a7);
      sfx5 += (repxa5a11 * wa5a1) + (repxa5a21 * wa5a2) + (repxa5a31 * wa5a3) + (repxa5a41 * wa5a4) + (repxa5a01 * wa5a0) + (repxa5a61 * wa5a6) + (repxa5a71 * wa5a7);
      sfy5 += (repya5a11 * wa5a1) + (repya5a21 * wa5a2) + (repya5a31 * wa5a3) + (repya5a41 * wa5a4) + (repya5a01 * wa5a0) + (repya5a61 * wa5a6) + (repya5a71 * wa5a7);
      sfx6 += (repxa6a11 * wa6a1) + (repxa6a21 * wa6a2) + (repxa6a31 * wa6a3) + (repxa6a41 * wa6a4) + (repxa6a51 * wa6a5) + (repxa6a01 * wa6a0) + (repxa6a71 * wa6a7);
      sfy6 += (repya6a11 * wa6a1) + (repya6a21 * wa6a2) + (repya6a31 * wa6a3) + (repya6a41 * wa6a4) + (repya6a51 * wa6a5) + (repya6a01 * wa6a0) + (repya6a71 * wa6a7);
      sfx7 += (repxa7a11 * wa7a1) + (repxa7a21 * wa7a2) + (repxa7a31 * wa7a3) + (repxa7a41 * wa7a4) + (repxa7a51 * wa7a5) + (repxa7a61 * wa7a6) + (repxa7a01 * wa7a0);
      sfy7 += (repya7a11 * wa7a1) + (repya7a21 * wa7a2) + (repya7a31 * wa7a3) + (repya7a41 * wa7a4) + (repya7a51 * wa7a5) + (repya7a61 * wa7a6) + (repya7a01 * wa7a0);

      social_force[IndexX(i)] = sfx0;
      social_force[IndexY(i, n)] = sfy0;
      social_force[IndexX(i + 1)] = sfx1;
      social_force[IndexY(i + 1, n)] = sfy1;
      social_force[IndexX(i + 2)] = sfx2;
      social_force[IndexY(i + 2, n)] = sfy2;
      social_force[IndexX(i + 3)] = sfx3;
      social_force[IndexY(i + 3, n)] = sfy3;
      social_force[IndexX(i + 4)] = sfx4;
      social_force[IndexY(i + 4, n)] = sfy4;
      social_force[IndexX(i + 5)] = sfx5;
      social_force[IndexY(i + 5, n)] = sfy5;
      social_force[IndexX(i + 6)] = sfx6;
      social_force[IndexY(i + 6, n)] = sfy6;
      social_force[IndexX(i + 7)] = sfx7;
      social_force[IndexY(i + 7, n)] = sfy7;
    }

    for (int i = 0; i < n - 7; i += 8)
    {

      /************************************************/
      // LOADS
      /************************************************/
      float sfx0 = social_force[IndexX(i)]; //social force x
      float sfy0 = social_force[IndexY(i, n)];
      float sfx1 = social_force[IndexX(i + 1)];
      float sfy1 = social_force[IndexY(i + 1, n)];
      float sfx2 = social_force[IndexX(i + 2)];
      float sfy2 = social_force[IndexY(i + 2, n)];
      float sfx3 = social_force[IndexX(i + 3)];
      float sfy3 = social_force[IndexY(i + 3, n)];
      float sfx4 = social_force[IndexX(i + 4)];
      float sfy4 = social_force[IndexY(i + 4, n)];
      float sfx5 = social_force[IndexX(i + 5)];
      float sfy5 = social_force[IndexY(i + 5, n)];
      float sfx6 = social_force[IndexX(i + 6)];
      float sfy6 = social_force[IndexY(i + 6, n)];
      float sfx7 = social_force[IndexX(i + 7)];
      float sfy7 = social_force[IndexY(i + 7, n)];

      float rxa0 = position[IndexX(i)];
      float rya0 = position[IndexY(i, n)];
      float rxa1 = position[IndexX(i + 1)];
      float rya1 = position[IndexY(i + 1, n)];
      float rxa2 = position[IndexX(i + 2)];
      float rya2 = position[IndexY(i + 2, n)];
      float rxa3 = position[IndexX(i + 3)];
      float rya3 = position[IndexY(i + 3, n)];
      float rxa4 = position[IndexX(i + 4)];
      float rya4 = position[IndexY(i + 4, n)];
      float rxa5 = position[IndexX(i + 5)];
      float rya5 = position[IndexY(i + 5, n)];
      float rxa6 = position[IndexX(i + 6)];
      float rya6 = position[IndexY(i + 6, n)];
      float rxa7 = position[IndexX(i + 7)];
      float rya7 = position[IndexY(i + 7, n)];

      float exa0 = desired_direction[IndexX(i)];
      float eya0 = desired_direction[IndexY(i, n)];
      float exa1 = desired_direction[IndexX(i + 1)];
      float eya1 = desired_direction[IndexY(i + 1, n)];
      float exa2 = desired_direction[IndexX(i + 2)];
      float eya2 = desired_direction[IndexY(i + 2, n)];
      float exa3 = desired_direction[IndexX(i + 3)];
      float eya3 = desired_direction[IndexY(i + 3, n)];
      float exa4 = desired_direction[IndexX(i + 4)];
      float eya4 = desired_direction[IndexY(i + 4, n)];
      float exa5 = desired_direction[IndexX(i + 5)];
      float eya5 = desired_direction[IndexY(i + 5, n)];
      float exa6 = desired_direction[IndexX(i + 6)];
      float eya6 = desired_direction[IndexY(i + 6, n)];
      float exa7 = desired_direction[IndexX(i + 7)];
      float eya7 = desired_direction[IndexY(i + 7, n)];

      float avx0 = actual_velocity[IndexX(i)];
      float avy0 = actual_velocity[IndexY(i, n)];
      float avx1 = actual_velocity[IndexX(i + 1)];
      float avy1 = actual_velocity[IndexY(i + 1, n)];
      float avx2 = actual_velocity[IndexX(i + 2)];
      float avy2 = actual_velocity[IndexY(i + 2, n)];
      float avx3 = actual_velocity[IndexX(i + 3)];
      float avy3 = actual_velocity[IndexY(i + 3, n)];
      float avx4 = actual_velocity[IndexX(i + 4)];
      float avy4 = actual_velocity[IndexY(i + 4, n)];
      float avx5 = actual_velocity[IndexX(i + 5)];
      float avy5 = actual_velocity[IndexY(i + 5, n)];
      float avx6 = actual_velocity[IndexX(i + 6)];
      float avy6 = actual_velocity[IndexY(i + 6, n)];
      float avx7 = actual_velocity[IndexX(i + 7)];
      float avy7 = actual_velocity[IndexY(i + 7, n)];

      float dsv0 = desired_speed[i];     //desired speed value
      float dsv1 = desired_speed[i + 1]; //desired speed value
      float dsv2 = desired_speed[i + 2]; //desired speed value
      float dsv3 = desired_speed[i + 3]; //desired speed value
      float dsv4 = desired_speed[i + 4]; //desired speed value
      float dsv5 = desired_speed[i + 5]; //desired speed value
      float dsv6 = desired_speed[i + 6]; //desired speed value
      float dsv7 = desired_speed[i + 7]; //desired speed value

      float da0 = speed[i];
      float da1 = speed[i + 1];
      float da2 = speed[i + 2];
      float da3 = speed[i + 3];
      float da4 = speed[i + 4];
      float da5 = speed[i + 5];
      float da6 = speed[i + 6];
      float da7 = speed[i + 7];

      __m256 rxa = _mm256_set_ps(rxa7, rxa6, rxa5, rxa4, rxa3, rxa2, rxa1, rxa0);
      __m256 rya = _mm256_set_ps(rya7, rya6, rya5, rya4, rya3, rya2, rya1, rya0);
      __m256 one = _mm256_set1_ps(1.0);
      __m256 minus_one = _mm256_set1_ps(-1.0);
      __m256 half = _mm256_set1_ps(0.5);
      __m256 exp_constant = _mm256_set1_ps(0.00006103515); // 1 / 16384
      __m256 inv_sigma_vec = _mm256_set1_ps(-inv_sigma);
      __m256 div_factor_vec = _mm256_set1_ps(DIV_FACTOR);
      __m256 cospsi_vec = _mm256_set1_ps(cospsi);
      __m256 influencer_vec = _mm256_set1_ps(INFLUENCE);
      __m256 da = _mm256_set_ps(da7, da6, da5, da4, da3, da2, da1, da0);
      __m256 exa = _mm256_set_ps(exa7, exa6, exa5, exa4, exa3, exa2, exa1, exa0);
      __m256 eya = _mm256_set_ps(eya7, eya6, eya5, eya4, eya3, eya2, eya1, eya0);
      __m256 sfx = _mm256_setzero_ps();
      __m256 sfy = _mm256_setzero_ps();

      //iterate over all people
      for (int j = i + 8; j < n; j++) // for (int j = i + 1; j < n; j++)
      {

        __m256 rxb_vec = _mm256_broadcast_ss(position + IndexX(j));
        __m256 ryb_vec = _mm256_broadcast_ss(position + IndexY(j, n));
        __m256 exb0_vec = _mm256_broadcast_ss(desired_direction + IndexX(j));
        __m256 eyb0_vec = _mm256_broadcast_ss(desired_direction + IndexY(j, n));
        __m256 db0_vec = _mm256_broadcast_ss(speed + IndexX(j));

        db0_vec = _mm256_mul_ps(db0_vec, minus_one);
        __m256 da_vec = _mm256_sub_ps(_mm256_setzero_ps(), da);

        __m256 rxab = _mm256_sub_ps(rxa, rxb_vec);
        __m256 ryab = _mm256_sub_ps(rya, ryb_vec);

        __m256 rxab_2 = _mm256_mul_ps(rxab, rxab);
        __m256 rabnorm = _mm256_fmadd_ps(ryab, ryab, rxab_2);
        rabnorm = _mm256_sqrt_ps(rabnorm);

        //float everything
        __m256 rxabme = _mm256_fmadd_ps(db0_vec, exb0_vec, rxab);
        __m256 ryabme = _mm256_fmadd_ps(db0_vec, eyb0_vec, ryab);

        __m256 rxbame = _mm256_fmsub_ps(da_vec, exa, rxab);
        __m256 rybame = _mm256_fmsub_ps(da_vec, eya, ryab);

        // printf("da %f %f\n", da0, da_vec[0]);

        // printf("rxab %f %f\n", rxab0, rxab[0]);
        // printf("ryab %f %f\n", ryab0, ryab[0]);

        // printf("rxabme %f %f\n", rxabmex0, rxabme[0]);
        // printf("ryabme %f %f\n", ryabmey0, ryabme[0]);

        // printf("rxbame %f %f\n", rxbamex0, rxbame[0]);
        // printf("rybame %f %f\n", rybamey0, rybame[0]);

        __m256 rxabme_2 = _mm256_mul_ps(rxabme, rxabme);
        __m256 rabmenorm = _mm256_fmadd_ps(ryabme, ryabme, rxabme_2);
        rabmenorm = _mm256_sqrt_ps(rabmenorm);

        __m256 rxbame_2 = _mm256_mul_ps(rxbame, rxbame);
        __m256 rbamenorm = _mm256_fmadd_ps(rybame, rybame, rxbame_2);
        rbamenorm = _mm256_sqrt_ps(rbamenorm);

        __m256 normsumab = _mm256_add_ps(rabnorm, rabmenorm);
        __m256 normsumba = _mm256_add_ps(rabnorm, rbamenorm); // maybe wrong use rbanorm

        // printf("rabmenorm %f %f\n", rabmenorm0, rabmenorm[0]);
        // printf("rbamenorm %f %f\n", rbamenorm0, rbamenorm[0]);

        __m256 div_1x = _mm256_div_ps(rxab, rabnorm);
        __m256 div_1y = _mm256_div_ps(ryab, rabnorm);
        __m256 div_2x = _mm256_div_ps(rxabme, rabmenorm);
        __m256 div_2y = _mm256_div_ps(ryabme, rabmenorm);

        __m256 repxab = _mm256_add_ps(div_1x, div_2x);
        __m256 repyab = _mm256_add_ps(div_1y, div_2y);

        __m256 div_2x_me = _mm256_div_ps(rxbame, rbamenorm);
        __m256 div_2y_me = _mm256_div_ps(rybame, rbamenorm);

        __m256 repxba = _mm256_sub_ps(div_2x_me, div_1x);
        __m256 repyba = _mm256_sub_ps(div_2y_me, div_1y);

        // printf("repxab00 %f %f\n", repxab00, repxab[0]);

        __m256 db_2 = _mm256_mul_ps(db0_vec, db0_vec);
        __m256 da_2 = _mm256_mul_ps(da_vec, da_vec);
        __m256 bab = _mm256_fmsub_ps(normsumab, normsumab, db_2);
        __m256 bba = _mm256_fmsub_ps(normsumba, normsumba, da_2);

        // printf("da2 %f %f\n", da0 * da0, da_2[0]);
        // printf("normsumab %f %f\n", normsumab0, normsumab[0]);
        // printf("normsumba %f %f\n", normsumba0, normsumba[0]);

        // printf("bab %f %f\n", bab0, bab[0]);
        // printf("bba %f %f\n", bba0, bba[0]);

        __m256 babnorm = _mm256_sqrt_ps(bab);
        __m256 bbanorm = _mm256_sqrt_ps(bba);
        babnorm = _mm256_mul_ps(babnorm, half);
        bbanorm = _mm256_mul_ps(bbanorm, half);

        __m256 cfab = exp_fast_vec_2_5_1(_mm256_mul_ps(babnorm, inv_sigma_vec), one, exp_constant);
        cfab = _mm256_mul_ps(cfab, div_factor_vec);
        cfab = _mm256_mul_ps(cfab, normsumab);
        cfab = _mm256_div_ps(cfab, babnorm);

        __m256 cfba = exp_fast_vec_2_5_1(_mm256_mul_ps(bbanorm, inv_sigma_vec), one, exp_constant);
        cfba = _mm256_mul_ps(cfba, div_factor_vec);
        cfba = _mm256_mul_ps(cfba, normsumba);
        cfba = _mm256_div_ps(cfba, bbanorm);

        repxab = _mm256_mul_ps(repxab, cfab);
        repyab = _mm256_mul_ps(repyab, cfab);
        repxba = _mm256_mul_ps(repxba, cfba);
        repyba = _mm256_mul_ps(repyba, cfba);

        // printf("cfab %f %f\n", cfab0, cfab[0]);
        // printf("cfba %f %f\n", cfba0, cfba[0]);

        // printf("babnorm %f %f\n", babnorm0, babnorm[0]);
        // printf("bbanorm %f %f\n", bbanorm0, bbanorm[0]);

        __m256 cab = _mm256_mul_ps(exa, repxab);
        cab = _mm256_fmadd_ps(eya, repyab, cab);

        __m256 cba = _mm256_mul_ps(exb0_vec, repxba);
        cba = _mm256_fmadd_ps(eyb0_vec, repyba, cba);

        // printf("cab %f %f\n", cab0, cab[0]);
        // printf("cba %f %f\n", cba0, cba[0]);

        __m256 repxab_2 = _mm256_mul_ps(repxab, repxab);
        __m256 thab = _mm256_sqrt_ps(_mm256_fmadd_ps(repyab, repyab, repxab_2));
        thab = _mm256_mul_ps(thab, cospsi_vec);

        __m256 repxba_2 = _mm256_mul_ps(repxba, repxba);
        __m256 thba = _mm256_sqrt_ps(_mm256_fmadd_ps(repyba, repyba, repxba_2));
        thba = _mm256_mul_ps(thba, cospsi_vec);

        // printf("thab %f %f\n", thab0, thab[0]);
        // printf("thba %f %f\n", thba0, thba[0]);

        __m256 maskab = _mm256_cmp_ps(_mm256_mul_ps(cab, minus_one), thab, _CMP_GE_OQ);
        __m256 wab = _mm256_blendv_ps(influencer_vec, one, maskab);

        __m256 maskba = _mm256_cmp_ps(_mm256_mul_ps(cba, minus_one), thba, _CMP_GE_OQ);
        __m256 wba = _mm256_blendv_ps(influencer_vec, one, maskba);

        sfx = _mm256_fmadd_ps(repxab, wab, sfx);
        sfy = _mm256_fmadd_ps(repyab, wab, sfy);

        __m256 temp_x = _mm256_mul_ps(repxba, wba);
        __m256 temp_y = _mm256_mul_ps(repyba, wba);

        temp_x = _mm256_hadd_ps(temp_x, _mm256_permute2f128_ps(temp_x, temp_x, 3));
        temp_x = _mm256_hadd_ps(temp_x, temp_x);
        temp_x = _mm256_hadd_ps(temp_x, temp_x);

        temp_y = _mm256_hadd_ps(temp_y, _mm256_permute2f128_ps(temp_y, temp_y, 3));
        temp_y = _mm256_hadd_ps(temp_y, temp_y);
        temp_y = _mm256_hadd_ps(temp_y, temp_y);

        social_force[IndexX(j)] += temp_x[0];
        social_force[IndexY(j, n)] += temp_y[0];
      } // n-1 * (12 adds, 18 mults, 6 divs, 1 exp, 4 sqrts)

      // add stuff to sfx
      sfx0 += sfx[0];
      sfy0 += sfy[0];
      sfx1 += sfx[1];
      sfy1 += sfy[1];
      sfx2 += sfx[2];
      sfy2 += sfy[2];
      sfx3 += sfx[3];
      sfy3 += sfy[3];
      sfx4 += sfx[4];
      sfy4 += sfy[4];
      sfx5 += sfx[5];
      sfy5 += sfy[5];
      sfx6 += sfx[6];
      sfy6 += sfy[6];
      sfx7 += sfx[7];
      sfy7 += sfy[7];

      /************************************************/
      //UPDATE ACCELERATION TERM
      /************************************************/
      // get actual velocity, desired direction, desired speed

      // compute velocity difference
      float vdx00 = dsv0 * exa0; // 1 mul, 1 flop
      float vdy00 = dsv0 * eya0; // 1 mul, 1 flop
      float vdx01 = dsv1 * exa1; // 1 mul, 1 flop
      float vdy01 = dsv1 * eya1; // 1 mul, 1 flop
      float vdx02 = dsv2 * exa2; // 1 mul, 1 flop
      float vdy02 = dsv2 * eya2; // 1 mul, 1 flop
      float vdx03 = dsv3 * exa3; // 1 mul, 1 flop
      float vdy03 = dsv3 * eya3; // 1 mul, 1 flop
      float vdx04 = dsv4 * exa4; // 1 mul, 1 flop
      float vdy04 = dsv4 * eya4; // 1 mul, 1 flop
      float vdx05 = dsv5 * exa5; // 1 mul, 1 flop
      float vdy05 = dsv5 * eya5; // 1 mul, 1 flop
      float vdx06 = dsv6 * exa6; // 1 mul, 1 flop
      float vdy06 = dsv6 * eya6; // 1 mul, 1 flop
      float vdx07 = dsv7 * exa7; // 1 mul, 1 flop
      float vdy07 = dsv7 * eya7; // 1 mul, 1 flop

      float vdx10 = vdx00 - avx0; // 1 add, 1 flop
      float vdy10 = vdy00 - avy0; // 1 add, 1 flop
      float vdx11 = vdx01 - avx1; // 1 add, 1 flop
      float vdy11 = vdy01 - avy1; // 1 add, 1 flop
      float vdx12 = vdx02 - avx2; // 1 add, 1 flop
      float vdy12 = vdy02 - avy2; // 1 add, 1 flop
      float vdx13 = vdx03 - avx3; // 1 add, 1 flop
      float vdy13 = vdy03 - avy3; // 1 add, 1 flop
      float vdx14 = vdx04 - avx4; // 1 add, 1 flop
      float vdy14 = vdy04 - avy4; // 1 add, 1 flop
      float vdx15 = vdx05 - avx5; // 1 add, 1 flop
      float vdy15 = vdy05 - avy5; // 1 add, 1 flop
      float vdx16 = vdx06 - avx6; // 1 add, 1 flop
      float vdy16 = vdy06 - avy6; // 1 add, 1 flop
      float vdx17 = vdx07 - avx7; // 1 add, 1 flop
      float vdy17 = vdy07 - avy7; // 1 add, 1 flop

      // apply realxation time
      sfx0 += INV_RELAX_TIME * vdx10; // 1 mul => 1 flops
      sfy0 += INV_RELAX_TIME * vdy10; // 1 mul => 1 flops
      sfx1 += INV_RELAX_TIME * vdx11; // 1 mul => 1 flops
      sfy1 += INV_RELAX_TIME * vdy11; // 1 mul => 1 flops
      sfx2 += INV_RELAX_TIME * vdx12; // 1 mul => 1 flops
      sfy2 += INV_RELAX_TIME * vdy12; // 1 mul => 1 flops
      sfx3 += INV_RELAX_TIME * vdx13; // 1 mul => 1 flops
      sfy3 += INV_RELAX_TIME * vdy13; // 1 mul => 1 flops
      sfx4 += INV_RELAX_TIME * vdx14; // 1 mul => 1 flops
      sfy4 += INV_RELAX_TIME * vdy14; // 1 mul => 1 flops
      sfx5 += INV_RELAX_TIME * vdx15; // 1 mul => 1 flops
      sfy5 += INV_RELAX_TIME * vdy15; // 1 mul => 1 flops
      sfx6 += INV_RELAX_TIME * vdx16; // 1 mul => 1 flops
      sfy6 += INV_RELAX_TIME * vdy16; // 1 mul => 1 flops
      sfx7 += INV_RELAX_TIME * vdx17; // 1 mul => 1 flops
      sfy7 += INV_RELAX_TIME * vdy17; // 1 mul => 1 flops

      social_force[IndexX(i)] = sfx0;
      social_force[IndexY(i, n)] = sfy0;
      social_force[IndexX(i + 1)] = sfx1;
      social_force[IndexY(i + 1, n)] = sfy1;
      social_force[IndexX(i + 2)] = sfx2;
      social_force[IndexY(i + 2, n)] = sfy2;
      social_force[IndexX(i + 3)] = sfx3;
      social_force[IndexY(i + 3, n)] = sfy3;
      social_force[IndexX(i + 4)] = sfx4;
      social_force[IndexY(i + 4, n)] = sfy4;
      social_force[IndexX(i + 5)] = sfx5;
      social_force[IndexY(i + 5, n)] = sfy5;
      social_force[IndexX(i + 6)] = sfx6;
      social_force[IndexY(i + 6, n)] = sfy6;
      social_force[IndexX(i + 7)] = sfx7;
      social_force[IndexY(i + 7, n)] = sfy7;
    } //n * (12*(n-1) + 3*(n_borders) + 2) ADDS
    //n * (18*(n-1) + 8*(n_borders) + 4) MULTS
    //n * (6*(n-1) + n_borders )         DIVS
    //n * (n-1 + n_borders )             EXPS
    //n * (4*(n-1))                      SQRTS
    for (int i = 0; i < n - 3; i += 8)
    {
      /************************************************/
      // LOADS
      /************************************************/
      float cx0 = position[IndexX(i)];
      float cy0 = position[IndexY(i, n)];
      float cx1 = position[IndexX(i + 1)];
      float cy1 = position[IndexY(i + 1, n)];
      float cx2 = position[IndexX(i + 2)];
      float cy2 = position[IndexY(i + 2, n)];
      float cx3 = position[IndexX(i + 3)];
      float cy3 = position[IndexY(i + 3, n)];
      float cx4 = position[IndexX(i + 4)];
      float cy4 = position[IndexY(i + 4, n)];
      float cx5 = position[IndexX(i + 5)];
      float cy5 = position[IndexY(i + 5, n)];
      float cx6 = position[IndexX(i + 6)];
      float cy6 = position[IndexY(i + 6, n)];
      float cx7 = position[IndexX(i + 7)];
      float cy7 = position[IndexY(i + 7, n)];

      float max0 = desired_max_speed[i];
      float max1 = desired_max_speed[i + 1];
      float max2 = desired_max_speed[i + 2];
      float max3 = desired_max_speed[i + 3];
      float max4 = desired_max_speed[i + 4];
      float max5 = desired_max_speed[i + 5];
      float max6 = desired_max_speed[i + 6];
      float max7 = desired_max_speed[i + 7];

      //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
      float pvx0 = actual_velocity[IndexX(i)] + social_force[IndexX(i)] * TIMESTEP;               // 1 add, 1 mult => 2 flops
      float pvy0 = actual_velocity[IndexY(i, n)] + social_force[IndexY(i, n)] * TIMESTEP;         // 1 add, 1 mult => 2 flops
      float pvx1 = actual_velocity[IndexX(i + 1)] + social_force[IndexX(i + 1)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      float pvy1 = actual_velocity[IndexY(i + 1, n)] + social_force[IndexY(i + 1, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      float pvx2 = actual_velocity[IndexX(i + 2)] + social_force[IndexX(i + 2)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      float pvy2 = actual_velocity[IndexY(i + 2, n)] + social_force[IndexY(i + 2, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      float pvx3 = actual_velocity[IndexX(i + 3)] + social_force[IndexX(i + 3)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      float pvy3 = actual_velocity[IndexY(i + 3, n)] + social_force[IndexY(i + 3, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      float pvx4 = actual_velocity[IndexX(i + 4)] + social_force[IndexX(i + 4)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      float pvy4 = actual_velocity[IndexY(i + 4, n)] + social_force[IndexY(i + 4, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      float pvx5 = actual_velocity[IndexX(i + 5)] + social_force[IndexX(i + 5)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      float pvy5 = actual_velocity[IndexY(i + 5, n)] + social_force[IndexY(i + 5, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      float pvx6 = actual_velocity[IndexX(i + 6)] + social_force[IndexX(i + 6)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      float pvy6 = actual_velocity[IndexY(i + 6, n)] + social_force[IndexY(i + 6, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      float pvx7 = actual_velocity[IndexX(i + 7)] + social_force[IndexX(i + 7)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      float pvy7 = actual_velocity[IndexY(i + 7, n)] + social_force[IndexY(i + 7, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops

      /************************************************/
      // UPDATE POSITION
      /************************************************/
      //compute the norm of the preferd velocity
      float xysq00 = pvx0 * pvx0; // 1 mult => 1 flop
      float xysq10 = pvx1 * pvx1; // 1 mult => 1 flop
      float xysq20 = pvx2 * pvx2; // 1 mult => 1 flop
      float xysq30 = pvx3 * pvx3; // 1 mult => 1 flop
      float xysq40 = pvx4 * pvx4; // 1 mult => 1 flop
      float xysq50 = pvx5 * pvx5; // 1 mult => 1 flop
      float xysq60 = pvx6 * pvx6; // 1 mult => 1 flop
      float xysq70 = pvx7 * pvx7; // 1 mult => 1 flop

      float xysq01 = xysq00 + (pvy0 * pvy0); // 1 add, 1 mult1 => 2 flops
      float xysq11 = xysq10 + (pvy1 * pvy1); // 1 add, 1 mult1 => 2 flops
      float xysq21 = xysq20 + (pvy2 * pvy2); // 1 add, 1 mult1 => 2 flops
      float xysq31 = xysq30 + (pvy3 * pvy3); // 1 add, 1 mult1 => 2 flops
      float xysq41 = xysq40 + (pvy4 * pvy4); // 1 add, 1 mult1 => 2 flops
      float xysq51 = xysq50 + (pvy5 * pvy5); // 1 add, 1 mult1 => 2 flops
      float xysq61 = xysq60 + (pvy6 * pvy6); // 1 add, 1 mult1 => 2 flops
      float xysq71 = xysq70 + (pvy7 * pvy7); // 1 add, 1 mult1 => 2 flops

      float nv0 = sqrt(xysq01); // 1 sqrt => 1 flops
      float nv1 = sqrt(xysq11); // 1 sqrt => 1 flops
      float nv2 = sqrt(xysq21); // 1 sqrt => 1 flops
      float nv3 = sqrt(xysq31); // 1 sqrt => 1 flops
      float nv4 = sqrt(xysq41); // 1 sqrt => 1 flops
      float nv5 = sqrt(xysq51); // 1 sqrt => 1 flops
      float nv6 = sqrt(xysq61); // 1 sqrt => 1 flops
      float nv7 = sqrt(xysq71); // 1 sqrt => 1 flops

      //formula 12 in the paper --> compute control_value according to norm
      float cv0 = nv0 > max0 ? (max0 / nv0) : 1.0; // 1 div => 1 flops
      float cv1 = nv1 > max1 ? (max1 / nv1) : 1.0; // 1 div => 1 flops
      float cv2 = nv2 > max2 ? (max2 / nv2) : 1.0; // 1 div => 1 flops
      float cv3 = nv3 > max3 ? (max3 / nv3) : 1.0; // 1 div => 1 flops
      float cv4 = nv4 > max4 ? (max4 / nv4) : 1.0; // 1 div => 1 flops
      float cv5 = nv5 > max5 ? (max5 / nv5) : 1.0; // 1 div => 1 flops
      float cv6 = nv6 > max6 ? (max6 / nv6) : 1.0; // 1 div => 1 flops
      float cv7 = nv7 > max7 ? (max7 / nv7) : 1.0; // 1 div => 1 flops

      //apply control value
      pvx0 *= cv0; // 1 mul, 1 flop
      pvy0 *= cv0; // 1 mul, 1 flop
      pvx1 *= cv1; // 1 mul, 1 flop
      pvy1 *= cv1; // 1 mul, 1 flop
      pvx2 *= cv2; // 1 mul, 1 flop
      pvy2 *= cv2; // 1 mul, 1 flop
      pvx3 *= cv3; // 1 mul, 1 flop
      pvy3 *= cv3; // 1 mul, 1 flop
      pvx4 *= cv4; // 1 mul, 1 flop
      pvy4 *= cv4; // 1 mul, 1 flop
      pvx5 *= cv5; // 1 mul, 1 flop
      pvy5 *= cv5; // 1 mul, 1 flop
      pvx6 *= cv6; // 1 mul, 1 flop
      pvy6 *= cv6; // 1 mul, 1 flop
      pvx7 *= cv7; // 1 mul, 1 flop
      pvy7 *= cv7; // 1 mul, 1 flop

      cx0 += pvx0 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cy0 += pvy0 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cx1 += pvx1 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cy1 += pvy1 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cx2 += pvx2 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cy2 += pvy2 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cx3 += pvx3 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cy3 += pvy3 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cx4 += pvx4 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cy4 += pvy4 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cx5 += pvx5 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cy5 += pvy5 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cx6 += pvx6 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cy6 += pvy6 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cx7 += pvx7 * TIMESTEP; // 1 add, 1 mul => 2 flops
      cy7 += pvy7 * TIMESTEP; // 1 add, 1 mul => 2 flops

      /************************************************/
      //UPDATE DESIRED DIRECTION
      /************************************************/
      // get current position and target

      // compute differences
      float dx0 = final_destination[IndexX(i)] - cx0;        // 1 add => 1 flop
      float dy0 = final_destination[IndexY(i, n)] - cy0;     // 1 add => 1 flop
      float dx1 = final_destination[IndexX(i + 1)] - cx1;    // 1 add => 1 flop
      float dy1 = final_destination[IndexY(i + 1, n)] - cy1; // 1 add => 1 flop
      float dx2 = final_destination[IndexX(i + 2)] - cx2;    // 1 add => 1 flop
      float dy2 = final_destination[IndexY(i + 2, n)] - cy2; // 1 add => 1 flop
      float dx3 = final_destination[IndexX(i + 3)] - cx3;    // 1 add => 1 flop
      float dy3 = final_destination[IndexY(i + 3, n)] - cy3; // 1 add => 1 flop
      float dx4 = final_destination[IndexX(i + 4)] - cx4;    // 1 add => 1 flop
      float dy4 = final_destination[IndexY(i + 4, n)] - cy4; // 1 add => 1 flop
      float dx5 = final_destination[IndexX(i + 5)] - cx5;    // 1 add => 1 flop
      float dy5 = final_destination[IndexY(i + 5, n)] - cy5; // 1 add => 1 flop
      float dx6 = final_destination[IndexX(i + 6)] - cx6;    // 1 add => 1 flop
      float dy6 = final_destination[IndexY(i + 6, n)] - cy6; // 1 add => 1 flop
      float dx7 = final_destination[IndexX(i + 7)] - cx7;    // 1 add => 1 flop
      float dy7 = final_destination[IndexY(i + 7, n)] - cy7; // 1 add => 1 flop

      float d0_0 = dx0 * dx0; // 1 add, 2 mult => 3 flops
      float d1_0 = dx1 * dx1; // 1 add, 2 mult => 3 flops
      float d2_0 = dx2 * dx2; // 1 add, 2 mult => 3 flops
      float d3_0 = dx3 * dx3; // 1 add, 2 mult => 3 flops
      float d4_0 = dx4 * dx4; // 1 add, 2 mult => 3 flops
      float d5_0 = dx5 * dx5; // 1 add, 2 mult => 3 flops
      float d6_0 = dx6 * dx6; // 1 add, 2 mult => 3 flops
      float d7_0 = dx7 * dx7; // 1 add, 2 mult => 3 flops

      float d0_1 = d0_0 + dy0 * dy0; // 1 add, 2 mult => 3 flops
      float d1_1 = d1_0 + dy1 * dy1; // 1 add, 2 mult => 3 flops
      float d2_1 = d2_0 + dy2 * dy2; // 1 add, 2 mult => 3 flops
      float d3_1 = d3_0 + dy3 * dy3; // 1 add, 2 mult => 3 flops
      float d4_1 = d4_0 + dy4 * dy4; // 1 add, 2 mult => 3 flops
      float d5_1 = d5_0 + dy5 * dy5; // 1 add, 2 mult => 3 flops
      float d6_1 = d6_0 + dy6 * dy6; // 1 add, 2 mult => 3 flops
      float d7_1 = d7_0 + dy7 * dy7; // 1 add, 2 mult => 3 flops

      float n0 = sqrt(d0_1); // 1 sqrt => 1 flop
      float n1 = sqrt(d1_1); // 1 sqrt => 1 flop
      float n2 = sqrt(d2_1); // 1 sqrt => 1 flop
      float n3 = sqrt(d3_1); // 1 sqrt => 1 flop
      float n4 = sqrt(d4_1); // 1 sqrt => 1 flop
      float n5 = sqrt(d5_1); // 1 sqrt => 1 flop
      float n6 = sqrt(d6_1); // 1 sqrt => 1 flop
      float n7 = sqrt(d7_1); // 1 sqrt => 1 flop

      social_force[IndexX(i)] = 0.0;
      social_force[IndexY(i, n)] = 0.0;
      social_force[IndexX(i + 1)] = 0.0;
      social_force[IndexY(i + 1, n)] = 0.0;
      social_force[IndexX(i + 2)] = 0.0;
      social_force[IndexY(i + 2, n)] = 0.0;
      social_force[IndexX(i + 3)] = 0.0;
      social_force[IndexY(i + 3, n)] = 0.0;
      social_force[IndexX(i + 4)] = 0.0;
      social_force[IndexY(i + 4, n)] = 0.0;
      social_force[IndexX(i + 5)] = 0.0;
      social_force[IndexY(i + 5, n)] = 0.0;
      social_force[IndexX(i + 6)] = 0.0;
      social_force[IndexY(i + 6, n)] = 0.0;
      social_force[IndexX(i + 7)] = 0.0;
      social_force[IndexY(i + 7, n)] = 0.0;

      //update position
      position[IndexX(i)] = cx0;
      position[IndexY(i, n)] = cy0;
      position[IndexX(i + 1)] = cx1;
      position[IndexY(i + 1, n)] = cy1;
      position[IndexX(i + 2)] = cx2;
      position[IndexY(i + 2, n)] = cy2;
      position[IndexX(i + 3)] = cx3;
      position[IndexY(i + 3, n)] = cy3;
      position[IndexX(i + 4)] = cx4;
      position[IndexY(i + 4, n)] = cy4;
      position[IndexX(i + 5)] = cx5;
      position[IndexY(i + 5, n)] = cy5;
      position[IndexX(i + 6)] = cx6;
      position[IndexY(i + 6, n)] = cy6;
      position[IndexX(i + 7)] = cx7;
      position[IndexY(i + 7, n)] = cy7;

      //update speed value, desire direction, actual_velocity
      speed[i] = cv0 * nv0;     // 1 mul, 1 flop
      speed[i + 1] = cv1 * nv1; // 1 mul, 1 flop
      speed[i + 2] = cv2 * nv2; // 1 mul, 1 flop
      speed[i + 3] = cv3 * nv3; // 1 mul, 1 flop
      speed[i + 4] = cv4 * nv4; // 1 mul, 1 flop
      speed[i + 5] = cv5 * nv5; // 1 mul, 1 flop
      speed[i + 6] = cv6 * nv6; // 1 mul, 1 flop
      speed[i + 7] = cv7 * nv7; // 1 mul, 1 flop
      actual_velocity[IndexX(i)] = pvx0;
      actual_velocity[IndexY(i, n)] = pvy0;
      actual_velocity[IndexX(i + 1)] = pvx1;
      actual_velocity[IndexY(i + 1, n)] = pvy1;
      actual_velocity[IndexX(i + 2)] = pvx2;
      actual_velocity[IndexY(i + 2, n)] = pvy2;
      actual_velocity[IndexX(i + 3)] = pvx3;
      actual_velocity[IndexY(i + 3, n)] = pvy3;
      actual_velocity[IndexX(i + 4)] = pvx4;
      actual_velocity[IndexY(i + 4, n)] = pvy4;
      actual_velocity[IndexX(i + 5)] = pvx5;
      actual_velocity[IndexY(i + 5, n)] = pvy5;
      actual_velocity[IndexX(i + 6)] = pvx6;
      actual_velocity[IndexY(i + 6, n)] = pvy6;
      actual_velocity[IndexX(i + 7)] = pvx7;
      actual_velocity[IndexY(i + 7, n)] = pvy7;

      // update desired_direction
      desired_direction[IndexX(i)] = dx0 / n0;        // 1 div => 1 flop
      desired_direction[IndexY(i, n)] = dy0 / n0;     // 1 div => 1 flop
      desired_direction[IndexX(i + 1)] = dx1 / n1;    // 1 div => 1 flop
      desired_direction[IndexY(i + 1, n)] = dy1 / n1; // 1 div => 1 flop
      desired_direction[IndexX(i + 2)] = dx2 / n2;    // 1 div => 1 flop
      desired_direction[IndexY(i + 2, n)] = dy2 / n2; // 1 div => 1 flop
      desired_direction[IndexX(i + 3)] = dx3 / n3;    // 1 div => 1 flop
      desired_direction[IndexY(i + 3, n)] = dy3 / n3; // 1 div => 1 flop
      desired_direction[IndexX(i + 4)] = dx4 / n4;    // 1 div => 1 flop
      desired_direction[IndexY(i + 4, n)] = dy4 / n4; // 1 div => 1 flop
      desired_direction[IndexX(i + 5)] = dx5 / n5;    // 1 div => 1 flop
      desired_direction[IndexY(i + 5, n)] = dy5 / n5; // 1 div => 1 flop
      desired_direction[IndexX(i + 6)] = dx6 / n6;    // 1 div => 1 flop
      desired_direction[IndexY(i + 6, n)] = dy6 / n6; // 1 div => 1 flop
      desired_direction[IndexX(i + 7)] = dx7 / n7;    // 1 div => 1 flop
      desired_direction[IndexY(i + 7, n)] = dy7 / n7; // 1 div => 1 flop
    }                                                 // n * 8    ADDS
                                                      // n * 11   MULTS
                                                      // n * 3    DIVS
                                                      // n * 2    SQRTS
  }                                                   // n_timesteps * [n * (12*(n-1) + 3*(n_borders) + 2) + n * 8]         ADDS
                                                      // n_timesteps * [n * (18*(n-1) + 8*(n_borders) + 4) + n + n * 11]    MULTS
                                                      // n_timesteps * [n * (6*(n-1) + n_borders) + n * 3]                  DIVS
                                                      // n_timesteps * [n * (n-1 + n_borders)]                              EXPS
                                                      // n_timesteps * [n * (4*(n-1)) + (n * 2)]                            SQRTD
  //printf("fma cycles 1: %llu  -  sqrt cycles 1: %llu\n",fma_cycles_1, sqrt_cycles_1);
  // printf("fma cycles 2: %llu  -  sqrt cycles 2: %llu\n",fma_cycles_2, sqrt_cycles_2);
}