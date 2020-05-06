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

    {

      __m256 position_i_x;
      __m256 position_i_y;

      __m256 position_j_x;
      __m256 position_j_y;

      __m256 r_ab_x;
      __m256 r_ab_y;
      __m256 e_a_x;
      __m256 e_a_y;
      __m256 e_b_x;
      __m256 e_b_y;

      __m256 vb;
      __m256 delta_b;

      __m256 r_ab_2_x;
      __m256 r_ab_2_y;
      __m256 r_ab_norm;
      __m256 r_ab_norm_inv;

      __m256 r_ab_me_x;
      __m256 r_ab_me_y;

      __m256 r_ab_me_2_x;
      __m256 r_ab_me_2_y;
      __m256 r_ab_me_norm;
      __m256 r_ab_me_norm_inv;

      __m256 norm_sum;
      __m256 repulsion_x;
      __m256 repulsion_y;

      __m256 norm_sum_2, delta_b_2;
      __m256 b;
      __m256 exp;
      __m256 common_factor;
      __m256 check;
      __m256 check_x;
      __m256 check_y;

      __m256 repulsion_2_x;
      __m256 repulsion_2_y;
      __m256 repulsion_norm;
      __m256 threshold;
      __m256 w;
      __m256 mask;

      __m256 timestep_vec = _mm256_set1_ps(-TIMESTEP);
      __m256 minus_sigma_inv_vec = _mm256_set1_ps(-1.0 / SIGMA);
      __m256 div_factor_vec = _mm256_set1_ps(DIV_FACTOR);
      __m256 projection_factor_vec = _mm256_set1_ps(PROJECTION_FACTOR);
      __m256 influencer_vec = _mm256_set1_ps(INFLUENCE);

      __m256 one = _mm256_set1_ps(1);
      __m256 half_vec = _mm256_set1_ps(0.5);
      __m256 minus1_vec = _mm256_set1_ps(-1);
      __m256 eps = _mm256_set1_ps(1e-12);
      __m256 exp_constant = _mm256_set1_ps(0.00006103515); // 1 / 16384

      __m256 two_vec = _mm256_set1_ps(2);
      __m256 sigma_vec = _mm256_set1_ps(SIGMA);

      __m256 current_mask;

      // diagonal
      int n_eight = n / 8;
      for (int k = 0; k < n_eight; k++)
      {
        position_j_x = _mm256_load_ps(position + 8 * k);
        position_j_y = _mm256_load_ps(position + n + 8 * k);

        e_b_x = _mm256_load_ps(desired_direction + 8 * k);
        e_b_y = _mm256_load_ps(desired_direction + n + 8 * k);

        vb = _mm256_load_ps(speed + 8 * k);

        for (int i = 8 * k; i < 8 * (k + 1); i++)
        {
          position_i_x = _mm256_broadcast_ss(position + i);
          position_i_y = _mm256_broadcast_ss(position + n + i);
          e_a_x = _mm256_broadcast_ss(desired_direction + i);
          e_a_y = _mm256_broadcast_ss(desired_direction + n + i);

          r_ab_x = _mm256_sub_ps(position_i_x, position_j_x);
          r_ab_y = _mm256_sub_ps(position_i_y, position_j_y);

          delta_b = vb;

          // compute norm r_ab
          r_ab_2_x = _mm256_mul_ps(r_ab_x, r_ab_x);
          r_ab_2_y = _mm256_mul_ps(r_ab_y, r_ab_y);
          r_ab_norm = _mm256_sqrt_ps(_mm256_add_ps(r_ab_2_x, r_ab_2_y));

          r_ab_me_x = _mm256_sub_ps(r_ab_x, _mm256_mul_ps(delta_b, e_b_x));
          r_ab_me_y = _mm256_sub_ps(r_ab_y, _mm256_mul_ps(delta_b, e_b_y));

          // compute norm r_ab_me
          r_ab_me_2_x = _mm256_mul_ps(r_ab_me_x, r_ab_me_x);
          r_ab_me_2_y = _mm256_mul_ps(r_ab_me_y, r_ab_me_y);
          r_ab_me_norm = _mm256_sqrt_ps(_mm256_add_ps(r_ab_me_2_x, r_ab_me_2_y));

          norm_sum = _mm256_add_ps(r_ab_norm, r_ab_me_norm);

          repulsion_x = _mm256_div_ps(r_ab_x, r_ab_norm);
          repulsion_x = _mm256_add_ps(repulsion_x, _mm256_div_ps(r_ab_me_x, r_ab_me_norm));

          repulsion_y = _mm256_div_ps(r_ab_y, r_ab_norm);
          repulsion_y = _mm256_add_ps(repulsion_y, _mm256_div_ps(r_ab_me_y, r_ab_me_norm));

          norm_sum_2 = _mm256_mul_ps(norm_sum, norm_sum);
          delta_b_2 = _mm256_mul_ps(delta_b, delta_b);
          b = _mm256_sqrt_ps(_mm256_sub_ps(norm_sum_2, delta_b_2));
          b = _mm256_div_ps(b, two_vec);

          exp = _mm256_div_ps(b, sigma_vec);
          exp = _mm256_mul_ps(exp, minus1_vec);
          exp = exp_fast_vec_2_5_1(exp, one, exp_constant);

          common_factor = _mm256_mul_ps(norm_sum, div_factor_vec);
          common_factor = _mm256_div_ps(common_factor, b);
          common_factor = _mm256_mul_ps(exp, common_factor);

          repulsion_x = _mm256_mul_ps(repulsion_x, common_factor);
          repulsion_y = _mm256_mul_ps(repulsion_y, common_factor);

          // compute norm r_ab
          repulsion_2_x = _mm256_mul_ps(repulsion_x, repulsion_x);
          repulsion_2_y = _mm256_mul_ps(repulsion_y, repulsion_y);
          threshold = _mm256_sqrt_ps(_mm256_add_ps(repulsion_2_x, repulsion_2_y));

          check_x = _mm256_mul_ps(e_a_x, repulsion_x);
          check_y = _mm256_mul_ps(e_a_y, repulsion_y);

          check = _mm256_add_ps(check_x, check_y);

          threshold = _mm256_mul_ps(threshold, projection_factor_vec);

          mask = _mm256_cmp_ps(_mm256_mul_ps(check, minus1_vec), threshold, _CMP_GE_OQ);

          w = _mm256_blendv_ps(influencer_vec, one, mask);

          repulsion_x = _mm256_mul_ps(w, repulsion_x);
          repulsion_y = _mm256_mul_ps(w, repulsion_y);

          switch (i % 8)
          {
          case 0:
            repulsion_x = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_x, 0b11111110);
            repulsion_y = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_y, 0b11111110);

            break;
          case 1:
            repulsion_x = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_x, 0b11111101);
            repulsion_y = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_y, 0b11111101);
            break;
          case 2:
            repulsion_x = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_x, 0b11111011);
            repulsion_y = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_y, 0b11111011);
            break;
          case 3:
            repulsion_x = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_x, 0b11110111);
            repulsion_y = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_y, 0b11110111);
            break;
          case 4:
            repulsion_x = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_x, 0b11101111);
            repulsion_y = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_y, 0b11101111);
            break;
          case 5:
            repulsion_x = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_x, 0b11011111);
            repulsion_y = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_y, 0b11011111);
            break;
          case 6:
            repulsion_x = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_x, 0b10111111);
            repulsion_y = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_y, 0b10111111);
            break;
          case 7:
            repulsion_x = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_x, 0b01111111);
            repulsion_y = _mm256_blend_ps(_mm256_setzero_ps(), repulsion_y, 0b01111111);
            break;
          }

          repulsion_x = _mm256_hadd_ps(repulsion_x, _mm256_permute2f128_ps(repulsion_x, repulsion_x, 3));
          repulsion_x = _mm256_hadd_ps(repulsion_x, repulsion_x);
          repulsion_x = _mm256_hadd_ps(repulsion_x, repulsion_x);

          repulsion_y = _mm256_hadd_ps(repulsion_y, _mm256_permute2f128_ps(repulsion_y, repulsion_y, 3));
          repulsion_y = _mm256_hadd_ps(repulsion_y, repulsion_y);
          repulsion_y = _mm256_hadd_ps(repulsion_y, repulsion_y);

          social_force[IndexX(i)] += repulsion_x[0];
          social_force[IndexY(i, n)] += repulsion_y[0];
        }
      }
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