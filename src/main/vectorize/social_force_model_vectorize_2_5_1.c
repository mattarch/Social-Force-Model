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

      __m256 rxa = _mm256_load_ps(position + IndexX(i));
      __m256 rya = _mm256_load_ps(position + IndexY(i, n));
      __m256 one = _mm256_set1_ps(1.0);
      __m256 minus_one = _mm256_set1_ps(-1.0);
      __m256 two = _mm256_set1_ps(2.0);
      __m256 exp_constant = _mm256_set1_ps(0.00006103515); // 1 / 16384
      __m256 inv_sigma_vec = _mm256_set1_ps(-inv_sigma);
      __m256 div_factor_vec = _mm256_set1_ps(DIV_FACTOR);
      __m256 cospsi_vec = _mm256_set1_ps(cospsi);
      __m256 influencer_vec = _mm256_set1_ps(INFLUENCE);
      __m256 utimesr_vec = _mm256_set1_ps(UTIMESR);
      __m256 invr_vec = _mm256_set1_ps(1 / R);

      __m256 da = _mm256_load_ps(speed + IndexX(i));
      __m256 exa = _mm256_load_ps(desired_direction + IndexX(i));
      __m256 eya = _mm256_load_ps(desired_direction + IndexY(i, n));
      __m256 sfx = _mm256_setzero_ps();
      __m256 sfy = _mm256_setzero_ps();
      __m256 fb_vec = _mm256_set1_ps(fb);
      __m256 sb_vec = _mm256_set1_ps(sb);

      /************************************************/
      // UPDATE BORDER REPULSION TERM
      /************************************************/
      __m256 ryaB0 = _mm256_sub_ps(rya, fb_vec);
      __m256 se0 = exp_fast_vec_2_5_1(_mm256_mul_ps(ryaB0, invr_vec), one, exp_constant);
      se0 = _mm256_mul_ps(se0, utimesr_vec);
      se0 = _mm256_div_ps(se0, ryaB0);
      se0 = _mm256_mul_ps(se0, minus_one);
      __m256 rb0 = _mm256_mul_ps(se0, ryaB0);

      __m256 ryaB1 = _mm256_sub_ps(rya, sb_vec);
      __m256 se1 = exp_fast_vec_2_5_1(_mm256_mul_ps(_mm256_mul_ps(ryaB1,minus_one), invr_vec), one, exp_constant);
      se1 = _mm256_mul_ps(se1, utimesr_vec);
      se1 = _mm256_div_ps(se1, ryaB1);
      __m256 rb1 = _mm256_mul_ps(se1, ryaB1);

      sfy = _mm256_add_ps(rb0, rb1);

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
        rabnorm = _mm256_rsqrt_ps(rabnorm);

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
        rabmenorm = _mm256_rsqrt_ps(rabmenorm);

        __m256 rxbame_2 = _mm256_mul_ps(rxbame, rxbame);
        __m256 rbamenorm = _mm256_fmadd_ps(rybame, rybame, rxbame_2);
        rbamenorm = _mm256_rsqrt_ps(rbamenorm);

        __m256 normsumab = _mm256_add_ps(_mm256_rcp_ps(rabnorm), _mm256_rcp_ps(rabmenorm));
        __m256 normsumba = _mm256_add_ps(_mm256_rcp_ps(rabnorm), _mm256_rcp_ps(rbamenorm)); // maybe wrong use rbanorm

        // printf("rabmenorm %f %f\n", rabmenorm0, rabmenorm[0]);
        // printf("rbamenorm %f %f\n", rbamenorm0, rbamenorm[0]);

        __m256 div_1x = _mm256_mul_ps(rxab, rabnorm);
        __m256 div_1y = _mm256_mul_ps(ryab, rabnorm);
        __m256 div_2x = _mm256_mul_ps(rxabme, rabmenorm);
        __m256 div_2y = _mm256_mul_ps(ryabme, rabmenorm);

        __m256 repxab = _mm256_add_ps(div_1x, div_2x);
        __m256 repyab = _mm256_add_ps(div_1y, div_2y);

        __m256 div_2x_me = _mm256_mul_ps(rxbame, rbamenorm);
        __m256 div_2y_me = _mm256_mul_ps(rybame, rbamenorm);

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

        __m256 babnorm = _mm256_rsqrt_ps(bab);
        __m256 bbanorm = _mm256_rsqrt_ps(bba);
        babnorm = _mm256_mul_ps(babnorm, two);
        bbanorm = _mm256_mul_ps(bbanorm, two);

        __m256 cfab = exp_fast_vec_2_5_1(_mm256_mul_ps(_mm256_rcp_ps(babnorm), inv_sigma_vec), one, exp_constant);
        cfab = _mm256_mul_ps(cfab, div_factor_vec);
        cfab = _mm256_mul_ps(cfab, normsumab);
        cfab = _mm256_mul_ps(cfab, babnorm);

        __m256 cfba = exp_fast_vec_2_5_1(_mm256_mul_ps(_mm256_rcp_ps(bbanorm), inv_sigma_vec), one, exp_constant);
        cfba = _mm256_mul_ps(cfba, div_factor_vec);
        cfba = _mm256_mul_ps(cfba, normsumba);
        cfba = _mm256_mul_ps(cfba, bbanorm);

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
        __m256 thab = _mm256_rsqrt_ps(_mm256_fmadd_ps(repyab, repyab, repxab_2));
        thab = _mm256_mul_ps(_mm256_rcp_ps(thab), cospsi_vec);

        __m256 repxba_2 = _mm256_mul_ps(repxba, repxba);
        __m256 thba = _mm256_rsqrt_ps(_mm256_fmadd_ps(repyba, repyba, repxba_2));
        thba = _mm256_mul_ps(_mm256_rcp_ps(thba), cospsi_vec);

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

      /************************************************/
      //UPDATE ACCELERATION TERM
      /************************************************/
      // get actual velocity, desired direction, desired speed

      // compute velocity difference

      __m256 vdx;
      __m256 vdy;
      __m256 dsv = _mm256_load_ps(desired_speed + IndexX(i));
      __m256 avx = _mm256_load_ps(actual_velocity + IndexX(i));
      __m256 avy = _mm256_load_ps(actual_velocity + IndexY(i, n));
      __m256 inv_relax_time_vec = _mm256_set1_ps(INV_RELAX_TIME);

      vdx = _mm256_mul_ps(dsv, exa);
      vdy = _mm256_mul_ps(dsv, eya);

      vdx = _mm256_sub_ps(vdx, avx);
      vdy = _mm256_sub_ps(vdy, avy);

      sfx = _mm256_fmadd_ps(inv_relax_time_vec, vdx, sfx);
      sfy = _mm256_fmadd_ps(inv_relax_time_vec, vdy, sfy);

      __m256 timestep_vec = _mm256_set1_ps(TIMESTEP);

      /************************************************/
      // LOADS
      /************************************************/

      __m256 cx = _mm256_load_ps(position + IndexX(i));
      __m256 cy = _mm256_load_ps(position + IndexY(i, n));
      __m256 max = _mm256_load_ps(desired_max_speed + IndexX(i));

      __m256 social_force_x = _mm256_load_ps(social_force + IndexX(i));
      __m256 social_force_y = _mm256_load_ps(social_force + IndexY(i, n));
      social_force_x = _mm256_add_ps(social_force_x, sfx);
      social_force_y = _mm256_add_ps(social_force_y, sfy);

      __m256 pvx = _mm256_load_ps(actual_velocity + IndexX(i));
      __m256 pvy = _mm256_load_ps(actual_velocity + IndexY(i, n));

      //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
      pvx = _mm256_fmadd_ps(social_force_x, timestep_vec, pvx);
      pvy = _mm256_fmadd_ps(social_force_y, timestep_vec, pvy);

      /************************************************/
      // UPDATE POSITION
      /************************************************/
      //compute the norm of the preferd velocity

      __m256 pvx_2 = _mm256_mul_ps(pvx, pvx);
      __m256 nv_2 = _mm256_fmadd_ps(pvy, pvy, pvx_2);
      __m256 nv = _mm256_sqrt_ps(nv_2);

      __m256 mask = _mm256_cmp_ps(nv, max, _CMP_GT_OQ);
      __m256 cv = _mm256_blendv_ps(one, _mm256_div_ps(max, nv), mask);

      pvx = _mm256_mul_ps(pvx, cv);
      pvy = _mm256_mul_ps(pvy, cv);

      cx = _mm256_fmadd_ps(pvx, timestep_vec, cx);
      cy = _mm256_fmadd_ps(pvy, timestep_vec, cy);

      /************************************************/
      //UPDATE DESIRED DIRECTION
      /************************************************/
      // get current position and target
      __m256 final_destination_x = _mm256_load_ps(final_destination + IndexX(i));
      __m256 final_destination_y = _mm256_load_ps(final_destination + IndexY(i, n));

      // compute differences

      __m256 dx = _mm256_sub_ps(final_destination_x, cx);
      __m256 dy = _mm256_sub_ps(final_destination_y, cy);

      __m256 dx_2 = _mm256_mul_ps(dx, dx);
      __m256 d_2 = _mm256_fmadd_ps(dy, dy, dx_2);
      __m256 d = _mm256_sqrt_ps(d_2);

      _mm256_store_ps(social_force + IndexX(i), _mm256_setzero_ps());
      _mm256_store_ps(social_force + IndexY(i, n), _mm256_setzero_ps());
      _mm256_store_ps(position + IndexX(i), cx);
      _mm256_store_ps(position + IndexY(i, n), cy);
      _mm256_store_ps(actual_velocity + IndexX(i), pvx);
      _mm256_store_ps(actual_velocity + IndexY(i, n), pvy);
      _mm256_store_ps(desired_direction + IndexX(i), _mm256_div_ps(dx, d));
      _mm256_store_ps(desired_direction + IndexY(i, n), _mm256_div_ps(dy, d));
      _mm256_store_ps(speed + IndexX(i), _mm256_mul_ps(cv, nv));
    } // n * 8    ADDS
      // n * 11   MULTS
      // n * 3    DIVS
      // n * 2    SQRTS
  }   // n_timesteps * [n * (12*(n-1) + 3*(n_borders) + 2) + n * 8]         ADDS
  // n_timesteps * [n * (18*(n-1) + 8*(n_borders) + 4) + n + n * 11]    MULTS
  // n_timesteps * [n * (6*(n-1) + n_borders) + n * 3]                  DIVS
  // n_timesteps * [n * (n-1 + n_borders)]                              EXPS
  // n_timesteps * [n * (4*(n-1)) + (n * 2)]                            SQRTD
  //printf("fma cycles 1: %llu  -  sqrt cycles 1: %llu\n",fma_cycles_1, sqrt_cycles_1);
  // printf("fma cycles 2: %llu  -  sqrt cycles 2: %llu\n",fma_cycles_2, sqrt_cycles_2);
}