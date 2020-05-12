/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <immintrin.h>

#include "social_force_model_vectorize_4.h"
#include "../social_force.h"
#include "../utility.h"

extern char filename_global[80];

__m256 exp_fast_vec_4(__m256 x, __m256 one, __m256 exp_constant)
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


void simulation_basic_vectorize_4(int number_of_people, int n_timesteps, float *position, float *speed, float *desired_direction, float *final_destination, float *borders, float *actual_velocity, float *acceleration_term,
                                  float *people_repulsion_term, float *border_repulsion_term, float *social_force, float *desired_speed, float *desired_max_speed)
{
  // start simulation
  CONSOLE_PRINT(("Start simulation with %d persons\n", number_of_people));

  // simulate steps
  for (int step = 0; step < n_timesteps; step++)
  {
    {
      __m256 current_x;
      __m256 current_y;
      __m256 target_x;
      __m256 target_y;
      __m256 delta_x;
      __m256 delta_y;
      __m256 delta_x_2;
      __m256 delta_y_2;
      __m256 normalizer;

      __m256 result_x;
      __m256 result_y;

      __m256 one = _mm256_set1_ps(1);

      // iterate over all persons and update desired_direction
      for (int i = 0; i < number_of_people - 7; i += 8)
      {
        // get current position and target
        current_x = _mm256_load_ps(position + i); // now xy positions for two persons in register
        target_x = _mm256_load_ps(final_destination + i);
        current_y = _mm256_load_ps(position + number_of_people + i); // now xy positions for two persons in register
        target_y = _mm256_load_ps(final_destination + number_of_people + i);

        // compute differences
        delta_x = _mm256_sub_ps(target_x, current_x);
        delta_y = _mm256_sub_ps(target_y, current_y);

        // compute norm
        delta_x_2 = _mm256_mul_ps(delta_x, delta_x); // square each entry
        normalizer = _mm256_rsqrt_ps(_mm256_fmadd_ps(delta_y, delta_y, delta_x_2));

        result_x = _mm256_mul_ps(delta_x, normalizer);
        result_y = _mm256_mul_ps(delta_y, normalizer);

        _mm256_store_ps(desired_direction + i, result_x);
        _mm256_store_ps(desired_direction + number_of_people + i, result_y);
      }
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

      __m256 current_mask;

      // diagonal
      int n_eight = number_of_people / 8;
      for (int k = 0; k < n_eight; k++)
      {
        position_j_x = _mm256_load_ps(position + 8 * k);
        position_j_y = _mm256_load_ps(position + number_of_people + 8 * k);

        e_b_x = _mm256_load_ps(desired_direction + 8 * k);
        e_b_y = _mm256_load_ps(desired_direction + number_of_people + 8 * k);

        vb = _mm256_load_ps(speed + 8 * k);

        for (int i = 8 * k; i < 8 * (k + 1); i++)
        {
          position_i_x = _mm256_broadcast_ss(position + i);
          position_i_y = _mm256_broadcast_ss(position + number_of_people + i);
          e_a_x = _mm256_broadcast_ss(desired_direction + i);
          e_a_y = _mm256_broadcast_ss(desired_direction + number_of_people + i);

          r_ab_x = _mm256_sub_ps(position_i_x, position_j_x);
          r_ab_y = _mm256_sub_ps(position_i_y, position_j_y);

          delta_b = _mm256_mul_ps(vb, timestep_vec);

          // compute norm r_ab
          r_ab_2_x = _mm256_mul_ps(r_ab_x, r_ab_x);
          r_ab_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_y, r_ab_y, r_ab_2_x));

          r_ab_me_x = _mm256_fmadd_ps(delta_b, e_b_x, r_ab_x);
          r_ab_me_y = _mm256_fmadd_ps(delta_b, e_b_y, r_ab_y);

          // compute norm r_ab_me
          r_ab_me_2_x = _mm256_mul_ps(r_ab_me_x, r_ab_me_x);
          r_ab_me_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

          // sum up norms
          norm_sum = _mm256_add_ps(_mm256_rcp_ps(r_ab_norm), _mm256_rcp_ps(r_ab_me_norm));

          repulsion_x = _mm256_mul_ps(r_ab_x, r_ab_norm);
          repulsion_x = _mm256_fmadd_ps(r_ab_me_x, r_ab_me_norm, repulsion_x);

          repulsion_y = _mm256_mul_ps(r_ab_y, r_ab_norm);
          repulsion_y = _mm256_fmadd_ps(r_ab_me_y, r_ab_me_norm, repulsion_y);

          delta_b_2 = _mm256_mul_ps(delta_b, delta_b);
          b = _mm256_rsqrt_ps(_mm256_fmsub_ps(norm_sum, norm_sum, delta_b_2));
          b = _mm256_rcp_ps(b);
          b = _mm256_mul_ps(b, half_vec);

          exp = _mm256_mul_ps(b, minus_sigma_inv_vec);
          exp = exp_fast_vec_4(exp, one, exp_constant);

          common_factor = _mm256_mul_ps(norm_sum, div_factor_vec);
          common_factor = _mm256_div_ps(common_factor, b);
          common_factor = _mm256_mul_ps(exp, common_factor);

          repulsion_x = _mm256_mul_ps(repulsion_x, common_factor);
          repulsion_y = _mm256_mul_ps(repulsion_y, common_factor);

          // compute norm r_ab
          repulsion_2_y = _mm256_mul_ps(repulsion_y, repulsion_y);
          threshold = _mm256_rsqrt_ps(_mm256_fmadd_ps(repulsion_x, repulsion_x, repulsion_2_y));
          threshold = _mm256_rcp_ps(threshold);

          check_y = _mm256_mul_ps(e_a_y, repulsion_y);
          check = _mm256_fmadd_ps(e_a_x, repulsion_x, check_y);

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

          _mm256_store_ps(people_repulsion_term + 8 * i, repulsion_x);
          _mm256_store_ps(people_repulsion_term + 8 * number_of_people + 8 * i, repulsion_y);
        }
      }
      // first block row

      for (int k = 1; k < 2; k++)
      {
        for (int i = (k - 1) * 8; i < k * 8; i++)
        {
          position_i_x = _mm256_broadcast_ss(position + i);
          position_i_y = _mm256_broadcast_ss(position + number_of_people + i);
          e_a_x = _mm256_broadcast_ss(desired_direction + i);
          e_a_y = _mm256_broadcast_ss(desired_direction + number_of_people + i);

          __m256 people_repulsion_x = _mm256_load_ps(people_repulsion_term + 8 * i);
          __m256 people_repulsion_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * i);

          for (int j = k * 8; j < number_of_people - 7; j += 8)
          {
            position_j_x = _mm256_load_ps(position + j);
            position_j_y = _mm256_load_ps(position + number_of_people + j);

            r_ab_x = _mm256_sub_ps(position_i_x, position_j_x);
            r_ab_y = _mm256_sub_ps(position_i_y, position_j_y);

            e_b_x = _mm256_load_ps(desired_direction + j);
            e_b_y = _mm256_load_ps(desired_direction + number_of_people + j);

            vb = _mm256_load_ps(speed + j);

            delta_b = _mm256_mul_ps(vb, timestep_vec);

            // compute norm r_ab
            r_ab_2_x = _mm256_mul_ps(r_ab_x, r_ab_x);
            r_ab_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_y, r_ab_y, r_ab_2_x));

            r_ab_me_x = _mm256_fmadd_ps(delta_b, e_b_x, r_ab_x);
            r_ab_me_y = _mm256_fmadd_ps(delta_b, e_b_y, r_ab_y);

            // compute norm r_ab_me
            r_ab_me_2_x = _mm256_mul_ps(r_ab_me_x, r_ab_me_x);
            r_ab_me_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

            // sum up norms
            norm_sum = _mm256_add_ps(_mm256_rcp_ps(r_ab_norm), _mm256_rcp_ps(r_ab_me_norm));

            repulsion_x = _mm256_mul_ps(r_ab_x, r_ab_norm);
            repulsion_x = _mm256_fmadd_ps(r_ab_me_x, r_ab_me_norm, repulsion_x);

            repulsion_y = _mm256_mul_ps(r_ab_y, r_ab_norm);
            repulsion_y = _mm256_fmadd_ps(r_ab_me_y, r_ab_me_norm, repulsion_y);

            delta_b_2 = _mm256_mul_ps(delta_b, delta_b);
            b = _mm256_rsqrt_ps(_mm256_fmsub_ps(norm_sum, norm_sum, delta_b_2));
            b = _mm256_rcp_ps(b);
            b = _mm256_mul_ps(b, half_vec);

            exp = _mm256_mul_ps(b, minus_sigma_inv_vec);
            exp = exp_fast_vec_4(exp, one, exp_constant);

            common_factor = _mm256_mul_ps(norm_sum, div_factor_vec);
            common_factor = _mm256_div_ps(common_factor, b);
            common_factor = _mm256_mul_ps(exp, common_factor);

            repulsion_x = _mm256_mul_ps(repulsion_x, common_factor);
            repulsion_y = _mm256_mul_ps(repulsion_y, common_factor);

            // compute norm r_ab
            repulsion_2_y = _mm256_mul_ps(repulsion_y, repulsion_y);
            threshold = _mm256_rsqrt_ps(_mm256_fmadd_ps(repulsion_x, repulsion_x, repulsion_2_y));
            threshold = _mm256_rcp_ps(threshold);

            check_y = _mm256_mul_ps(e_a_y, repulsion_y);
            check = _mm256_fmadd_ps(e_a_x, repulsion_x, check_y);

            threshold = _mm256_mul_ps(threshold, projection_factor_vec);

            mask = _mm256_cmp_ps(_mm256_mul_ps(check, minus1_vec), threshold, _CMP_GE_OQ);

            w = _mm256_blendv_ps(influencer_vec, one, mask);

            repulsion_x = _mm256_mul_ps(w, repulsion_x);
            repulsion_y = _mm256_mul_ps(w, repulsion_y);

            people_repulsion_x = _mm256_add_ps(people_repulsion_x, repulsion_x);
            people_repulsion_y = _mm256_add_ps(people_repulsion_y, repulsion_y);
          }
          _mm256_store_ps(people_repulsion_term + 8 * i, people_repulsion_x);
          _mm256_store_ps(people_repulsion_term + 8 * number_of_people + 8 * i, people_repulsion_y);
        }
      }

      // upper triangle

      for (int k = 2; k < n_eight; k++)
      {
        for (int i = (k - 1) * 8; i < k * 8; i++)
        {
          position_i_x = _mm256_broadcast_ss(position + i);
          position_i_y = _mm256_broadcast_ss(position + number_of_people + i);
          e_a_x = _mm256_broadcast_ss(desired_direction + i);
          e_a_y = _mm256_broadcast_ss(desired_direction + number_of_people + i);

          __m256 people_repulsion_x = _mm256_load_ps(people_repulsion_term + 8 * i);
          __m256 people_repulsion_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * i);

          for (int j = 0; j < (k - 1) * 8; j += 8)
          {
            position_j_x = _mm256_load_ps(position + j);
            position_j_y = _mm256_load_ps(position + number_of_people + j);

            r_ab_x = _mm256_sub_ps(position_i_x, position_j_x);
            r_ab_y = _mm256_sub_ps(position_i_y, position_j_y);

            e_b_x = _mm256_load_ps(desired_direction + j);
            e_b_y = _mm256_load_ps(desired_direction + number_of_people + j);

            vb = _mm256_load_ps(speed + j);

            delta_b = _mm256_mul_ps(vb, timestep_vec);

            // compute norm r_ab
            r_ab_2_x = _mm256_mul_ps(r_ab_x, r_ab_x);
            r_ab_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_y, r_ab_y, r_ab_2_x));

            r_ab_me_x = _mm256_fmadd_ps(delta_b, e_b_x, r_ab_x);
            r_ab_me_y = _mm256_fmadd_ps(delta_b, e_b_y, r_ab_y);

            // compute norm r_ab_me
            r_ab_me_2_x = _mm256_mul_ps(r_ab_me_x, r_ab_me_x);
            r_ab_me_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

            // sum up norms
            norm_sum = _mm256_add_ps(_mm256_rcp_ps(r_ab_norm), _mm256_rcp_ps(r_ab_me_norm));

            repulsion_x = _mm256_mul_ps(r_ab_x, r_ab_norm);
            repulsion_x = _mm256_fmadd_ps(r_ab_me_x, r_ab_me_norm, repulsion_x);

            repulsion_y = _mm256_mul_ps(r_ab_y, r_ab_norm);
            repulsion_y = _mm256_fmadd_ps(r_ab_me_y, r_ab_me_norm, repulsion_y);

            delta_b_2 = _mm256_mul_ps(delta_b, delta_b);
            b = _mm256_rsqrt_ps(_mm256_fmsub_ps(norm_sum, norm_sum, delta_b_2));
            b = _mm256_rcp_ps(b);
            b = _mm256_mul_ps(b, half_vec);

            exp = _mm256_mul_ps(b, minus_sigma_inv_vec);
            exp = exp_fast_vec_4(exp, one, exp_constant);

            common_factor = _mm256_mul_ps(norm_sum, div_factor_vec);
            common_factor = _mm256_div_ps(common_factor, b);
            common_factor = _mm256_mul_ps(exp, common_factor);

            repulsion_x = _mm256_mul_ps(repulsion_x, common_factor);
            repulsion_y = _mm256_mul_ps(repulsion_y, common_factor);

            // compute norm r_ab
            repulsion_2_y = _mm256_mul_ps(repulsion_y, repulsion_y);
            threshold = _mm256_rsqrt_ps(_mm256_fmadd_ps(repulsion_x, repulsion_x, repulsion_2_y));
            threshold = _mm256_rcp_ps(threshold);

            check_y = _mm256_mul_ps(e_a_y, repulsion_y);
            check = _mm256_fmadd_ps(e_a_x, repulsion_x, check_y);

            threshold = _mm256_mul_ps(threshold, projection_factor_vec);

            mask = _mm256_cmp_ps(_mm256_mul_ps(check, minus1_vec), threshold, _CMP_GE_OQ);

            w = _mm256_blendv_ps(influencer_vec, one, mask);

            repulsion_x = _mm256_mul_ps(w, repulsion_x);
            repulsion_y = _mm256_mul_ps(w, repulsion_y);

            people_repulsion_x = _mm256_add_ps(people_repulsion_x, repulsion_x);
            people_repulsion_y = _mm256_add_ps(people_repulsion_y, repulsion_y);
          }

          for (int j = k * 8; j < number_of_people - 7; j += 8)
          {
            position_j_x = _mm256_load_ps(position + j);
            position_j_y = _mm256_load_ps(position + number_of_people + j);

            r_ab_x = _mm256_sub_ps(position_i_x, position_j_x);
            r_ab_y = _mm256_sub_ps(position_i_y, position_j_y);

            e_b_x = _mm256_load_ps(desired_direction + j);
            e_b_y = _mm256_load_ps(desired_direction + number_of_people + j);

            vb = _mm256_load_ps(speed + j);

            delta_b = _mm256_mul_ps(vb, timestep_vec);

            // compute norm r_ab
            r_ab_2_x = _mm256_mul_ps(r_ab_x, r_ab_x);
            r_ab_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_y, r_ab_y, r_ab_2_x));

            r_ab_me_x = _mm256_fmadd_ps(delta_b, e_b_x, r_ab_x);
            r_ab_me_y = _mm256_fmadd_ps(delta_b, e_b_y, r_ab_y);

            // compute norm r_ab_me
            r_ab_me_2_x = _mm256_mul_ps(r_ab_me_x, r_ab_me_x);
            r_ab_me_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

            // sum up norms
            norm_sum = _mm256_add_ps(_mm256_rcp_ps(r_ab_norm), _mm256_rcp_ps(r_ab_me_norm));

            repulsion_x = _mm256_mul_ps(r_ab_x, r_ab_norm);
            repulsion_x = _mm256_fmadd_ps(r_ab_me_x, r_ab_me_norm, repulsion_x);

            repulsion_y = _mm256_mul_ps(r_ab_y, r_ab_norm);
            repulsion_y = _mm256_fmadd_ps(r_ab_me_y, r_ab_me_norm, repulsion_y);

            delta_b_2 = _mm256_mul_ps(delta_b, delta_b);
            b = _mm256_rsqrt_ps(_mm256_fmsub_ps(norm_sum, norm_sum, delta_b_2));
            b = _mm256_rcp_ps(b);
            b = _mm256_mul_ps(b, half_vec);

            exp = _mm256_mul_ps(b, minus_sigma_inv_vec);
            exp = exp_fast_vec_4(exp, one, exp_constant);

            common_factor = _mm256_mul_ps(norm_sum, div_factor_vec);
            common_factor = _mm256_div_ps(common_factor, b);
            common_factor = _mm256_mul_ps(exp, common_factor);

            repulsion_x = _mm256_mul_ps(repulsion_x, common_factor);
            repulsion_y = _mm256_mul_ps(repulsion_y, common_factor);

            // compute norm r_ab
            repulsion_2_y = _mm256_mul_ps(repulsion_y, repulsion_y);
            threshold = _mm256_rsqrt_ps(_mm256_fmadd_ps(repulsion_x, repulsion_x, repulsion_2_y));
            threshold = _mm256_rcp_ps(threshold);

            check_y = _mm256_mul_ps(e_a_y, repulsion_y);
            check = _mm256_fmadd_ps(e_a_x, repulsion_x, check_y);

            threshold = _mm256_mul_ps(threshold, projection_factor_vec);

            mask = _mm256_cmp_ps(_mm256_mul_ps(check, minus1_vec), threshold, _CMP_GE_OQ);

            w = _mm256_blendv_ps(influencer_vec, one, mask);

            repulsion_x = _mm256_mul_ps(w, repulsion_x);
            repulsion_y = _mm256_mul_ps(w, repulsion_y);

            people_repulsion_x = _mm256_add_ps(people_repulsion_x, repulsion_x);
            people_repulsion_y = _mm256_add_ps(people_repulsion_y, repulsion_y);
          }

          _mm256_store_ps(people_repulsion_term + 8 * i, people_repulsion_x);
          _mm256_store_ps(people_repulsion_term + 8 * number_of_people + 8 * i, people_repulsion_y);
        }
      }
      // last row of block
      for (int k = n_eight - 1; k < n_eight; k++)
      {
        for (int i = k * 8; i < (k + 1) * 8; i++)
        {
          position_i_x = _mm256_broadcast_ss(position + i);
          position_i_y = _mm256_broadcast_ss(position + number_of_people + i);
          e_a_x = _mm256_broadcast_ss(desired_direction + i);
          e_a_y = _mm256_broadcast_ss(desired_direction + number_of_people + i);

          __m256 people_repulsion_x = _mm256_load_ps(people_repulsion_term + 8 * i);
          __m256 people_repulsion_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * i);

          for (int j = 0; j < k * 8; j += 8)
          {
            position_j_x = _mm256_load_ps(position + j);
            position_j_y = _mm256_load_ps(position + number_of_people + j);

            r_ab_x = _mm256_sub_ps(position_i_x, position_j_x);
            r_ab_y = _mm256_sub_ps(position_i_y, position_j_y);

            e_b_x = _mm256_load_ps(desired_direction + j);
            e_b_y = _mm256_load_ps(desired_direction + number_of_people + j);

            vb = _mm256_load_ps(speed + j);

            delta_b = _mm256_mul_ps(vb, timestep_vec);

            // compute norm r_ab
            r_ab_2_x = _mm256_mul_ps(r_ab_x, r_ab_x);
            r_ab_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_y, r_ab_y, r_ab_2_x));

            r_ab_me_x = _mm256_fmadd_ps(delta_b, e_b_x, r_ab_x);
            r_ab_me_y = _mm256_fmadd_ps(delta_b, e_b_y, r_ab_y);

            // compute norm r_ab_me
            r_ab_me_2_x = _mm256_mul_ps(r_ab_me_x, r_ab_me_x);
            r_ab_me_norm = _mm256_rsqrt_ps(_mm256_fmadd_ps(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

            // sum up norms
            norm_sum = _mm256_add_ps(_mm256_rcp_ps(r_ab_norm), _mm256_rcp_ps(r_ab_me_norm));

            repulsion_x = _mm256_mul_ps(r_ab_x, r_ab_norm);
            repulsion_x = _mm256_fmadd_ps(r_ab_me_x, r_ab_me_norm, repulsion_x);

            repulsion_y = _mm256_mul_ps(r_ab_y, r_ab_norm);
            repulsion_y = _mm256_fmadd_ps(r_ab_me_y, r_ab_me_norm, repulsion_y);

            delta_b_2 = _mm256_mul_ps(delta_b, delta_b);
            b = _mm256_rsqrt_ps(_mm256_fmsub_ps(norm_sum, norm_sum, delta_b_2));
            b = _mm256_rcp_ps(b);
            b = _mm256_mul_ps(b, half_vec);

            exp = _mm256_mul_ps(b, minus_sigma_inv_vec);
            exp = exp_fast_vec_4(exp, one, exp_constant);

            common_factor = _mm256_mul_ps(norm_sum, div_factor_vec);
            common_factor = _mm256_div_ps(common_factor, b);
            common_factor = _mm256_mul_ps(exp, common_factor);

            repulsion_x = _mm256_mul_ps(repulsion_x, common_factor);
            repulsion_y = _mm256_mul_ps(repulsion_y, common_factor);

            // compute norm r_ab
            repulsion_2_y = _mm256_mul_ps(repulsion_y, repulsion_y);
            threshold = _mm256_rsqrt_ps(_mm256_fmadd_ps(repulsion_x, repulsion_x, repulsion_2_y));
            threshold = _mm256_rcp_ps(threshold);

            check_y = _mm256_mul_ps(e_a_y, repulsion_y);
            check = _mm256_fmadd_ps(e_a_x, repulsion_x, check_y);

            threshold = _mm256_mul_ps(threshold, projection_factor_vec);

            mask = _mm256_cmp_ps(_mm256_mul_ps(check, minus1_vec), threshold, _CMP_GE_OQ);

            w = _mm256_blendv_ps(influencer_vec, one, mask);

            repulsion_x = _mm256_mul_ps(w, repulsion_x);
            repulsion_y = _mm256_mul_ps(w, repulsion_y);

            people_repulsion_x = _mm256_add_ps(people_repulsion_x, repulsion_x);
            people_repulsion_y = _mm256_add_ps(people_repulsion_y, repulsion_y);
          }
          _mm256_store_ps(people_repulsion_term + 8 * i, people_repulsion_x);
          _mm256_store_ps(people_repulsion_term + 8 * number_of_people + 8 * i, people_repulsion_y);
        }
      }
    }

    {
      __m256 social_force_x;
      __m256 social_force_y;

      // acceleration term
      __m256 actual_velocity_x;
      __m256 actual_velocity_y;
      __m256 desired_direction_x;
      __m256 desired_direction_y;

      __m256 desired_speed_value;
      __m256 v_delta_x;
      __m256 v_delta_y;
      __m256 inv_relax_time_vec = _mm256_set1_ps(INV_RELAX_TIME);
      // acceleration term

      // border repulsion term
      __m256 border;
      __m256 r_a_y;

      __m256 r_aB_y;

      __m256 r_aB_minus_y;

      __m256 r_aB_norm;

      __m256 mask_y;

      __m256 common_factor;
      __m256 exp;

      __m256 border_1;
      __m256 r_a_y_1;

      __m256 r_aB_y_1;

      __m256 r_aB_minus_y_1;

      __m256 r_aB_norm_1;

      __m256 mask_y_1;

      __m256 common_factor_1;
      __m256 exp_1;

      // update position
      __m256 prefered_velocity_x;
      __m256 prefered_velocity_y;

      __m256 position_x;
      __m256 position_y;

      __m256 prefered_velocity_2_x;
      __m256 prefered_velocity_2_y;

      __m256 prefered_velocity_norm_inv;
      __m256 prefered_velocity_norm;
      __m256 max_speed;
      __m256 control_value;
      __m256 mask;

      __m128d actual_speed_x;
      __m128d actual_speed_y;

      __m256 actual_speed_xy_256;

      __m256 control_value_prefered_velocity_xy_norm;
      __m256 control_value_prefered_velocity_xy_norm_inv;

      __m256 timestep_vec = _mm256_set1_ps(TIMESTEP);
      __m256 one = _mm256_set1_ps(1);

      //

      __m256 zero = _mm256_set1_ps(0);

      __m256 minus1 = _mm256_set1_ps(-1);

      __m256 r_vec_inv = _mm256_set1_ps(1 / R);
      __m256 minus_r_vec_inv = _mm256_set1_ps(-1 / R);
      __m256 u_alpha_b_vec = _mm256_set1_ps(U_ALPHA_B / R);

      __m256 exp_constant = _mm256_set1_ps(0.00006103515); // 1 / 16384

      border = _mm256_broadcast_ss(borders + 0);
      border_1 = _mm256_broadcast_ss(borders + 1);

      // border repulsion term

      __m256 unpacked;
      __m256 unpacked1;

      __m256 people_repulsion_x;
      __m256 people_repulsion_y;

      __m256 row_add0_x;
      __m256 row_add1_x;
      __m256 row_add2_x;
      __m256 row_add3_x;
      __m256 row_add4_x;
      __m256 row_add5_x;
      __m256 row_add6_x;
      __m256 row_add7_x;

      __m256 row_add0_y;
      __m256 row_add1_y;
      __m256 row_add2_y;
      __m256 row_add3_y;
      __m256 row_add4_y;
      __m256 row_add5_y;
      __m256 row_add6_y;
      __m256 row_add7_y;

      for (int p = 0; p < number_of_people - 7; p += 8)
      {

        /*************************************
        // acceleration term
        *************************************/

        actual_velocity_x = _mm256_load_ps(actual_velocity + p);
        actual_velocity_y = _mm256_load_ps(actual_velocity + number_of_people + p);
        desired_direction_x = _mm256_load_ps(desired_direction + p);
        desired_direction_y = _mm256_load_ps(desired_direction + number_of_people + p);
        desired_speed_value = _mm256_load_ps(desired_speed + p);

        v_delta_x = _mm256_fmsub_ps(desired_speed_value, desired_direction_x, actual_velocity_x);
        v_delta_y = _mm256_fmsub_ps(desired_speed_value, desired_direction_y, actual_velocity_y);

        v_delta_x = _mm256_mul_ps(v_delta_x, inv_relax_time_vec);
        v_delta_y = _mm256_mul_ps(v_delta_y, inv_relax_time_vec);

        social_force_x = v_delta_x;
        social_force_y = v_delta_y;

        // acceleration term

        /*************************************
        // border repulsion term
        *************************************/

        r_a_y = _mm256_load_ps(position + number_of_people + p);

        r_aB_y = _mm256_sub_ps(r_a_y, border);

        r_aB_minus_y = _mm256_mul_ps(r_aB_y, minus1);

        mask_y = _mm256_cmp_ps(r_aB_y, zero, _CMP_GE_OQ);
        r_aB_norm = _mm256_blendv_ps(r_aB_minus_y, r_aB_y, mask_y);

        exp = _mm256_mul_ps(r_aB_norm, minus_r_vec_inv);
        exp = exp_fast_vec_4(exp, one, exp_constant);

        common_factor = _mm256_mul_ps(u_alpha_b_vec, _mm256_rcp_ps(r_aB_norm));
        common_factor = _mm256_mul_ps(exp, common_factor);

        common_factor = _mm256_mul_ps(r_aB_y, common_factor);
        //
        r_aB_y_1 = _mm256_sub_ps(r_a_y, border_1);

        r_aB_minus_y_1 = _mm256_mul_ps(r_aB_y_1, minus1);

        mask_y_1 = _mm256_cmp_ps(r_aB_y_1, zero, _CMP_GE_OQ);
        r_aB_norm_1 = _mm256_blendv_ps(r_aB_minus_y_1, r_aB_y_1, mask_y_1);

        exp_1 = _mm256_mul_ps(r_aB_norm_1, minus_r_vec_inv);
        exp_1 = exp_fast_vec_4(exp_1, one, exp_constant);

        common_factor_1 = _mm256_mul_ps(u_alpha_b_vec, _mm256_rcp_ps(r_aB_norm_1));
        common_factor_1 = _mm256_mul_ps(exp_1, common_factor_1);

        common_factor_1 = _mm256_mul_ps(r_aB_y_1, common_factor_1);

        social_force_y = _mm256_add_ps(social_force_y, common_factor);
        social_force_y = _mm256_add_ps(social_force_y, common_factor_1);

        // border repulsion term

        /*************************************
        // add up people repulsion term
        *************************************/

        row_add0_x = _mm256_load_ps(people_repulsion_term + 8 * p);
        row_add1_x = _mm256_load_ps(people_repulsion_term + 8 * (p + 1));
        row_add2_x = _mm256_load_ps(people_repulsion_term + 8 * (p + 2));
        row_add3_x = _mm256_load_ps(people_repulsion_term + 8 * (p + 3));
        row_add4_x = _mm256_load_ps(people_repulsion_term + 8 * (p + 4));
        row_add5_x = _mm256_load_ps(people_repulsion_term + 8 * (p + 5));
        row_add6_x = _mm256_load_ps(people_repulsion_term + 8 * (p + 6));
        row_add7_x = _mm256_load_ps(people_repulsion_term + 8 * (p + 7));

        row_add0_x = _mm256_hadd_ps(row_add0_x, _mm256_permute2f128_ps(row_add0_x, row_add0_x, 3));
        row_add1_x = _mm256_hadd_ps(row_add1_x, _mm256_permute2f128_ps(row_add1_x, row_add1_x, 3));
        row_add2_x = _mm256_hadd_ps(row_add2_x, _mm256_permute2f128_ps(row_add2_x, row_add2_x, 3));
        row_add3_x = _mm256_hadd_ps(row_add3_x, _mm256_permute2f128_ps(row_add3_x, row_add3_x, 3));
        row_add4_x = _mm256_hadd_ps(row_add4_x, _mm256_permute2f128_ps(row_add4_x, row_add4_x, 3));
        row_add5_x = _mm256_hadd_ps(row_add5_x, _mm256_permute2f128_ps(row_add5_x, row_add5_x, 3));
        row_add6_x = _mm256_hadd_ps(row_add6_x, _mm256_permute2f128_ps(row_add6_x, row_add6_x, 3));
        row_add7_x = _mm256_hadd_ps(row_add7_x, _mm256_permute2f128_ps(row_add7_x, row_add7_x, 3));

        row_add0_x = _mm256_hadd_ps(row_add0_x, row_add0_x);
        row_add1_x = _mm256_hadd_ps(row_add1_x, row_add1_x);
        row_add2_x = _mm256_hadd_ps(row_add2_x, row_add2_x);
        row_add3_x = _mm256_hadd_ps(row_add3_x, row_add3_x);
        row_add4_x = _mm256_hadd_ps(row_add4_x, row_add4_x);
        row_add5_x = _mm256_hadd_ps(row_add5_x, row_add5_x);
        row_add6_x = _mm256_hadd_ps(row_add6_x, row_add6_x);
        row_add7_x = _mm256_hadd_ps(row_add7_x, row_add7_x);

        row_add0_x = _mm256_hadd_ps(row_add0_x, row_add0_x);
        row_add1_x = _mm256_hadd_ps(row_add1_x, row_add1_x);
        row_add2_x = _mm256_hadd_ps(row_add2_x, row_add2_x);
        row_add3_x = _mm256_hadd_ps(row_add3_x, row_add3_x);
        row_add4_x = _mm256_hadd_ps(row_add4_x, row_add4_x);
        row_add5_x = _mm256_hadd_ps(row_add5_x, row_add5_x);
        row_add6_x = _mm256_hadd_ps(row_add6_x, row_add6_x);
        row_add7_x = _mm256_hadd_ps(row_add7_x, row_add7_x);

        unpacked = _mm256_setzero_ps();
        unpacked = _mm256_blend_ps(unpacked, row_add0_x, 0b00000001);
        unpacked = _mm256_blend_ps(unpacked, row_add1_x, 0b00000010);
        unpacked = _mm256_blend_ps(unpacked, row_add2_x, 0b00000100);
        unpacked = _mm256_blend_ps(unpacked, row_add3_x, 0b00001000);
        unpacked = _mm256_blend_ps(unpacked, row_add4_x, 0b00010000);
        unpacked = _mm256_blend_ps(unpacked, row_add5_x, 0b00100000);
        unpacked = _mm256_blend_ps(unpacked, row_add6_x, 0b01000000);
        unpacked = _mm256_blend_ps(unpacked, row_add7_x, 0b10000000);

        social_force_x = _mm256_add_ps(unpacked, social_force_x);

        row_add0_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * p);
        row_add1_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * (p + 1));
        row_add2_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * (p + 2));
        row_add3_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * (p + 3));
        row_add4_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * (p + 4));
        row_add5_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * (p + 5));
        row_add6_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * (p + 6));
        row_add7_y = _mm256_load_ps(people_repulsion_term + 8 * number_of_people + 8 * (p + 7));

        row_add0_y = _mm256_hadd_ps(row_add0_y, _mm256_permute2f128_ps(row_add0_y, row_add0_y, 3));
        row_add1_y = _mm256_hadd_ps(row_add1_y, _mm256_permute2f128_ps(row_add1_y, row_add1_y, 3));
        row_add2_y = _mm256_hadd_ps(row_add2_y, _mm256_permute2f128_ps(row_add2_y, row_add2_y, 3));
        row_add3_y = _mm256_hadd_ps(row_add3_y, _mm256_permute2f128_ps(row_add3_y, row_add3_y, 3));
        row_add4_y = _mm256_hadd_ps(row_add4_y, _mm256_permute2f128_ps(row_add4_y, row_add4_y, 3));
        row_add5_y = _mm256_hadd_ps(row_add5_y, _mm256_permute2f128_ps(row_add5_y, row_add5_y, 3));
        row_add6_y = _mm256_hadd_ps(row_add6_y, _mm256_permute2f128_ps(row_add6_y, row_add6_y, 3));
        row_add7_y = _mm256_hadd_ps(row_add7_y, _mm256_permute2f128_ps(row_add7_y, row_add7_y, 3));

        row_add0_y = _mm256_hadd_ps(row_add0_y, row_add0_y);
        row_add1_y = _mm256_hadd_ps(row_add1_y, row_add1_y);
        row_add2_y = _mm256_hadd_ps(row_add2_y, row_add2_y);
        row_add3_y = _mm256_hadd_ps(row_add3_y, row_add3_y);
        row_add4_y = _mm256_hadd_ps(row_add4_y, row_add4_y);
        row_add5_y = _mm256_hadd_ps(row_add5_y, row_add5_y);
        row_add6_y = _mm256_hadd_ps(row_add6_y, row_add6_y);
        row_add7_y = _mm256_hadd_ps(row_add7_y, row_add7_y);

        row_add0_y = _mm256_hadd_ps(row_add0_y, row_add0_y);
        row_add1_y = _mm256_hadd_ps(row_add1_y, row_add1_y);
        row_add2_y = _mm256_hadd_ps(row_add2_y, row_add2_y);
        row_add3_y = _mm256_hadd_ps(row_add3_y, row_add3_y);
        row_add4_y = _mm256_hadd_ps(row_add4_y, row_add4_y);
        row_add5_y = _mm256_hadd_ps(row_add5_y, row_add5_y);
        row_add6_y = _mm256_hadd_ps(row_add6_y, row_add6_y);
        row_add7_y = _mm256_hadd_ps(row_add7_y, row_add7_y);

        unpacked = _mm256_setzero_ps();
        unpacked = _mm256_blend_ps(unpacked, row_add0_y, 0b00000001);
        unpacked = _mm256_blend_ps(unpacked, row_add1_y, 0b00000010);
        unpacked = _mm256_blend_ps(unpacked, row_add2_y, 0b00000100);
        unpacked = _mm256_blend_ps(unpacked, row_add3_y, 0b00001000);
        unpacked = _mm256_blend_ps(unpacked, row_add4_y, 0b00010000);
        unpacked = _mm256_blend_ps(unpacked, row_add5_y, 0b00100000);
        unpacked = _mm256_blend_ps(unpacked, row_add6_y, 0b01000000);
        unpacked = _mm256_blend_ps(unpacked, row_add7_y, 0b10000000);

        social_force_y = _mm256_add_ps(unpacked, social_force_y);

        // people repulsion term

        /*************************************
        // update position
        *************************************/

        prefered_velocity_x = _mm256_load_ps(actual_velocity + p);
        prefered_velocity_y = _mm256_load_ps(actual_velocity + number_of_people + p);

        prefered_velocity_x = _mm256_fmadd_ps(social_force_x, timestep_vec, prefered_velocity_x);
        prefered_velocity_y = _mm256_fmadd_ps(social_force_y, timestep_vec, prefered_velocity_y);

        // norm
        prefered_velocity_2_x = _mm256_mul_ps(prefered_velocity_x, prefered_velocity_x);
        prefered_velocity_norm_inv = _mm256_rsqrt_ps(_mm256_fmadd_ps(prefered_velocity_y, prefered_velocity_y, prefered_velocity_2_x));
        prefered_velocity_norm = _mm256_rcp_ps(prefered_velocity_norm_inv);

        max_speed = _mm256_load_ps(desired_max_speed + p);

        mask = _mm256_cmp_ps(prefered_velocity_norm, max_speed, _CMP_GT_OQ);
        control_value = _mm256_blendv_ps(one, _mm256_mul_ps(max_speed, prefered_velocity_norm_inv), mask);

        prefered_velocity_x = _mm256_mul_ps(prefered_velocity_x, control_value);
        prefered_velocity_y = _mm256_mul_ps(prefered_velocity_y, control_value);

        _mm256_store_ps(actual_velocity + p, prefered_velocity_x);
        _mm256_store_ps(actual_velocity + number_of_people + p, prefered_velocity_y);

        control_value_prefered_velocity_xy_norm = _mm256_mul_ps(control_value, prefered_velocity_norm);
        control_value_prefered_velocity_xy_norm_inv = _mm256_rcp_ps(control_value_prefered_velocity_xy_norm);

        _mm256_store_ps(speed + p, control_value_prefered_velocity_xy_norm);

        _mm256_store_ps(desired_direction + p, _mm256_mul_ps(prefered_velocity_x, control_value_prefered_velocity_xy_norm_inv));
        _mm256_store_ps(desired_direction + number_of_people + p, _mm256_mul_ps(prefered_velocity_y, control_value_prefered_velocity_xy_norm_inv));

        position_x = _mm256_load_ps(position + p);
        position_y = _mm256_load_ps(position + number_of_people + p);

        _mm256_store_ps(position + p, _mm256_fmadd_ps(prefered_velocity_x, timestep_vec, position_x));
        _mm256_store_ps(position + number_of_people + p, _mm256_fmadd_ps(prefered_velocity_y, timestep_vec, position_y));
      }
    }
    CONSOLE_PRINT(("Finished iteration %d\n", (step + 1)));
  }

  CONSOLE_PRINT(("Simulation terminated\n"));
}
