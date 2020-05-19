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


void simulation_basic_vectorize_4_double(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination, double *borders, double *actual_velocity, double *acceleration_term,
                                         double *people_repulsion_term, double *border_repulsion_term, double *social_force, double *desired_speed, double *desired_max_speed)
{
  // start simulation
  CONSOLE_PRINT(("Start simulation with %d persons\n", number_of_people));

  // simulate steps
  for (int step = 0; step < n_timesteps; step++)
  {
    {
      __m256d current_x;
      __m256d current_y;
      __m256d target_x;
      __m256d target_y;
      __m256d delta_x;
      __m256d delta_y;
      __m256d delta_x_2;
      __m256d delta_y_2;
      __m256d normalizer;

      __m256d result_x;
      __m256d result_y;

      __m256d one = _mm256_set1_pd(1);

      // iterate over all persons and update desired_direction
      for (int i = 0; i < number_of_people - 3; i += 4)
      {
        // get current position and target
        current_x = _mm256_load_pd(position + i); // now xy positions for two persons in register
        target_x = _mm256_load_pd(final_destination + i);
        current_y = _mm256_load_pd(position + number_of_people + i); // now xy positions for two persons in register
        target_y = _mm256_load_pd(final_destination + number_of_people + i);

        // compute differences
        delta_x = _mm256_sub_pd(target_x, current_x);
        delta_y = _mm256_sub_pd(target_y, current_y);

        // compute norm
        delta_x_2 = _mm256_mul_pd(delta_x, delta_x); // square each entry
        normalizer = _mm256_sqrt_pd(_mm256_fmadd_pd(delta_y, delta_y, delta_x_2));

        result_x = _mm256_div_pd(delta_x, normalizer);
        result_y = _mm256_div_pd(delta_y, normalizer);

        _mm256_store_pd(desired_direction + i, result_x);
        _mm256_store_pd(desired_direction + number_of_people + i, result_y);
      }
    }
    {
      __m256d position_i_x;
      __m256d position_i_y;

      __m256d position_j_x;
      __m256d position_j_y;

      __m256d r_ab_x;
      __m256d r_ab_y;
      __m256d e_a_x;
      __m256d e_a_y;
      __m256d e_b_x;
      __m256d e_b_y;

      __m256d vb;
      __m256d delta_b;

      __m256d r_ab_2_x;
      __m256d r_ab_2_y;
      __m256d r_ab_norm;
      __m256d r_ab_norm_inv;

      __m256d r_ab_me_x;
      __m256d r_ab_me_y;

      __m256d r_ab_me_2_x;
      __m256d r_ab_me_2_y;
      __m256d r_ab_me_norm;
      __m256d r_ab_me_norm_inv;

      __m256d norm_sum;
      __m256d repulsion_x;
      __m256d repulsion_y;

      __m256d norm_sum_2, delta_b_2;
      __m256d b;
      __m256d exp;
      __m256d common_factor;
      __m256d check;
      __m256d check_x;
      __m256d check_y;

      __m256d repulsion_2_x;
      __m256d repulsion_2_y;
      __m256d repulsion_norm;
      __m256d threshold;
      __m256d w;
      __m256d mask;

      __m256d timestep_vec = _mm256_set1_pd(-TIMESTEP);
      __m256d minus_sigma_inv_vec = _mm256_set1_pd(-1.0 / SIGMA);
      __m256d div_factor_vec = _mm256_set1_pd(DIV_FACTOR);
      __m256d projection_factor_vec = _mm256_set1_pd(PROJECTION_FACTOR);
      __m256d influencer_vec = _mm256_set1_pd(INFLUENCE);

      __m256d one = _mm256_set1_pd(1);
      __m256d half_vec = _mm256_set1_pd(0.5);
      __m256d minus1_vec = _mm256_set1_pd(-1);
      __m256d eps = _mm256_set1_pd(1e-12);

      __m256d current_mask;

      // diagonal
      int n_eight = number_of_people / 4;
      for (int k = 0; k < n_eight; k++)
      {
        position_j_x = _mm256_load_pd(position + 4 * k);
        position_j_y = _mm256_load_pd(position + number_of_people + 4 * k);

        e_b_x = _mm256_load_pd(desired_direction + 4 * k);
        e_b_y = _mm256_load_pd(desired_direction + number_of_people + 4 * k);

        vb = _mm256_load_pd(speed + 4 * k);

        for (int i = 4 * k; i < 4 * (k + 1); i++)
        {
          position_i_x = _mm256_broadcast_sd(position + i);
          position_i_y = _mm256_broadcast_sd(position + number_of_people + i);
          e_a_x = _mm256_broadcast_sd(desired_direction + i);
          e_a_y = _mm256_broadcast_sd(desired_direction + number_of_people + i);

          r_ab_x = _mm256_sub_pd(position_i_x, position_j_x);
          r_ab_y = _mm256_sub_pd(position_i_y, position_j_y);

          delta_b = _mm256_mul_pd(vb, timestep_vec);

          // compute norm r_ab
          r_ab_2_x = _mm256_mul_pd(r_ab_x, r_ab_x);
          r_ab_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_y, r_ab_y, r_ab_2_x));

          r_ab_me_x = _mm256_fmadd_pd(delta_b, e_b_x, r_ab_x);
          r_ab_me_y = _mm256_fmadd_pd(delta_b, e_b_y, r_ab_y);

          // compute norm r_ab_me
          r_ab_me_2_x = _mm256_mul_pd(r_ab_me_x, r_ab_me_x);
          r_ab_me_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

          // sum up norms
          norm_sum = _mm256_add_pd((r_ab_norm), (r_ab_me_norm));

          repulsion_x = _mm256_div_pd(r_ab_x, r_ab_norm);
          repulsion_x = _mm256_fmadd_pd(r_ab_me_x, _mm256_div_pd(one, r_ab_me_norm), repulsion_x);

          repulsion_y = _mm256_div_pd(r_ab_y, r_ab_norm);
          repulsion_y = _mm256_fmadd_pd(r_ab_me_y, _mm256_div_pd(one, r_ab_me_norm), repulsion_y);

          delta_b_2 = _mm256_mul_pd(delta_b, delta_b);
          b = _mm256_sqrt_pd(_mm256_fmsub_pd(norm_sum, norm_sum, delta_b_2));
          b = _mm256_mul_pd(b, half_vec);

          exp = _mm256_mul_pd(b, minus_sigma_inv_vec);
          exp = exp_fast_vec_double(exp);

          common_factor = _mm256_mul_pd(norm_sum, div_factor_vec);
          common_factor = _mm256_div_pd(common_factor, b);
          common_factor = _mm256_mul_pd(exp, common_factor);

          repulsion_x = _mm256_mul_pd(repulsion_x, common_factor);
          repulsion_y = _mm256_mul_pd(repulsion_y, common_factor);

          // compute norm r_ab
          repulsion_2_y = _mm256_mul_pd(repulsion_y, repulsion_y);
          threshold = _mm256_sqrt_pd(_mm256_fmadd_pd(repulsion_x, repulsion_x, repulsion_2_y));

          check_y = _mm256_mul_pd(e_a_y, repulsion_y);
          check = _mm256_fmadd_pd(e_a_x, repulsion_x, check_y);

          threshold = _mm256_mul_pd(threshold, projection_factor_vec);

          mask = _mm256_cmp_pd(_mm256_mul_pd(check, minus1_vec), threshold, _CMP_GE_OQ);

          w = _mm256_blendv_pd(influencer_vec, one, mask);

          repulsion_x = _mm256_mul_pd(w, repulsion_x);
          repulsion_y = _mm256_mul_pd(w, repulsion_y);

          switch (i % 4)
          {
          case 0:
            repulsion_x = _mm256_blend_pd(_mm256_setzero_pd(), repulsion_x, 0b1110);
            repulsion_y = _mm256_blend_pd(_mm256_setzero_pd(), repulsion_y, 0b1110);

            break;
          case 1:
            repulsion_x = _mm256_blend_pd(_mm256_setzero_pd(), repulsion_x, 0b1101);
            repulsion_y = _mm256_blend_pd(_mm256_setzero_pd(), repulsion_y, 0b1101);
            break;
          case 2:
            repulsion_x = _mm256_blend_pd(_mm256_setzero_pd(), repulsion_x, 0b1011);
            repulsion_y = _mm256_blend_pd(_mm256_setzero_pd(), repulsion_y, 0b1011);
            break;
          case 3:
            repulsion_x = _mm256_blend_pd(_mm256_setzero_pd(), repulsion_x, 0b0111);
            repulsion_y = _mm256_blend_pd(_mm256_setzero_pd(), repulsion_y, 0b0111);
            break;
            }

          _mm256_store_pd(people_repulsion_term + 4 * i, repulsion_x);
          _mm256_store_pd(people_repulsion_term + 4 * number_of_people + 4 * i, repulsion_y);
        }
      }
      // first block row

      for (int k = 1; k < 2; k++)
      {
        for (int i = (k - 1) * 4; i < k * 4; i++)
        {
          position_i_x = _mm256_broadcast_sd(position + i);
          position_i_y = _mm256_broadcast_sd(position + number_of_people + i);
          e_a_x = _mm256_broadcast_sd(desired_direction + i);
          e_a_y = _mm256_broadcast_sd(desired_direction + number_of_people + i);

          __m256d people_repulsion_x = _mm256_load_pd(people_repulsion_term + 4 * i);
          __m256d people_repulsion_y = _mm256_load_pd(people_repulsion_term + 4 * number_of_people + 4 * i);

          for (int j = k * 4; j < number_of_people - 3; j += 4)
          {
            position_j_x = _mm256_load_pd(position + j);
            position_j_y = _mm256_load_pd(position + number_of_people + j);

            r_ab_x = _mm256_sub_pd(position_i_x, position_j_x);
            r_ab_y = _mm256_sub_pd(position_i_y, position_j_y);

            e_b_x = _mm256_load_pd(desired_direction + j);
            e_b_y = _mm256_load_pd(desired_direction + number_of_people + j);

            vb = _mm256_load_pd(speed + j);

            delta_b = _mm256_mul_pd(vb, timestep_vec);

            // compute norm r_ab
            r_ab_2_x = _mm256_mul_pd(r_ab_x, r_ab_x);
            r_ab_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_y, r_ab_y, r_ab_2_x));

            r_ab_me_x = _mm256_fmadd_pd(delta_b, e_b_x, r_ab_x);
            r_ab_me_y = _mm256_fmadd_pd(delta_b, e_b_y, r_ab_y);

            // compute norm r_ab_me
            r_ab_me_2_x = _mm256_mul_pd(r_ab_me_x, r_ab_me_x);
            r_ab_me_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

            // sum up norms
            norm_sum = _mm256_add_pd((r_ab_norm), (r_ab_me_norm));

            repulsion_x = _mm256_div_pd(r_ab_x, r_ab_norm);
            repulsion_x = _mm256_fmadd_pd(r_ab_me_x, _mm256_div_pd(one, r_ab_me_norm), repulsion_x);

            repulsion_y = _mm256_div_pd(r_ab_y, r_ab_norm);
            repulsion_y = _mm256_fmadd_pd(r_ab_me_y, _mm256_div_pd(one, r_ab_me_norm), repulsion_y);

            delta_b_2 = _mm256_mul_pd(delta_b, delta_b);
            b = _mm256_sqrt_pd(_mm256_fmsub_pd(norm_sum, norm_sum, delta_b_2));
            b = _mm256_mul_pd(b, half_vec);

            exp = _mm256_mul_pd(b, minus_sigma_inv_vec);
            exp = exp_fast_vec_double(exp);

            common_factor = _mm256_mul_pd(norm_sum, div_factor_vec);
            common_factor = _mm256_div_pd(common_factor, b);
            common_factor = _mm256_mul_pd(exp, common_factor);

            repulsion_x = _mm256_mul_pd(repulsion_x, common_factor);
            repulsion_y = _mm256_mul_pd(repulsion_y, common_factor);

            // compute norm r_ab
            repulsion_2_y = _mm256_mul_pd(repulsion_y, repulsion_y);
            threshold = _mm256_sqrt_pd(_mm256_fmadd_pd(repulsion_x, repulsion_x, repulsion_2_y));

            check_y = _mm256_mul_pd(e_a_y, repulsion_y);
            check = _mm256_fmadd_pd(e_a_x, repulsion_x, check_y);

            threshold = _mm256_mul_pd(threshold, projection_factor_vec);

            mask = _mm256_cmp_pd(_mm256_mul_pd(check, minus1_vec), threshold, _CMP_GE_OQ);

            w = _mm256_blendv_pd(influencer_vec, one, mask);

            repulsion_x = _mm256_mul_pd(w, repulsion_x);
            repulsion_y = _mm256_mul_pd(w, repulsion_y);

            people_repulsion_x = _mm256_add_pd(people_repulsion_x, repulsion_x);
            people_repulsion_y = _mm256_add_pd(people_repulsion_y, repulsion_y);
          }
          _mm256_store_pd(people_repulsion_term + 4 * i, people_repulsion_x);
          _mm256_store_pd(people_repulsion_term + 4 * number_of_people + 4 * i, people_repulsion_y);
        }
      }

      // upper triangle

      for (int k = 2; k < n_eight; k++)
      {
        for (int i = (k - 1) * 4; i < k * 4; i++)
        {
          position_i_x = _mm256_broadcast_sd(position + i);
          position_i_y = _mm256_broadcast_sd(position + number_of_people + i);
          e_a_x = _mm256_broadcast_sd(desired_direction + i);
          e_a_y = _mm256_broadcast_sd(desired_direction + number_of_people + i);

          __m256d people_repulsion_x = _mm256_load_pd(people_repulsion_term + 4 * i);
          __m256d people_repulsion_y = _mm256_load_pd(people_repulsion_term + 4 * number_of_people + 4 * i);

          for (int j = 0; j < (k - 1) * 4; j += 4)
          {
            position_j_x = _mm256_load_pd(position + j);
            position_j_y = _mm256_load_pd(position + number_of_people + j);

            r_ab_x = _mm256_sub_pd(position_i_x, position_j_x);
            r_ab_y = _mm256_sub_pd(position_i_y, position_j_y);

            e_b_x = _mm256_load_pd(desired_direction + j);
            e_b_y = _mm256_load_pd(desired_direction + number_of_people + j);

            vb = _mm256_load_pd(speed + j);

            delta_b = _mm256_mul_pd(vb, timestep_vec);

            // compute norm r_ab
            r_ab_2_x = _mm256_mul_pd(r_ab_x, r_ab_x);
            r_ab_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_y, r_ab_y, r_ab_2_x));

            r_ab_me_x = _mm256_fmadd_pd(delta_b, e_b_x, r_ab_x);
            r_ab_me_y = _mm256_fmadd_pd(delta_b, e_b_y, r_ab_y);

            // compute norm r_ab_me
            r_ab_me_2_x = _mm256_mul_pd(r_ab_me_x, r_ab_me_x);
            r_ab_me_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

            // sum up norms
            norm_sum = _mm256_add_pd((r_ab_norm), (r_ab_me_norm));

            repulsion_x = _mm256_div_pd(r_ab_x, r_ab_norm);
            repulsion_x = _mm256_fmadd_pd(r_ab_me_x, _mm256_div_pd(one, r_ab_me_norm), repulsion_x);

            repulsion_y = _mm256_div_pd(r_ab_y, r_ab_norm);
            repulsion_y = _mm256_fmadd_pd(r_ab_me_y, _mm256_div_pd(one, r_ab_me_norm), repulsion_y);

            delta_b_2 = _mm256_mul_pd(delta_b, delta_b);
            b = _mm256_sqrt_pd(_mm256_fmsub_pd(norm_sum, norm_sum, delta_b_2));
            b = _mm256_mul_pd(b, half_vec);

            exp = _mm256_mul_pd(b, minus_sigma_inv_vec);
            exp = exp_fast_vec_double(exp);

            common_factor = _mm256_mul_pd(norm_sum, div_factor_vec);
            common_factor = _mm256_div_pd(common_factor, b);
            common_factor = _mm256_mul_pd(exp, common_factor);

            repulsion_x = _mm256_mul_pd(repulsion_x, common_factor);
            repulsion_y = _mm256_mul_pd(repulsion_y, common_factor);

            // compute norm r_ab
            repulsion_2_y = _mm256_mul_pd(repulsion_y, repulsion_y);
            threshold = _mm256_sqrt_pd(_mm256_fmadd_pd(repulsion_x, repulsion_x, repulsion_2_y));

            check_y = _mm256_mul_pd(e_a_y, repulsion_y);
            check = _mm256_fmadd_pd(e_a_x, repulsion_x, check_y);

            threshold = _mm256_mul_pd(threshold, projection_factor_vec);

            mask = _mm256_cmp_pd(_mm256_mul_pd(check, minus1_vec), threshold, _CMP_GE_OQ);

            w = _mm256_blendv_pd(influencer_vec, one, mask);

            repulsion_x = _mm256_mul_pd(w, repulsion_x);
            repulsion_y = _mm256_mul_pd(w, repulsion_y);

            people_repulsion_x = _mm256_add_pd(people_repulsion_x, repulsion_x);
            people_repulsion_y = _mm256_add_pd(people_repulsion_y, repulsion_y);
          }

          for (int j = k * 4; j < number_of_people - 3; j += 4)
          {
            position_j_x = _mm256_load_pd(position + j);
            position_j_y = _mm256_load_pd(position + number_of_people + j);

            r_ab_x = _mm256_sub_pd(position_i_x, position_j_x);
            r_ab_y = _mm256_sub_pd(position_i_y, position_j_y);

            e_b_x = _mm256_load_pd(desired_direction + j);
            e_b_y = _mm256_load_pd(desired_direction + number_of_people + j);

            vb = _mm256_load_pd(speed + j);

            delta_b = _mm256_mul_pd(vb, timestep_vec);

            // compute norm r_ab
            r_ab_2_x = _mm256_mul_pd(r_ab_x, r_ab_x);
            r_ab_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_y, r_ab_y, r_ab_2_x));

            r_ab_me_x = _mm256_fmadd_pd(delta_b, e_b_x, r_ab_x);
            r_ab_me_y = _mm256_fmadd_pd(delta_b, e_b_y, r_ab_y);

            // compute norm r_ab_me
            r_ab_me_2_x = _mm256_mul_pd(r_ab_me_x, r_ab_me_x);
            r_ab_me_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

            // sum up norms
            norm_sum = _mm256_add_pd((r_ab_norm), (r_ab_me_norm));

            repulsion_x = _mm256_div_pd(r_ab_x, r_ab_norm);
            repulsion_x = _mm256_fmadd_pd(r_ab_me_x, _mm256_div_pd(one, r_ab_me_norm), repulsion_x);

            repulsion_y = _mm256_div_pd(r_ab_y, r_ab_norm);
            repulsion_y = _mm256_fmadd_pd(r_ab_me_y, _mm256_div_pd(one, r_ab_me_norm), repulsion_y);

            delta_b_2 = _mm256_mul_pd(delta_b, delta_b);
            b = _mm256_sqrt_pd(_mm256_fmsub_pd(norm_sum, norm_sum, delta_b_2));
            b = _mm256_mul_pd(b, half_vec);

            exp = _mm256_mul_pd(b, minus_sigma_inv_vec);
            exp = exp_fast_vec_double(exp);

            common_factor = _mm256_mul_pd(norm_sum, div_factor_vec);
            common_factor = _mm256_div_pd(common_factor, b);
            common_factor = _mm256_mul_pd(exp, common_factor);

            repulsion_x = _mm256_mul_pd(repulsion_x, common_factor);
            repulsion_y = _mm256_mul_pd(repulsion_y, common_factor);

            // compute norm r_ab
            repulsion_2_y = _mm256_mul_pd(repulsion_y, repulsion_y);
            threshold = _mm256_sqrt_pd(_mm256_fmadd_pd(repulsion_x, repulsion_x, repulsion_2_y));

            check_y = _mm256_mul_pd(e_a_y, repulsion_y);
            check = _mm256_fmadd_pd(e_a_x, repulsion_x, check_y);

            threshold = _mm256_mul_pd(threshold, projection_factor_vec);

            mask = _mm256_cmp_pd(_mm256_mul_pd(check, minus1_vec), threshold, _CMP_GE_OQ);

            w = _mm256_blendv_pd(influencer_vec, one, mask);

            repulsion_x = _mm256_mul_pd(w, repulsion_x);
            repulsion_y = _mm256_mul_pd(w, repulsion_y);

            people_repulsion_x = _mm256_add_pd(people_repulsion_x, repulsion_x);
            people_repulsion_y = _mm256_add_pd(people_repulsion_y, repulsion_y);
          }

          _mm256_store_pd(people_repulsion_term + 4 * i, people_repulsion_x);
          _mm256_store_pd(people_repulsion_term + 4 * number_of_people + 4 * i, people_repulsion_y);
        }
      }
      // last row of block
      for (int k = n_eight - 1; k < n_eight; k++)
      {
        for (int i = k * 4; i < (k + 1) * 4; i++)
        {
          position_i_x = _mm256_broadcast_sd(position + i);
          position_i_y = _mm256_broadcast_sd(position + number_of_people + i);
          e_a_x = _mm256_broadcast_sd(desired_direction + i);
          e_a_y = _mm256_broadcast_sd(desired_direction + number_of_people + i);

          __m256d people_repulsion_x = _mm256_load_pd(people_repulsion_term + 4 * i);
          __m256d people_repulsion_y = _mm256_load_pd(people_repulsion_term + 4 * number_of_people + 4 * i);

          for (int j = 0; j < k * 4; j += 4)
          {
            position_j_x = _mm256_load_pd(position + j);
            position_j_y = _mm256_load_pd(position + number_of_people + j);

            r_ab_x = _mm256_sub_pd(position_i_x, position_j_x);
            r_ab_y = _mm256_sub_pd(position_i_y, position_j_y);

            e_b_x = _mm256_load_pd(desired_direction + j);
            e_b_y = _mm256_load_pd(desired_direction + number_of_people + j);

            vb = _mm256_load_pd(speed + j);

            delta_b = _mm256_mul_pd(vb, timestep_vec);

            // compute norm r_ab
            r_ab_2_x = _mm256_mul_pd(r_ab_x, r_ab_x);
            r_ab_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_y, r_ab_y, r_ab_2_x));

            r_ab_me_x = _mm256_fmadd_pd(delta_b, e_b_x, r_ab_x);
            r_ab_me_y = _mm256_fmadd_pd(delta_b, e_b_y, r_ab_y);

            // compute norm r_ab_me
            r_ab_me_2_x = _mm256_mul_pd(r_ab_me_x, r_ab_me_x);
            r_ab_me_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(r_ab_me_y, r_ab_me_y, r_ab_me_2_x));

            // sum up norms
            norm_sum = _mm256_add_pd((r_ab_norm), (r_ab_me_norm));

            repulsion_x = _mm256_div_pd(r_ab_x, r_ab_norm);
            repulsion_x = _mm256_fmadd_pd(r_ab_me_x, _mm256_div_pd(one, r_ab_me_norm), repulsion_x);

            repulsion_y = _mm256_div_pd(r_ab_y, r_ab_norm);
            repulsion_y = _mm256_fmadd_pd(r_ab_me_y, _mm256_div_pd(one, r_ab_me_norm), repulsion_y);

            delta_b_2 = _mm256_mul_pd(delta_b, delta_b);
            b = _mm256_sqrt_pd(_mm256_fmsub_pd(norm_sum, norm_sum, delta_b_2));
            b = _mm256_mul_pd(b, half_vec);

            exp = _mm256_mul_pd(b, minus_sigma_inv_vec);
            exp = exp_fast_vec_double(exp);

            common_factor = _mm256_mul_pd(norm_sum, div_factor_vec);
            common_factor = _mm256_div_pd(common_factor, b);
            common_factor = _mm256_mul_pd(exp, common_factor);

            repulsion_x = _mm256_mul_pd(repulsion_x, common_factor);
            repulsion_y = _mm256_mul_pd(repulsion_y, common_factor);

            // compute norm r_ab
            repulsion_2_y = _mm256_mul_pd(repulsion_y, repulsion_y);
            threshold = _mm256_sqrt_pd(_mm256_fmadd_pd(repulsion_x, repulsion_x, repulsion_2_y));

            check_y = _mm256_mul_pd(e_a_y, repulsion_y);
            check = _mm256_fmadd_pd(e_a_x, repulsion_x, check_y);

            threshold = _mm256_mul_pd(threshold, projection_factor_vec);

            mask = _mm256_cmp_pd(_mm256_mul_pd(check, minus1_vec), threshold, _CMP_GE_OQ);

            w = _mm256_blendv_pd(influencer_vec, one, mask);

            repulsion_x = _mm256_mul_pd(w, repulsion_x);
            repulsion_y = _mm256_mul_pd(w, repulsion_y);

            people_repulsion_x = _mm256_add_pd(people_repulsion_x, repulsion_x);
            people_repulsion_y = _mm256_add_pd(people_repulsion_y, repulsion_y);
          }
          _mm256_store_pd(people_repulsion_term + 4 * i, people_repulsion_x);
          _mm256_store_pd(people_repulsion_term + 4 * number_of_people + 4 * i, people_repulsion_y);
        }
      }
    }

    {
      __m256d social_force_x;
      __m256d social_force_y;

      // acceleration term
      __m256d actual_velocity_x;
      __m256d actual_velocity_y;
      __m256d desired_direction_x;
      __m256d desired_direction_y;

      __m256d desired_speed_value;
      __m256d v_delta_x;
      __m256d v_delta_y;
      __m256d inv_relax_time_vec = _mm256_set1_pd(INV_RELAX_TIME);
      // acceleration term

      // border repulsion term
      __m256d border;
      __m256d r_a_y;

      __m256d r_aB_y;

      __m256d r_aB_minus_y;

      __m256d r_aB_norm;

      __m256d mask_y;

      __m256d common_factor;
      __m256d exp;

      __m256d border_1;
      __m256d r_a_y_1;

      __m256d r_aB_y_1;

      __m256d r_aB_minus_y_1;

      __m256d r_aB_norm_1;

      __m256d mask_y_1;

      __m256d common_factor_1;
      __m256d exp_1;

      // update position
      __m256d prefered_velocity_x;
      __m256d prefered_velocity_y;

      __m256d position_x;
      __m256d position_y;

      __m256d prefered_velocity_2_x;
      __m256d prefered_velocity_2_y;

      __m256d prefered_velocity_norm_inv;
      __m256d prefered_velocity_norm;
      __m256d max_speed;
      __m256d control_value;
      __m256d mask;

      __m128d actual_speed_x;
      __m128d actual_speed_y;

      __m256d actual_speed_xy_256;

      __m256d control_value_prefered_velocity_xy_norm;
      __m256d control_value_prefered_velocity_xy_norm_inv;

      __m256d timestep_vec = _mm256_set1_pd(TIMESTEP);
      __m256d one = _mm256_set1_pd(1);

      //

      __m256d zero = _mm256_set1_pd(0);

      __m256d minus1 = _mm256_set1_pd(-1);

      __m256d r_vec_inv = _mm256_set1_pd(1 / R);
      __m256d minus_r_vec_inv = _mm256_set1_pd(-1 / R);
      __m256d u_alpha_b_vec = _mm256_set1_pd(U_ALPHA_B / R);


      border = _mm256_broadcast_sd(borders + 0);
      border_1 = _mm256_broadcast_sd(borders + 1);

      // border repulsion term

      __m256d unpacked;
      __m256d unpacked1;

      __m256d people_repulsion_x;
      __m256d people_repulsion_y;

      __m256d row_add0_x;
      __m256d row_add1_x;
      __m256d row_add2_x;
      __m256d row_add3_x;
      __m256d row_add4_x;
      __m256d row_add5_x;
      __m256d row_add6_x;
      __m256d row_add7_x;

      __m256d row_add0_y;
      __m256d row_add1_y;
      __m256d row_add2_y;
      __m256d row_add3_y;
      __m256d row_add4_y;
      __m256d row_add5_y;
      __m256d row_add6_y;
      __m256d row_add7_y;

      for (int p = 0; p < number_of_people - 3; p += 4)
      {

        /*************************************
        // acceleration term
        *************************************/

        actual_velocity_x = _mm256_load_pd(actual_velocity + p);
        actual_velocity_y = _mm256_load_pd(actual_velocity + number_of_people + p);
        desired_direction_x = _mm256_load_pd(desired_direction + p);
        desired_direction_y = _mm256_load_pd(desired_direction + number_of_people + p);
        desired_speed_value = _mm256_load_pd(desired_speed + p);

        v_delta_x = _mm256_fmsub_pd(desired_speed_value, desired_direction_x, actual_velocity_x);
        v_delta_y = _mm256_fmsub_pd(desired_speed_value, desired_direction_y, actual_velocity_y);

        v_delta_x = _mm256_mul_pd(v_delta_x, inv_relax_time_vec);
        v_delta_y = _mm256_mul_pd(v_delta_y, inv_relax_time_vec);

        social_force_x = v_delta_x;
        social_force_y = v_delta_y;

        // acceleration term

        /*************************************
        // border repulsion term
        *************************************/

        r_a_y = _mm256_load_pd(position + number_of_people + p);

        r_aB_y = _mm256_sub_pd(r_a_y, border);

        r_aB_minus_y = _mm256_mul_pd(r_aB_y, minus1);

        mask_y = _mm256_cmp_pd(r_aB_y, zero, _CMP_GE_OQ);
        r_aB_norm = _mm256_blendv_pd(r_aB_minus_y, r_aB_y, mask_y);

        exp = _mm256_mul_pd(r_aB_norm, minus_r_vec_inv);
        exp = exp_fast_vec_double(exp);

        common_factor = _mm256_div_pd(u_alpha_b_vec, (r_aB_norm));
        common_factor = _mm256_mul_pd(exp, common_factor);

        common_factor = _mm256_mul_pd(r_aB_y, common_factor);
        //
        r_aB_y_1 = _mm256_sub_pd(r_a_y, border_1);

        r_aB_minus_y_1 = _mm256_mul_pd(r_aB_y_1, minus1);

        mask_y_1 = _mm256_cmp_pd(r_aB_y_1, zero, _CMP_GE_OQ);
        r_aB_norm_1 = _mm256_blendv_pd(r_aB_minus_y_1, r_aB_y_1, mask_y_1);

        exp_1 = _mm256_mul_pd(r_aB_norm_1, minus_r_vec_inv);
        exp_1 = exp_fast_vec_double(exp_1);

        common_factor_1 = _mm256_div_pd(u_alpha_b_vec, (r_aB_norm_1));
        common_factor_1 = _mm256_mul_pd(exp_1, common_factor_1);

        common_factor_1 = _mm256_mul_pd(r_aB_y_1, common_factor_1);

        social_force_y = _mm256_add_pd(social_force_y, common_factor);
        social_force_y = _mm256_add_pd(social_force_y, common_factor_1);

        // border repulsion term

        /*************************************
        // add up people repulsion term
        *************************************/

        row_add0_x = _mm256_load_pd(people_repulsion_term + 4 * p);
        row_add1_x = _mm256_load_pd(people_repulsion_term + 4 * (p + 1));
        row_add2_x = _mm256_load_pd(people_repulsion_term + 4 * (p + 2));
        row_add3_x = _mm256_load_pd(people_repulsion_term + 4 * (p + 3));
       
        row_add0_x = _mm256_hadd_pd(row_add0_x, _mm256_permute2f128_pd(row_add0_x, row_add0_x, 3));
        row_add1_x = _mm256_hadd_pd(row_add1_x, _mm256_permute2f128_pd(row_add1_x, row_add1_x, 3));
        row_add2_x = _mm256_hadd_pd(row_add2_x, _mm256_permute2f128_pd(row_add2_x, row_add2_x, 3));
        row_add3_x = _mm256_hadd_pd(row_add3_x, _mm256_permute2f128_pd(row_add3_x, row_add3_x, 3));
      
        row_add0_x = _mm256_hadd_pd(row_add0_x, row_add0_x);
        row_add1_x = _mm256_hadd_pd(row_add1_x, row_add1_x);
        row_add2_x = _mm256_hadd_pd(row_add2_x, row_add2_x);
        row_add3_x = _mm256_hadd_pd(row_add3_x, row_add3_x);
      
        unpacked = _mm256_setzero_pd();
        unpacked = _mm256_blend_pd(unpacked, row_add0_x, 0b0001);
        unpacked = _mm256_blend_pd(unpacked, row_add1_x, 0b0010);
        unpacked = _mm256_blend_pd(unpacked, row_add2_x, 0b0100);
        unpacked = _mm256_blend_pd(unpacked, row_add3_x, 0b1000);
       

        social_force_x = _mm256_add_pd(unpacked, social_force_x);

        row_add0_y = _mm256_load_pd(people_repulsion_term + 4 * number_of_people + 4 * p);
        row_add1_y = _mm256_load_pd(people_repulsion_term + 4 * number_of_people + 4 * (p + 1));
        row_add2_y = _mm256_load_pd(people_repulsion_term + 4 * number_of_people + 4 * (p + 2));
        row_add3_y = _mm256_load_pd(people_repulsion_term + 4 * number_of_people + 4 * (p + 3));

        row_add0_y = _mm256_hadd_pd(row_add0_y, _mm256_permute2f128_pd(row_add0_y, row_add0_y, 3));
        row_add1_y = _mm256_hadd_pd(row_add1_y, _mm256_permute2f128_pd(row_add1_y, row_add1_y, 3));
        row_add2_y = _mm256_hadd_pd(row_add2_y, _mm256_permute2f128_pd(row_add2_y, row_add2_y, 3));
        row_add3_y = _mm256_hadd_pd(row_add3_y, _mm256_permute2f128_pd(row_add3_y, row_add3_y, 3));

        row_add0_y = _mm256_hadd_pd(row_add0_y, row_add0_y);
        row_add1_y = _mm256_hadd_pd(row_add1_y, row_add1_y);
        row_add2_y = _mm256_hadd_pd(row_add2_y, row_add2_y);
        row_add3_y = _mm256_hadd_pd(row_add3_y, row_add3_y);

        unpacked = _mm256_setzero_pd();
        unpacked = _mm256_blend_pd(unpacked, row_add0_y, 0b0001);
        unpacked = _mm256_blend_pd(unpacked, row_add1_y, 0b0010);
        unpacked = _mm256_blend_pd(unpacked, row_add2_y, 0b0100);
        unpacked = _mm256_blend_pd(unpacked, row_add3_y, 0b1000);


        social_force_y = _mm256_add_pd(unpacked, social_force_y);

        // people repulsion term

        /*************************************
        // update position
        *************************************/

        prefered_velocity_x = _mm256_load_pd(actual_velocity + p);
        prefered_velocity_y = _mm256_load_pd(actual_velocity + number_of_people + p);

        prefered_velocity_x = _mm256_fmadd_pd(social_force_x, timestep_vec, prefered_velocity_x);
        prefered_velocity_y = _mm256_fmadd_pd(social_force_y, timestep_vec, prefered_velocity_y);

        // norm
        prefered_velocity_2_x = _mm256_mul_pd(prefered_velocity_x, prefered_velocity_x);
        prefered_velocity_norm = _mm256_sqrt_pd(_mm256_fmadd_pd(prefered_velocity_y, prefered_velocity_y, prefered_velocity_2_x));
        prefered_velocity_norm_inv = _mm256_div_pd(one,prefered_velocity_norm);

        max_speed = _mm256_load_pd(desired_max_speed + p);

        mask = _mm256_cmp_pd(prefered_velocity_norm, max_speed, _CMP_GT_OQ);
        control_value = _mm256_blendv_pd(one, _mm256_mul_pd(max_speed, prefered_velocity_norm_inv), mask);

        prefered_velocity_x = _mm256_mul_pd(prefered_velocity_x, control_value);
        prefered_velocity_y = _mm256_mul_pd(prefered_velocity_y, control_value);

        _mm256_store_pd(actual_velocity + p, prefered_velocity_x);
        _mm256_store_pd(actual_velocity + number_of_people + p, prefered_velocity_y);

        control_value_prefered_velocity_xy_norm = _mm256_mul_pd(control_value, prefered_velocity_norm);
        control_value_prefered_velocity_xy_norm_inv = _mm256_div_pd(one,control_value_prefered_velocity_xy_norm);

        _mm256_store_pd(speed + p, control_value_prefered_velocity_xy_norm);

        _mm256_store_pd(desired_direction + p, _mm256_mul_pd(prefered_velocity_x, control_value_prefered_velocity_xy_norm_inv));
        _mm256_store_pd(desired_direction + number_of_people + p, _mm256_mul_pd(prefered_velocity_y, control_value_prefered_velocity_xy_norm_inv));

        position_x = _mm256_load_pd(position + p);
        position_y = _mm256_load_pd(position + number_of_people + p);

        _mm256_store_pd(position + p, _mm256_fmadd_pd(prefered_velocity_x, timestep_vec, position_x));
        _mm256_store_pd(position + number_of_people + p, _mm256_fmadd_pd(prefered_velocity_y, timestep_vec, position_y));
      }
    }
    CONSOLE_PRINT(("Finished iteration %d\n", (step + 1)));
  }

  CONSOLE_PRINT(("Simulation terminated\n"));
}
