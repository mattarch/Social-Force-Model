/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

// system includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>

double exp_fast_custom_double(double x, int b)
{
    double denom = 1;
    for (int i = 0; i < b; i++)
    {
        denom *= 2;
    }
    x = 1.0 + x / denom;
    for (int i = 0; i < b; i++)
    {
        x *= x;
    }
    return x;
}

float exp_fast_custom_float(float x, int b)
{
    float denom = 1;
    for (int i = 0; i < b; i++)
    {
        denom *= 2;
    }
    x = 1.0 + x / denom;
    for (int i = 0; i < b; i++)
    {
        x *= x;
    }
    return x;
}

__m256d exp_fast_vec_double(__m256d x, int b)
{

    __m256d denom = _mm256_set1_pd(1.0);
    __m256d two = _mm256_set1_pd(2.0);

    for (int i = 0; i < b; i++)
    {
        denom = _mm256_mul_pd(denom, two);
    }
    x = _mm256_add_pd(_mm256_set1_pd(1.0), _mm256_div_pd(x, denom));
    for (int i = 0; i < b; i++)
    {
        x = _mm256_mul_pd(x, x);
    }

    return x;
}

__m256 exp_fast_vec_float(__m256 x, int b)
{

    __m256 denom = _mm256_set1_ps(1.0);
    __m256 two = _mm256_set1_ps(2.0);

    for (int i = 0; i < b; i++)
    {
        denom = _mm256_mul_ps(denom, two);
    }
    x = _mm256_add_ps(_mm256_set1_ps(1.0), _mm256_div_ps(x, denom));
    for (int i = 0; i < b; i++)
    {
        x = _mm256_mul_ps(x, x);
    }

    return x;
}

int choose = 2;

// main function
int main(int argc, char *argv[])
{
    if (argc <= 1)
    {
        printf("You did not feed me arguments, I will die now :( ...\n");
        exit(1);
    } //otherwise continue on our merry way....
    choose = atoi(argv[1]);

    if (choose == 0)
    {
        /* case stdc floats */

        for (int j = 5; j < 20; j++)
        {
            float input = 0.0;

            for (int i = 0; i < 1000; i++)
            {
                input -= 0.05;
                float current_input = input + rand() * 0.05 / RAND_MAX;
                float exp_real = exp(current_input);
                float exp_approx = exp_fast_custom_float(current_input, j);
                float diff = exp_real - exp_approx;
                printf("exp_fast_%d %d %f %f %f\n", j, j, exp_real, exp_approx, diff);
            }
        }
    }
    else if (choose == 1)
    {
        /* case stdc doubles */

        for (int j = 5; j < 50; j++)
        {
            double input = 0.0;

            for (int i = 0; i < 1000; i++)
            {
                input -= 0.05;
                double current_input = input + rand() * 0.05 / RAND_MAX;
                double exp_real = exp(current_input);
                double exp_approx = exp_fast_custom_double(current_input, j);
                double diff = exp_real - exp_approx;
                printf("exp_fast_%d %d %f %f %f\n", j, j, exp_real, exp_approx, diff);
            }
        }
    }
    else if (choose == 2)
    {
        /* case SIMD floats */

        for (int j = 5; j < 20; j++)
        {
            float input_f = 0.0;
            __m256 input = _mm256_setzero_ps();
            __m256 constant = _mm256_set1_ps(0.05);

            for (int i = 0; i < 1000; i++)
            {
                input_f -= 0.05;
                float rand_val = rand() * 0.05 / RAND_MAX;
                float current_input_f = input_f + rand_val;

                float exp_real = exp(current_input_f);

                input = _mm256_sub_ps(input, constant);
                __m256 random_value = _mm256_set1_ps(rand_val);
                __m256 current_input = _mm256_add_ps(input, random_value);

                __m256 exp_approx = exp_fast_vec_float(current_input, j);

                float diff = exp_real - exp_approx[0];
                printf("exp_fast_%d %d %f %f %f\n", j, j, exp_real, exp_approx[0], diff);
            }
        }
    }
    else if (choose == 3)
    {

        /* case SIMD doubles */

        for (int j = 5; j < 50; j++)
        {
            double input_f = 0.0;
            __m256d input = _mm256_setzero_pd();
            __m256d constant = _mm256_set1_pd(0.05);

            for (int i = 0; i < 1000; i++)
            {
                input_f -= 0.05;
                double rand_val = rand() * 0.05 / RAND_MAX;
                double current_input_f = input_f + rand_val;

                double exp_real = exp(current_input_f);

                input = _mm256_sub_pd(input, constant);
                __m256d random_value = _mm256_set1_pd(rand_val);
                __m256d current_input = _mm256_add_pd(input, random_value);

                __m256d exp_approx = exp_fast_vec_double(current_input, j);

                double diff = exp_real - exp_approx[0];
                printf("exp_fast_%d %d %f %f %f\n", j, j, exp_real, exp_approx[0], diff);
            }
        }
    }
    return 0;
}
