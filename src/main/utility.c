/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>

#include "utility.h"
#include "social_force.h"
#include "parse_args.h"

// fix aligenment
#include "aligned_free.h"
#include "aligned_malloc.h"

extern struct arguments arguments;
extern char filename_global[80];

/*
    This function returns a sample point from a normal distribution with mean "mu" and std.deviation sqrt of "sigma"

    Assumptions: The parameters are positive
  Parameters: 
              sigma: (float) : squared value of the std.deviation of the distribution
                 mu: (float) : mean of the distribution            
*/
float sampleNormal(float sigma, float mu)
{
    float u = ((float)rand() / (RAND_MAX)) * 2 - 1;
    float v = ((float)rand() / (RAND_MAX)) * 2 - 1;
    float r = u * u + v * v;
    if (r == 0 || r > 1)
        return sampleNormal(sigma, mu);
    float c = sqrt(-2 * log(r) / r);
    return ((u * c * sigma) + mu);
}

/*
  This function outputs the initial state of the Person matrix to the given filename
  Parameters: filename: name of the .txt file
                people: array of people
                     n: number of people
            n_features: number of features per person
*/
void output_to_file_initial_state(char *filename, float *position, float *actual_speed, float *desired_direction, float *final_destination, int n, int n_timestep)
{
    FILE *fptr;

    fptr = fopen(filename, "a"); // -a option: data is appended to end of the file; If the file does not exist, it will be created.

    if (!fptr)
    {
        perror("Error reading file in output_to_file_initial_state");
    }

    output_to_file_constants(filename);

    output_to_file_persons(filename, position, actual_speed, desired_direction, final_destination, n, n_timestep);
}

void output_to_file_initial_state_double(char *filename, double *position, double *actual_speed, double *desired_direction, double *final_destination, int n, int n_timestep)
{
    FILE *fptr;

    fptr = fopen(filename, "a"); // -a option: data is appended to end of the file; If the file does not exist, it will be created.

    if (!fptr)
    {
        perror("Error reading file in output_to_file_initial_state");
    }

    output_to_file_constants(filename);

    output_to_file_persons_double(filename, position, actual_speed, desired_direction, final_destination, n, n_timestep);
}

/*
  This function outputs the defined constants to the given filename
  Parameters: filename: name of the .txt file
*/
void output_to_file_constants(char *filename)
{
    FILE *fptr;

    fptr = fopen(filename, "a"); // -a option: data is appended to end of the file; If the file does not exist, it will be created.

    if (!fptr)
    {
        perror("Error reading file in output_to_file_constants");
        return;
    }

    fprintf(fptr, "14\n"); // output number of Variables
    fprintf(fptr, "AVG_SPEED %f\n", AVG_SPEED);
    fprintf(fptr, "RELAX_TIME %f\n", RELAX_TIME);
    fprintf(fptr, "WALK_WAY_LENGTH %f\n", arguments.walkway_length);
    fprintf(fptr, "WALK_WAY_WIDTH %f\n", arguments.walkway_width);
    fprintf(fptr, "NUMBER_OF_PEOPLE %d\n", arguments.n_people);
    fprintf(fptr, "N_BORDERS %d\n", (int)N_BORDERS);
    fprintf(fptr, "TIMESTEP %f\n", TIMESTEP);
    fprintf(fptr, "N_TIMESTEP %d\n", arguments.n_timesteps);

    fprintf(fptr, "V_ALPHA_BETA %f\n", V_ALPHA_BETA);
    fprintf(fptr, "SIGMA %f\n", SIGMA);
    fprintf(fptr, "U_ALPHA_B %f\n", U_ALPHA_B);
    fprintf(fptr, "R %f\n", R);
    fprintf(fptr, "PSI %f\n", PSI);
    fprintf(fptr, "INFLUENCE %f\n", INFLUENCE);

    fclose(fptr);
}

/*
  This function output the Person matrix to the given filename (needed after initialization)
  Parameters: filename: name of the .txt file
                people: array of people
                     n: number of people
            n_features: number of features per person
*/
void output_to_file_persons(char *filename, float *position, float *speed, float *desired_direction, float *final_destination, int n, int n_timestep)
{
    FILE *fptr;

    fptr = fopen(filename, "a"); // -a option: data is appended to end of the file; If the file does not exist, it will be created.

    if (!fptr)
    {
        perror("Error reading file in output_to_file_persons");
    }

    for (int i = 0; i < n; i++)
    {
        float position_x = position[IndexX(i)];
        float position_y = position[IndexY(i, n)];
        float speed_c = speed[i];
        float desired_direction_x = desired_direction[IndexX(i)];
        float desired_direction_y = desired_direction[IndexY(i, n)];
        float final_destination_x = final_destination[IndexX(i)];
        float final_destination_y = final_destination[IndexY(i, n)];
        float vel_x = desired_direction_x * speed_c;
        float vel_y = desired_direction_y * speed_c;

        fprintf(fptr, "%f ", position_x);
        fprintf(fptr, "%f ", position_y);
        fprintf(fptr, "%f ", speed_c);
        fprintf(fptr, "%f ", desired_direction_x);
        fprintf(fptr, "%f ", desired_direction_y);
        fprintf(fptr, "%f ", final_destination_x);
        fprintf(fptr, "%f ", final_destination_y);
        fprintf(fptr, "%f ", vel_x);
        fprintf(fptr, "%f ", vel_y);
        fprintf(fptr, "\n");
    }

    fclose(fptr);
}

void output_to_file_persons_double(char *filename, double *position, double *speed, double *desired_direction, double *final_destination, int n, int n_timestep)
{
    FILE *fptr;

    fptr = fopen(filename, "a"); // -a option: data is appended to end of the file; If the file does not exist, it will be created.

    if (!fptr)
    {
        perror("Error reading file in output_to_file_persons");
    }

    for (int i = 0; i < n; i++)
    {
        double position_x = position[IndexX(i)];
        double position_y = position[IndexY(i, n)];
        double speed_c = speed[i];
        double desired_direction_x = desired_direction[IndexX(i)];
        double desired_direction_y = desired_direction[IndexY(i, n)];
        double final_destination_x = final_destination[IndexX(i)];
        double final_destination_y = final_destination[IndexY(i, n)];
        double vel_x = desired_direction_x * speed_c;
        double vel_y = desired_direction_y * speed_c;

        fprintf(fptr, "%f ", position_x);
        fprintf(fptr, "%f ", position_y);
        fprintf(fptr, "%f ", speed_c);
        fprintf(fptr, "%f ", desired_direction_x);
        fprintf(fptr, "%f ", desired_direction_y);
        fprintf(fptr, "%f ", final_destination_x);
        fprintf(fptr, "%f ", final_destination_y);
        fprintf(fptr, "%f ", vel_x);
        fprintf(fptr, "%f ", vel_y);
        fprintf(fptr, "\n");
    }

    fclose(fptr);
}
/*
  This function generates a filename with the current timestamp, accessible through filename_global
*/
void get_filename()
{

    if (arguments.filename == "\n")
    {
        time_t t = time(NULL);
        struct tm tm = *localtime(&t);
        sprintf(filename_global, "../../test/basic_%d-%02d-%02d_%02d-%02d-%02d%s", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, ".txt");
    }
    else if ('_' == arguments.filename[0])
    {
        strcpy(filename_global, "../../test/basic");
        strcat(filename_global, arguments.filename);
    }
    else
    {
        strcpy(filename_global, "../../test/");
        strcpy(filename_global, arguments.filename);
    }
}

/*
*   Function that frees all the pointers to float passed as an argument
*/
void free_all(int n, ...)
{
    int i;
    float **cur;

    va_list list;
    va_start(list, n);

    for (int i = 0; i < n; i++)
    {
        cur = va_arg(list, float **);
        aligned_free(*cur);
    }

    va_end(list);
}

int contains_substring(char *haystack, char *needle)
{
    return strstr(haystack, needle) != NULL;
}

float exp_fast(float x)
{
    x = 1.0 + x / 16384;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    return x;
}

float exp_taylor(float x)
{
    float sum = 1.0f; // initialize sum of series

    for (int i = 5; i > 0; --i)
        sum = 1 + x * sum / i;

    return sum;
}

void set_zero(float *p, int size)
{
    for (int i = 0; i < size; i++)
    {
        p[i] = 0;
    }
}

void set_zero_double(double *p, int size)
{
    for (int i = 0; i < size; i++)
    {
        p[i] = 0;
    }
}