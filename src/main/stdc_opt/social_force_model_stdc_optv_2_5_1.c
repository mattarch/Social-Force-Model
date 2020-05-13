/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#include <stdio.h>
#include <math.h>

#include "social_force_model_stdc_optv_2_5_1.h"
#include "../social_force.h"
#include "../utility.h"

void simulation_basic_optv_2_5_1(int number_of_people, int n_timesteps, double *position, double *speed, double *desired_direction, double *final_destination,
                                 double *borders, double *actual_velocity, double *acceleration_term, double *people_repulsion_term, double *border_repulsion_term,
                                 double *social_force, double *desired_speed, double *desired_max_speed)
{
  const double inv_sigma = 1 / SIGMA; // 1 div -> 1 flop
  const double cospsi = cos(PSI);
  int n = number_of_people;
  const double fb = borders[0]; // first border, this is positive and equal to walkway width
  const double sb = borders[1]; // second border, this is 0
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
      double sfx0 = social_force[IndexX(i)];
      double sfy0 = social_force[IndexY(i, n)];
      double sfx1 = social_force[IndexX(i + 1)];
      double sfy1 = social_force[IndexY(i + 1, n)];
      double sfx2 = social_force[IndexX(i + 2)];
      double sfy2 = social_force[IndexY(i + 2, n)];
      double sfx3 = social_force[IndexX(i + 3)];
      double sfy3 = social_force[IndexY(i + 3, n)];
      double sfx4 = social_force[IndexX(i + 4)];
      double sfy4 = social_force[IndexY(i + 4, n)];
      double sfx5 = social_force[IndexX(i + 5)];
      double sfy5 = social_force[IndexY(i + 5, n)];
      double sfx6 = social_force[IndexX(i + 6)];
      double sfy6 = social_force[IndexY(i + 6, n)];
      double sfx7 = social_force[IndexX(i + 7)];
      double sfy7 = social_force[IndexY(i + 7, n)];

      double rxa0 = position[IndexX(i)];
      double rya0 = position[IndexY(i, n)];
      double rxa1 = position[IndexX(i + 1)];
      double rya1 = position[IndexY(i + 1, n)];
      double rxa2 = position[IndexX(i + 2)];
      double rya2 = position[IndexY(i + 2, n)];
      double rxa3 = position[IndexX(i + 3)];
      double rya3 = position[IndexY(i + 3, n)];
      double rxa4 = position[IndexX(i + 4)];
      double rya4 = position[IndexY(i + 4, n)];
      double rxa5 = position[IndexX(i + 5)];
      double rya5 = position[IndexY(i + 5, n)];
      double rxa6 = position[IndexX(i + 6)];
      double rya6 = position[IndexY(i + 6, n)];
      double rxa7 = position[IndexX(i + 7)];
      double rya7 = position[IndexY(i + 7, n)];

      double exa0 = desired_direction[IndexX(i)];
      double eya0 = desired_direction[IndexY(i, n)];
      double exa1 = desired_direction[IndexX(i + 1)];
      double eya1 = desired_direction[IndexY(i + 1, n)];
      double exa2 = desired_direction[IndexX(i + 2)];
      double eya2 = desired_direction[IndexY(i + 2, n)];
      double exa3 = desired_direction[IndexX(i + 3)];
      double eya3 = desired_direction[IndexY(i + 3, n)];
      double exa4 = desired_direction[IndexX(i + 4)];
      double eya4 = desired_direction[IndexY(i + 4, n)];
      double exa5 = desired_direction[IndexX(i + 5)];
      double eya5 = desired_direction[IndexY(i + 5, n)];
      double exa6 = desired_direction[IndexX(i + 6)];
      double eya6 = desired_direction[IndexY(i + 6, n)];
      double exa7 = desired_direction[IndexX(i + 7)];
      double eya7 = desired_direction[IndexY(i + 7, n)];

      double avx0 = actual_velocity[IndexX(i)];
      double avy0 = actual_velocity[IndexY(i, n)];
      double avx1 = actual_velocity[IndexX(i + 1)];
      double avy1 = actual_velocity[IndexY(i + 1, n)];
      double avx2 = actual_velocity[IndexX(i + 2)];
      double avy2 = actual_velocity[IndexY(i + 2, n)];
      double avx3 = actual_velocity[IndexX(i + 3)];
      double avy3 = actual_velocity[IndexY(i + 3, n)];
      double avx4 = actual_velocity[IndexX(i + 4)];
      double avy4 = actual_velocity[IndexY(i + 4, n)];
      double avx5 = actual_velocity[IndexX(i + 5)];
      double avy5 = actual_velocity[IndexY(i + 5, n)];
      double avx6 = actual_velocity[IndexX(i + 6)];
      double avy6 = actual_velocity[IndexY(i + 6, n)];
      double avx7 = actual_velocity[IndexX(i + 7)];
      double avy7 = actual_velocity[IndexY(i + 7, n)];

      double dsv0 = desired_speed[i];     //desired speed value
      double dsv1 = desired_speed[i + 1]; //desired speed value
      double dsv2 = desired_speed[i + 2]; //desired speed value
      double dsv3 = desired_speed[i + 3]; //desired speed value
      double dsv4 = desired_speed[i + 4]; //desired speed value
      double dsv5 = desired_speed[i + 5]; //desired speed value
      double dsv6 = desired_speed[i + 6]; //desired speed value
      double dsv7 = desired_speed[i + 7]; //desired speed value

      double da0 = speed[i];
      double da1 = speed[i + 1];
      double da2 = speed[i + 2];
      double da3 = speed[i + 3];
      double da4 = speed[i + 4];
      double da5 = speed[i + 5];
      double da6 = speed[i + 6];
      double da7 = speed[i + 7];

      /************************************************/
      // UPDATE BORDER REPULSION TERM
      /************************************************/

      double ryaB00 = rya0 - fb; 
      double ryaB10 = rya1 - fb; 
      double ryaB20 = rya2 - fb; 
      double ryaB30 = rya3 - fb; 
      double ryaB40 = rya4 - fb; 
      double ryaB50 = rya5 - fb; 
      double ryaB60 = rya6 - fb; 
      double ryaB70 = rya7 - fb; 

      double se00 = exp_fast(ryaB00 * INV_R) * UTIMESR / (-ryaB00); //1 exp, 2 mult, 1 div
      double se10 = exp_fast(ryaB10 * INV_R) * UTIMESR / (-ryaB10); //1 exp, 2 mult, 1 div
      double se20 = exp_fast(ryaB20 * INV_R) * UTIMESR / (-ryaB20); //1 exp, 2 mult, 1 div
      double se30 = exp_fast(ryaB30 * INV_R) * UTIMESR / (-ryaB30); //1 exp, 2 mult, 1 div
      double se40 = exp_fast(ryaB40 * INV_R) * UTIMESR / (-ryaB40); //1 exp, 2 mult, 1 div
      double se50 = exp_fast(ryaB50 * INV_R) * UTIMESR / (-ryaB50); //1 exp, 2 mult, 1 div
      double se60 = exp_fast(ryaB60 * INV_R) * UTIMESR / (-ryaB60); //1 exp, 2 mult, 1 div
      double se70 = exp_fast(ryaB70 * INV_R) * UTIMESR / (-ryaB70); //1 exp, 2 mult, 1 div

      double rb00 = se00 * ryaB00; // repulsion border 0 , 1 mult
      double rb10 = se10 * ryaB10; // repulsion border 0 , 1 mult
      double rb20 = se20 * ryaB20; // repulsion border 0 , 1 mult
      double rb30 = se30 * ryaB30; // repulsion border 0 , 1 mult
      double rb40 = se40 * ryaB40; // repulsion border 0 , 1 mult
      double rb50 = se50 * ryaB50; // repulsion border 0 , 1 mult
      double rb60 = se60 * ryaB60; // repulsion border 0 , 1 mult
      double rb70 = se70 * ryaB70; // repulsion border 0 , 1 mult

      double ryaB01 = rya0 - sb; //1 add
      double ryaB11 = rya1 - sb; //1 add
      double ryaB21 = rya2 - sb; //1 add
      double ryaB31 = rya3 - sb; //1 add
      double ryaB41 = rya4 - sb; //1 add
      double ryaB51 = rya5 - sb; //1 add
      double ryaB61 = rya6 - sb; //1 add
      double ryaB71 = rya7 - sb; //1 add

      double se01 = exp_fast(-ryaB01 * INV_R) * UTIMESR / (ryaB01); //1 exp, 2 mult, 1 div
      double se11 = exp_fast(-ryaB11 * INV_R) * UTIMESR / (ryaB11); //1 exp, 2 mult, 1 div
      double se21 = exp_fast(-ryaB21 * INV_R) * UTIMESR / (ryaB21); //1 exp, 2 mult, 1 div
      double se31 = exp_fast(-ryaB31 * INV_R) * UTIMESR / (ryaB31); //1 exp, 2 mult, 1 div
      double se41 = exp_fast(-ryaB41 * INV_R) * UTIMESR / (ryaB41); //1 exp, 2 mult, 1 div
      double se51 = exp_fast(-ryaB51 * INV_R) * UTIMESR / (ryaB51); //1 exp, 2 mult, 1 div
      double se61 = exp_fast(-ryaB61 * INV_R) * UTIMESR / (ryaB61); //1 exp, 2 mult, 1 div
      double se71 = exp_fast(-ryaB71 * INV_R) * UTIMESR / (ryaB71); //1 exp, 2 mult, 1 div

      double rb01 = rb00 + se01 * ryaB01; //1 mult
      double rb11 = rb10 + se11 * ryaB11; //1 mult
      double rb21 = rb20 + se21 * ryaB21; //1 mult
      double rb31 = rb30 + se31 * ryaB31; //1 mult
      double rb41 = rb40 + se41 * ryaB41; //1 mult
      double rb51 = rb50 + se51 * ryaB51; //1 mult
      double rb61 = rb60 + se61 * ryaB61; //1 mult
      double rb71 = rb70 + se71 * ryaB71; //1 mult

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

      double rxa0a1 = rxa0 - rxa1;
      double rya0a1 = rya0 - rya1;
      double rxa0a2 = rxa0 - rxa2;
      double rya0a2 = rya0 - rya2;
      double rxa0a3 = rxa0 - rxa3;
      double rya0a3 = rya0 - rya3;
      double rxa0a4 = rxa0 - rxa4;
      double rya0a4 = rya0 - rya4;
      double rxa0a5 = rxa0 - rxa5;
      double rya0a5 = rya0 - rya5;
      double rxa0a6 = rxa0 - rxa6;
      double rya0a6 = rya0 - rya6;
      double rxa0a7 = rxa0 - rxa7;
      double rya0a7 = rya0 - rya7;

      double rxa1a2 = rxa1 - rxa2;
      double rya1a2 = rya1 - rya2;
      double rxa1a3 = rxa1 - rxa3;
      double rya1a3 = rya1 - rya3;
      double rxa1a4 = rxa1 - rxa4;
      double rya1a4 = rya1 - rya4;
      double rxa1a5 = rxa1 - rxa5;
      double rya1a5 = rya1 - rya5;
      double rxa1a6 = rxa1 - rxa6;
      double rya1a6 = rya1 - rya6;
      double rxa1a7 = rxa1 - rxa7;
      double rya1a7 = rya1 - rya7;

      double rxa2a3 = rxa2 - rxa3;
      double rya2a3 = rya2 - rya3;
      double rxa2a4 = rxa2 - rxa4;
      double rya2a4 = rya2 - rya4;
      double rxa2a5 = rxa2 - rxa5;
      double rya2a5 = rya2 - rya5;
      double rxa2a6 = rxa2 - rxa6;
      double rya2a6 = rya2 - rya6;
      double rxa2a7 = rxa2 - rxa7;
      double rya2a7 = rya2 - rya7;

      double rxa3a4 = rxa3 - rxa4;
      double rya3a4 = rya3 - rya4;
      double rxa3a5 = rxa3 - rxa5;
      double rya3a5 = rya3 - rya5;
      double rxa3a6 = rxa3 - rxa6;
      double rya3a6 = rya3 - rya6;
      double rxa3a7 = rxa3 - rxa7;
      double rya3a7 = rya3 - rya7;

      double rxa4a5 = rxa4 - rxa5;
      double rya4a5 = rya4 - rya5;
      double rxa4a6 = rxa4 - rxa6;
      double rya4a6 = rya4 - rya6;
      double rxa4a7 = rxa4 - rxa7;
      double rya4a7 = rya4 - rya7;

      double rxa5a6 = rxa5 - rxa6;
      double rya5a6 = rya5 - rya6;
      double rxa5a7 = rxa5 - rxa7;
      double rya5a7 = rya5 - rya7;

      double rxa6a7 = rxa6 - rxa7;
      double rya6a7 = rya6 - rya7;

      double ra0a1 = rxa0a1 * rxa0a1 + rya0a1 * rya0a1;
      double ra0a2 = rxa0a2 * rxa0a2 + rya0a2 * rya0a2;
      double ra0a3 = rxa0a3 * rxa0a3 + rya0a3 * rya0a3;
      double ra0a4 = rxa0a4 * rxa0a4 + rya0a4 * rya0a4;
      double ra0a5 = rxa0a5 * rxa0a5 + rya0a5 * rya0a5;
      double ra0a6 = rxa0a6 * rxa0a6 + rya0a6 * rya0a6;
      double ra0a7 = rxa0a7 * rxa0a7 + rya0a7 * rya0a7;

      double ra1a2 = rxa1a2 * rxa1a2 + rya1a2 * rya1a2;
      double ra1a3 = rxa1a3 * rxa1a3 + rya1a3 * rya1a3;
      double ra1a4 = rxa1a4 * rxa1a4 + rya1a4 * rya1a4;
      double ra1a5 = rxa1a5 * rxa1a5 + rya1a5 * rya1a5;
      double ra1a6 = rxa1a6 * rxa1a6 + rya1a6 * rya1a6;
      double ra1a7 = rxa1a7 * rxa1a7 + rya1a7 * rya1a7;

      double ra2a3 = rxa2a3 * rxa2a3 + rya2a3 * rya2a3;
      double ra2a4 = rxa2a4 * rxa2a4 + rya2a4 * rya2a4;
      double ra2a5 = rxa2a5 * rxa2a5 + rya2a5 * rya2a5;
      double ra2a6 = rxa2a6 * rxa2a6 + rya2a6 * rya2a6;
      double ra2a7 = rxa2a7 * rxa2a7 + rya2a7 * rya2a7;

      double ra3a4 = rxa3a4 * rxa3a4 + rya3a4 * rya3a4;
      double ra3a5 = rxa3a5 * rxa3a5 + rya3a5 * rya3a5;
      double ra3a6 = rxa3a6 * rxa3a6 + rya3a6 * rya3a6;
      double ra3a7 = rxa3a7 * rxa3a7 + rya3a7 * rya3a7;

      double ra4a5 = rxa4a5 * rxa4a5 + rya4a5 * rya4a5;
      double ra4a6 = rxa4a6 * rxa4a6 + rya4a6 * rya4a6;
      double ra4a7 = rxa4a7 * rxa4a7 + rya4a7 * rya4a7;

      double ra5a6 = rxa5a6 * rxa5a6 + rya5a6 * rya5a6;
      double ra5a7 = rxa5a7 * rxa5a7 + rya5a7 * rya5a7;

      double ra6a7 = rxa6a7 * rxa6a7 + rya6a7 * rya6a7;

      double ra0a1norm = sqrt(ra0a1);
      double ra0a2norm = sqrt(ra0a2);
      double ra0a3norm = sqrt(ra0a3);
      double ra0a4norm = sqrt(ra0a4);
      double ra0a5norm = sqrt(ra0a5);
      double ra0a6norm = sqrt(ra0a6);
      double ra0a7norm = sqrt(ra0a7);

      double ra1a2norm = sqrt(ra1a2);
      double ra1a3norm = sqrt(ra1a3);
      double ra1a4norm = sqrt(ra1a4);
      double ra1a5norm = sqrt(ra1a5);
      double ra1a6norm = sqrt(ra1a6);
      double ra1a7norm = sqrt(ra1a7);

      double ra2a3norm = sqrt(ra2a3);
      double ra2a4norm = sqrt(ra2a4);
      double ra2a5norm = sqrt(ra2a5);
      double ra2a6norm = sqrt(ra2a6);
      double ra2a7norm = sqrt(ra2a7);

      double ra3a4norm = sqrt(ra3a4);
      double ra3a5norm = sqrt(ra3a5);
      double ra3a6norm = sqrt(ra3a6);
      double ra3a7norm = sqrt(ra3a7);

      double ra4a5norm = sqrt(ra4a5);
      double ra4a6norm = sqrt(ra4a6);
      double ra4a7norm = sqrt(ra4a7);

      double ra5a6norm = sqrt(ra5a6);
      double ra5a7norm = sqrt(ra5a7);

      double ra6a7norm = sqrt(ra6a7);

      //everything twice

      double rxa0a1mex = rxa0a1 - da1 * exa1; //1 add, 1 mult
      double rya0a1mey = rya0a1 - da1 * eya1; //1 add, 1 mult
      double rxa0a2mex = rxa0a2 - da2 * exa2; //1 add, 1 mult
      double rya0a2mey = rya0a2 - da2 * eya2; //1 add, 1 mult
      double rxa0a3mex = rxa0a3 - da3 * exa3; //1 add, 1 mult
      double rya0a3mey = rya0a3 - da3 * eya3; //1 add, 1 mult
      double rxa0a4mex = rxa0a4 - da4 * exa4; //1 add, 1 mult
      double rya0a4mey = rya0a4 - da4 * eya4; //1 add, 1 mult
      double rxa0a5mex = rxa0a5 - da5 * exa5; //1 add, 1 mult
      double rya0a5mey = rya0a5 - da5 * eya5; //1 add, 1 mult
      double rxa0a6mex = rxa0a6 - da6 * exa6; //1 add, 1 mult
      double rya0a6mey = rya0a6 - da6 * eya6; //1 add, 1 mult
      double rxa0a7mex = rxa0a7 - da7 * exa7; //1 add, 1 mult
      double rya0a7mey = rya0a7 - da7 * eya7; //1 add, 1 mult

      double rxa1a0mex = -rxa0a1 - da0 * exa0; //1 add, 1 mult
      double rya1a0mey = -rya0a1 - da0 * eya0; //1 add, 1 mult
      double rxa2a0mex = -rxa0a2 - da0 * exa0; //1 add, 1 mult
      double rya2a0mey = -rya0a2 - da0 * eya0; //1 add, 1 mult
      double rxa3a0mex = -rxa0a3 - da0 * exa0; //1 add, 1 mult
      double rya3a0mey = -rya0a3 - da0 * eya0; //1 add, 1 mult
      double rxa4a0mex = -rxa0a4 - da0 * exa0; //1 add, 1 mult
      double rya4a0mey = -rya0a4 - da0 * eya0; //1 add, 1 mult
      double rxa5a0mex = -rxa0a5 - da0 * exa0; //1 add, 1 mult
      double rya5a0mey = -rya0a5 - da0 * eya0; //1 add, 1 mult
      double rxa6a0mex = -rxa0a6 - da0 * exa0; //1 add, 1 mult
      double rya6a0mey = -rya0a6 - da0 * eya0; //1 add, 1 mult
      double rxa7a0mex = -rxa0a7 - da0 * exa0; //1 add, 1 mult
      double rya7a0mey = -rya0a7 - da0 * eya0; //1 add, 1 mult

      double rxa1a2mex = rxa1a2 - da2 * exa2; //1 add, 1 mult
      double rya1a2mey = rya1a2 - da2 * eya2; //1 add, 1 mult
      double rxa1a3mex = rxa1a3 - da3 * exa3; //1 add, 1 mult
      double rya1a3mey = rya1a3 - da3 * eya3; //1 add, 1 mult
      double rxa1a4mex = rxa1a4 - da4 * exa4; //1 add, 1 mult
      double rya1a4mey = rya1a4 - da4 * eya4; //1 add, 1 mult
      double rxa1a5mex = rxa1a5 - da5 * exa5; //1 add, 1 mult
      double rya1a5mey = rya1a5 - da5 * eya5; //1 add, 1 mult
      double rxa1a6mex = rxa1a6 - da6 * exa6; //1 add, 1 mult
      double rya1a6mey = rya1a6 - da6 * eya6; //1 add, 1 mult
      double rxa1a7mex = rxa1a7 - da7 * exa7; //1 add, 1 mult
      double rya1a7mey = rya1a7 - da7 * eya7; //1 add, 1 mult

      double rxa2a1mex = -rxa1a2 - da1 * exa1; //1 add, 1 mult
      double rya2a1mey = -rya1a2 - da1 * eya1; //1 add, 1 mult
      double rxa3a1mex = -rxa1a3 - da1 * exa1; //1 add, 1 mult
      double rya3a1mey = -rya1a3 - da1 * eya1; //1 add, 1 mult
      double rxa4a1mex = -rxa1a4 - da1 * exa1; //1 add, 1 mult
      double rya4a1mey = -rya1a4 - da1 * eya1; //1 add, 1 mult
      double rxa5a1mex = -rxa1a5 - da1 * exa1; //1 add, 1 mult
      double rya5a1mey = -rya1a5 - da1 * eya1; //1 add, 1 mult
      double rxa6a1mex = -rxa1a6 - da1 * exa1; //1 add, 1 mult
      double rya6a1mey = -rya1a6 - da1 * eya1; //1 add, 1 mult
      double rxa7a1mex = -rxa1a7 - da1 * exa1; //1 add, 1 mult
      double rya7a1mey = -rya1a7 - da1 * eya1; //1 add, 1 mult

      double rxa2a3mex = rxa2a3 - da3 * exa3; //1 add, 1 mult
      double rya2a3mey = rya2a3 - da3 * eya3; //1 add, 1 mult
      double rxa2a4mex = rxa2a4 - da4 * exa4; //1 add, 1 mult
      double rya2a4mey = rya2a4 - da4 * eya4; //1 add, 1 mult
      double rxa2a5mex = rxa2a5 - da5 * exa5; //1 add, 1 mult
      double rya2a5mey = rya2a5 - da5 * eya5; //1 add, 1 mult
      double rxa2a6mex = rxa2a6 - da6 * exa6; //1 add, 1 mult
      double rya2a6mey = rya2a6 - da6 * eya6; //1 add, 1 mult
      double rxa2a7mex = rxa2a7 - da7 * exa7; //1 add, 1 mult
      double rya2a7mey = rya2a7 - da7 * eya7; //1 add, 1 mult

      double rxa3a2mex = -rxa2a3 - da2 * exa2; //1 add, 1 mult
      double rya3a2mey = -rya2a3 - da2 * eya2; //1 add, 1 mult
      double rxa4a2mex = -rxa2a4 - da2 * exa2; //1 add, 1 mult
      double rya4a2mey = -rya2a4 - da2 * eya2; //1 add, 1 mult
      double rxa5a2mex = -rxa2a5 - da2 * exa2; //1 add, 1 mult
      double rya5a2mey = -rya2a5 - da2 * eya2; //1 add, 1 mult
      double rxa6a2mex = -rxa2a6 - da2 * exa2; //1 add, 1 mult
      double rya6a2mey = -rya2a6 - da2 * eya2; //1 add, 1 mult
      double rxa7a2mex = -rxa2a7 - da2 * exa2; //1 add, 1 mult
      double rya7a2mey = -rya2a7 - da2 * eya2; //1 add, 1 mult

      double rxa3a4mex = rxa3a4 - da4 * exa4; //1 add, 1 mult
      double rya3a4mey = rya3a4 - da4 * eya4; //1 add, 1 mult
      double rxa3a5mex = rxa3a5 - da5 * exa5; //1 add, 1 mult
      double rya3a5mey = rya3a5 - da5 * eya5; //1 add, 1 mult
      double rxa3a6mex = rxa3a6 - da6 * exa6; //1 add, 1 mult
      double rya3a6mey = rya3a6 - da6 * eya6; //1 add, 1 mult
      double rxa3a7mex = rxa3a7 - da7 * exa7; //1 add, 1 mult
      double rya3a7mey = rya3a7 - da7 * eya7; //1 add, 1 mult

      double rxa4a3mex = -rxa3a4 - da3 * exa3; //1 add, 1 mult
      double rya4a3mey = -rya3a4 - da3 * eya3; //1 add, 1 mult
      double rxa5a3mex = -rxa3a5 - da3 * exa3; //1 add, 1 mult
      double rya5a3mey = -rya3a5 - da3 * eya3; //1 add, 1 mult
      double rxa6a3mex = -rxa3a6 - da3 * exa3; //1 add, 1 mult
      double rya6a3mey = -rya3a6 - da3 * eya3; //1 add, 1 mult
      double rxa7a3mex = -rxa3a7 - da3 * exa3; //1 add, 1 mult
      double rya7a3mey = -rya3a7 - da3 * eya3; //1 add, 1 mult

      double rxa4a5mex = rxa4a5 - da5 * exa5; //1 add, 1 mult
      double rya4a5mey = rya4a5 - da5 * eya5; //1 add, 1 mult
      double rxa4a6mex = rxa4a6 - da6 * exa6; //1 add, 1 mult
      double rya4a6mey = rya4a6 - da6 * eya6; //1 add, 1 mult
      double rxa4a7mex = rxa4a7 - da7 * exa7; //1 add, 1 mult
      double rya4a7mey = rya4a7 - da7 * eya7; //1 add, 1 mult

      double rxa5a4mex = -rxa4a5 - da4 * exa4; //1 add, 1 mult
      double rya5a4mey = -rya4a5 - da4 * eya4; //1 add, 1 mult
      double rxa6a4mex = -rxa4a6 - da4 * exa4; //1 add, 1 mult
      double rya6a4mey = -rya4a6 - da4 * eya4; //1 add, 1 mult
      double rxa7a4mex = -rxa4a7 - da4 * exa4; //1 add, 1 mult
      double rya7a4mey = -rya4a7 - da4 * eya4; //1 add, 1 mult

      double rxa5a6mex = rxa5a6 - da6 * exa6; //1 add, 1 mult
      double rya5a6mey = rya5a6 - da6 * eya6; //1 add, 1 mult
      double rxa5a7mex = rxa5a7 - da7 * exa7; //1 add, 1 mult
      double rya5a7mey = rya5a7 - da7 * eya7; //1 add, 1 mult

      double rxa6a5mex = -rxa5a6 - da5 * exa5; //1 add, 1 mult
      double rya6a5mey = -rya5a6 - da5 * eya5; //1 add, 1 mult
      double rxa7a5mex = -rxa5a7 - da5 * exa5; //1 add, 1 mult
      double rya7a5mey = -rya5a7 - da5 * eya5; //1 add, 1 mult

      double rxa6a7mex = rxa6a7 - da7 * exa7; //1 add, 1 mult
      double rya6a7mey = rya6a7 - da7 * eya7; //1 add, 1 mult

      double rxa7a6mex = -rxa6a7 - da6 * exa6; //1 add, 1 mult
      double rya7a6mey = -rya6a7 - da6 * eya6; //1 add, 1 mult

      double ra0a1me = rxa0a1mex * rxa0a1mex + rya0a1mey * rya0a1mey;
      double ra0a2me = rxa0a2mex * rxa0a2mex + rya0a2mey * rya0a2mey;
      double ra0a3me = rxa0a3mex * rxa0a3mex + rya0a3mey * rya0a3mey;
      double ra0a4me = rxa0a4mex * rxa0a4mex + rya0a4mey * rya0a4mey;
      double ra0a5me = rxa0a5mex * rxa0a5mex + rya0a5mey * rya0a5mey;
      double ra0a6me = rxa0a6mex * rxa0a6mex + rya0a6mey * rya0a6mey;
      double ra0a7me = rxa0a7mex * rxa0a7mex + rya0a7mey * rya0a7mey;

      double ra1a0me = rxa1a0mex * rxa1a0mex + rya1a0mey * rya1a0mey;
      double ra2a0me = rxa2a0mex * rxa2a0mex + rya2a0mey * rya2a0mey;
      double ra3a0me = rxa3a0mex * rxa3a0mex + rya3a0mey * rya3a0mey;
      double ra4a0me = rxa4a0mex * rxa4a0mex + rya4a0mey * rya4a0mey;
      double ra5a0me = rxa5a0mex * rxa5a0mex + rya5a0mey * rya5a0mey;
      double ra6a0me = rxa6a0mex * rxa6a0mex + rya6a0mey * rya6a0mey;
      double ra7a0me = rxa7a0mex * rxa7a0mex + rya7a0mey * rya7a0mey;

      double ra1a2me = rxa1a2mex * rxa1a2mex + rya1a2mey * rya1a2mey;
      double ra1a3me = rxa1a3mex * rxa1a3mex + rya1a3mey * rya1a3mey;
      double ra1a4me = rxa1a4mex * rxa1a4mex + rya1a4mey * rya1a4mey;
      double ra1a5me = rxa1a5mex * rxa1a5mex + rya1a5mey * rya1a5mey;
      double ra1a6me = rxa1a6mex * rxa1a6mex + rya1a6mey * rya1a6mey;
      double ra1a7me = rxa1a7mex * rxa1a7mex + rya1a7mey * rya1a7mey;

      double ra2a1me = rxa2a1mex * rxa2a1mex + rya2a1mey * rya2a1mey;
      double ra3a1me = rxa3a1mex * rxa3a1mex + rya3a1mey * rya3a1mey;
      double ra4a1me = rxa4a1mex * rxa4a1mex + rya4a1mey * rya4a1mey;
      double ra5a1me = rxa5a1mex * rxa5a1mex + rya5a1mey * rya5a1mey;
      double ra6a1me = rxa6a1mex * rxa6a1mex + rya6a1mey * rya6a1mey;
      double ra7a1me = rxa7a1mex * rxa7a1mex + rya7a1mey * rya7a1mey;

      double ra2a3me = rxa2a3mex * rxa2a3mex + rya2a3mey * rya2a3mey;
      double ra2a4me = rxa2a4mex * rxa2a4mex + rya2a4mey * rya2a4mey;
      double ra2a5me = rxa2a5mex * rxa2a5mex + rya2a5mey * rya2a5mey;
      double ra2a6me = rxa2a6mex * rxa2a6mex + rya2a6mey * rya2a6mey;
      double ra2a7me = rxa2a7mex * rxa2a7mex + rya2a7mey * rya2a7mey;

      double ra3a2me = rxa3a2mex * rxa3a2mex + rya3a2mey * rya3a2mey;
      double ra4a2me = rxa4a2mex * rxa4a2mex + rya4a2mey * rya4a2mey;
      double ra5a2me = rxa5a2mex * rxa5a2mex + rya5a2mey * rya5a2mey;
      double ra6a2me = rxa6a2mex * rxa6a2mex + rya6a2mey * rya6a2mey;
      double ra7a2me = rxa7a2mex * rxa7a2mex + rya7a2mey * rya7a2mey;

      double ra3a4me = rxa3a4mex * rxa3a4mex + rya3a4mey * rya3a4mey;
      double ra3a5me = rxa3a5mex * rxa3a5mex + rya3a5mey * rya3a5mey;
      double ra3a6me = rxa3a6mex * rxa3a6mex + rya3a6mey * rya3a6mey;
      double ra3a7me = rxa3a7mex * rxa3a7mex + rya3a7mey * rya3a7mey;

      double ra4a3me = rxa4a3mex * rxa4a3mex + rya4a3mey * rya4a3mey;
      double ra5a3me = rxa5a3mex * rxa5a3mex + rya5a3mey * rya5a3mey;
      double ra6a3me = rxa6a3mex * rxa6a3mex + rya6a3mey * rya6a3mey;
      double ra7a3me = rxa7a3mex * rxa7a3mex + rya7a3mey * rya7a3mey;

      double ra4a5me = rxa4a5mex * rxa4a5mex + rya4a5mey * rya4a5mey;
      double ra4a6me = rxa4a6mex * rxa4a6mex + rya4a6mey * rya4a6mey;
      double ra4a7me = rxa4a7mex * rxa4a7mex + rya4a7mey * rya4a7mey;

      double ra5a4me = rxa5a4mex * rxa5a4mex + rya5a4mey * rya5a4mey;
      double ra6a4me = rxa6a4mex * rxa6a4mex + rya6a4mey * rya6a4mey;
      double ra7a4me = rxa7a4mex * rxa7a4mex + rya7a4mey * rya7a4mey;

      double ra5a6me = rxa5a6mex * rxa5a6mex + rya5a6mey * rya5a6mey;
      double ra5a7me = rxa5a7mex * rxa5a7mex + rya5a7mey * rya5a7mey;

      double ra6a5me = rxa6a5mex * rxa6a5mex + rya6a5mey * rya6a5mey;
      double ra7a5me = rxa7a5mex * rxa7a5mex + rya7a5mey * rya7a5mey;

      double ra6a7me = rxa6a7mex * rxa6a7mex + rya6a7mey * rya6a7mey;

      double ra7a6me = rxa7a6mex * rxa7a6mex + rya7a6mey * rya7a6mey;

      double ra0a1menorm = sqrt(ra0a1me);
      double ra0a2menorm = sqrt(ra0a2me);
      double ra0a3menorm = sqrt(ra0a3me);
      double ra0a4menorm = sqrt(ra0a4me);
      double ra0a5menorm = sqrt(ra0a5me);
      double ra0a6menorm = sqrt(ra0a6me);
      double ra0a7menorm = sqrt(ra0a7me);

      double ra1a0menorm = sqrt(ra1a0me);
      double ra2a0menorm = sqrt(ra2a0me);
      double ra3a0menorm = sqrt(ra3a0me);
      double ra4a0menorm = sqrt(ra4a0me);
      double ra5a0menorm = sqrt(ra5a0me);
      double ra6a0menorm = sqrt(ra6a0me);
      double ra7a0menorm = sqrt(ra7a0me);

      double ra1a2menorm = sqrt(ra1a2me);
      double ra1a3menorm = sqrt(ra1a3me);
      double ra1a4menorm = sqrt(ra1a4me);
      double ra1a5menorm = sqrt(ra1a5me);
      double ra1a6menorm = sqrt(ra1a6me);
      double ra1a7menorm = sqrt(ra1a7me);

      double ra2a1menorm = sqrt(ra2a1me);
      double ra3a1menorm = sqrt(ra3a1me);
      double ra4a1menorm = sqrt(ra4a1me);
      double ra5a1menorm = sqrt(ra5a1me);
      double ra6a1menorm = sqrt(ra6a1me);
      double ra7a1menorm = sqrt(ra7a1me);

      double ra2a3menorm = sqrt(ra2a3me);
      double ra2a4menorm = sqrt(ra2a4me);
      double ra2a5menorm = sqrt(ra2a5me);
      double ra2a6menorm = sqrt(ra2a6me);
      double ra2a7menorm = sqrt(ra2a7me);

      double ra3a2menorm = sqrt(ra3a2me);
      double ra4a2menorm = sqrt(ra4a2me);
      double ra5a2menorm = sqrt(ra5a2me);
      double ra6a2menorm = sqrt(ra6a2me);
      double ra7a2menorm = sqrt(ra7a2me);

      double ra3a4menorm = sqrt(ra3a4me);
      double ra3a5menorm = sqrt(ra3a5me);
      double ra3a6menorm = sqrt(ra3a6me);
      double ra3a7menorm = sqrt(ra3a7me);

      double ra4a3menorm = sqrt(ra4a3me);
      double ra5a3menorm = sqrt(ra5a3me);
      double ra6a3menorm = sqrt(ra6a3me);
      double ra7a3menorm = sqrt(ra7a3me);

      double ra4a5menorm = sqrt(ra4a5me);
      double ra4a6menorm = sqrt(ra4a6me);
      double ra4a7menorm = sqrt(ra4a7me);

      double ra5a4menorm = sqrt(ra5a4me);
      double ra6a4menorm = sqrt(ra6a4me);
      double ra7a4menorm = sqrt(ra7a4me);

      double ra5a6menorm = sqrt(ra5a6me);
      double ra5a7menorm = sqrt(ra5a7me);

      double ra6a5menorm = sqrt(ra6a5me);
      double ra7a5menorm = sqrt(ra7a5me);

      double ra6a7menorm = sqrt(ra6a7me);

      double ra7a6menorm = sqrt(ra7a6me);

      double normsuma0a1 = ra0a1norm + ra0a1menorm;
      double normsuma0a2 = ra0a2norm + ra0a2menorm;
      double normsuma0a3 = ra0a3norm + ra0a3menorm;
      double normsuma0a4 = ra0a4norm + ra0a4menorm;
      double normsuma0a5 = ra0a5norm + ra0a5menorm;
      double normsuma0a6 = ra0a6norm + ra0a6menorm;
      double normsuma0a7 = ra0a7norm + ra0a7menorm;

      double normsuma1a0 = ra0a1norm + ra1a0menorm;
      double normsuma2a0 = ra0a2norm + ra2a0menorm;
      double normsuma3a0 = ra0a3norm + ra3a0menorm;
      double normsuma4a0 = ra0a4norm + ra4a0menorm;
      double normsuma5a0 = ra0a5norm + ra5a0menorm;
      double normsuma6a0 = ra0a6norm + ra6a0menorm;
      double normsuma7a0 = ra0a7norm + ra7a0menorm;

      double normsuma1a2 = ra1a2norm + ra1a2menorm;
      double normsuma1a3 = ra1a3norm + ra1a3menorm;
      double normsuma1a4 = ra1a4norm + ra1a4menorm;
      double normsuma1a5 = ra1a5norm + ra1a5menorm;
      double normsuma1a6 = ra1a6norm + ra1a6menorm;
      double normsuma1a7 = ra1a7norm + ra1a7menorm;

      double normsuma2a1 = ra1a2norm + ra2a1menorm;
      double normsuma3a1 = ra1a3norm + ra3a1menorm;
      double normsuma4a1 = ra1a4norm + ra4a1menorm;
      double normsuma5a1 = ra1a5norm + ra5a1menorm;
      double normsuma6a1 = ra1a6norm + ra6a1menorm;
      double normsuma7a1 = ra1a7norm + ra7a1menorm;

      double normsuma2a3 = ra2a3norm + ra2a3menorm;
      double normsuma2a4 = ra2a4norm + ra2a4menorm;
      double normsuma2a5 = ra2a5norm + ra2a5menorm;
      double normsuma2a6 = ra2a6norm + ra2a6menorm;
      double normsuma2a7 = ra2a7norm + ra2a7menorm;

      double normsuma3a2 = ra2a3norm + ra3a2menorm;
      double normsuma4a2 = ra2a4norm + ra4a2menorm;
      double normsuma5a2 = ra2a5norm + ra5a2menorm;
      double normsuma6a2 = ra2a6norm + ra6a2menorm;
      double normsuma7a2 = ra2a7norm + ra7a2menorm;

      double normsuma3a4 = ra3a4norm + ra3a4menorm;
      double normsuma3a5 = ra3a5norm + ra3a5menorm;
      double normsuma3a6 = ra3a6norm + ra3a6menorm;
      double normsuma3a7 = ra3a7norm + ra3a7menorm;

      double normsuma4a3 = ra3a4norm + ra4a3menorm;
      double normsuma5a3 = ra3a5norm + ra5a3menorm;
      double normsuma6a3 = ra3a6norm + ra6a3menorm;
      double normsuma7a3 = ra3a7norm + ra7a3menorm;

      double normsuma4a5 = ra4a5norm + ra4a5menorm;
      double normsuma4a6 = ra4a6norm + ra4a6menorm;
      double normsuma4a7 = ra4a7norm + ra4a7menorm;

      double normsuma5a4 = ra4a5norm + ra5a4menorm;
      double normsuma6a4 = ra4a6norm + ra6a4menorm;
      double normsuma7a4 = ra4a7norm + ra7a4menorm;

      double normsuma5a6 = ra5a6norm + ra5a6menorm;
      double normsuma5a7 = ra5a7norm + ra5a7menorm;

      double normsuma6a5 = ra5a6norm + ra6a5menorm;
      double normsuma7a5 = ra5a7norm + ra7a5menorm;

      double normsuma6a7 = ra6a7norm + ra6a7menorm;

      double normsuma7a6 = ra6a7norm + ra7a6menorm;

      double repxa0a1 = rxa0a1 / ra0a1norm + rxa0a1mex / ra0a1menorm;
      double repya0a1 = rya0a1 / ra0a1norm + rya0a1mey / ra0a1menorm;
      double repxa0a2 = rxa0a2 / ra0a2norm + rxa0a2mex / ra0a2menorm;
      double repya0a2 = rya0a2 / ra0a2norm + rya0a2mey / ra0a2menorm;
      double repxa0a3 = rxa0a3 / ra0a3norm + rxa0a3mex / ra0a3menorm;
      double repya0a3 = rya0a3 / ra0a3norm + rya0a3mey / ra0a3menorm;
      double repxa0a4 = rxa0a4 / ra0a4norm + rxa0a4mex / ra0a4menorm;
      double repya0a4 = rya0a4 / ra0a4norm + rya0a4mey / ra0a4menorm;
      double repxa0a5 = rxa0a5 / ra0a5norm + rxa0a5mex / ra0a5menorm;
      double repya0a5 = rya0a5 / ra0a5norm + rya0a5mey / ra0a5menorm;
      double repxa0a6 = rxa0a6 / ra0a6norm + rxa0a6mex / ra0a6menorm;
      double repya0a6 = rya0a6 / ra0a6norm + rya0a6mey / ra0a6menorm;
      double repxa0a7 = rxa0a7 / ra0a7norm + rxa0a7mex / ra0a7menorm;
      double repya0a7 = rya0a7 / ra0a7norm + rya0a7mey / ra0a7menorm;

      double repxa1a0 = -rxa0a1 / ra0a1norm + rxa1a0mex / ra1a0menorm;
      double repya1a0 = -rya0a1 / ra0a1norm + rya1a0mey / ra1a0menorm;
      double repxa2a0 = -rxa0a2 / ra0a2norm + rxa2a0mex / ra2a0menorm;
      double repya2a0 = -rya0a2 / ra0a2norm + rya2a0mey / ra2a0menorm;
      double repxa3a0 = -rxa0a3 / ra0a3norm + rxa3a0mex / ra3a0menorm;
      double repya3a0 = -rya0a3 / ra0a3norm + rya3a0mey / ra3a0menorm;
      double repxa4a0 = -rxa0a4 / ra0a4norm + rxa4a0mex / ra4a0menorm;
      double repya4a0 = -rya0a4 / ra0a4norm + rya4a0mey / ra4a0menorm;
      double repxa5a0 = -rxa0a5 / ra0a5norm + rxa5a0mex / ra5a0menorm;
      double repya5a0 = -rya0a5 / ra0a5norm + rya5a0mey / ra5a0menorm;
      double repxa6a0 = -rxa0a6 / ra0a6norm + rxa6a0mex / ra6a0menorm;
      double repya6a0 = -rya0a6 / ra0a6norm + rya6a0mey / ra6a0menorm;
      double repxa7a0 = -rxa0a7 / ra0a7norm + rxa7a0mex / ra7a0menorm;
      double repya7a0 = -rya0a7 / ra0a7norm + rya7a0mey / ra7a0menorm;

      double repxa1a2 = rxa1a2 / ra1a2norm + rxa1a2mex / ra1a2menorm;
      double repya1a2 = rya1a2 / ra1a2norm + rya1a2mey / ra1a2menorm;
      double repxa1a3 = rxa1a3 / ra1a3norm + rxa1a3mex / ra1a3menorm;
      double repya1a3 = rya1a3 / ra1a3norm + rya1a3mey / ra1a3menorm;
      double repxa1a4 = rxa1a4 / ra1a4norm + rxa1a4mex / ra1a4menorm;
      double repya1a4 = rya1a4 / ra1a4norm + rya1a4mey / ra1a4menorm;
      double repxa1a5 = rxa1a5 / ra1a5norm + rxa1a5mex / ra1a5menorm;
      double repya1a5 = rya1a5 / ra1a5norm + rya1a5mey / ra1a5menorm;
      double repxa1a6 = rxa1a6 / ra1a6norm + rxa1a6mex / ra1a6menorm;
      double repya1a6 = rya1a6 / ra1a6norm + rya1a6mey / ra1a6menorm;
      double repxa1a7 = rxa1a7 / ra1a7norm + rxa1a7mex / ra1a7menorm;
      double repya1a7 = rya1a7 / ra1a7norm + rya1a7mey / ra1a7menorm;

      double repxa2a1 = -rxa1a2 / ra1a2norm + rxa2a1mex / ra2a1menorm;
      double repya2a1 = -rya1a2 / ra1a2norm + rya2a1mey / ra2a1menorm;
      double repxa3a1 = -rxa1a3 / ra1a3norm + rxa3a1mex / ra3a1menorm;
      double repya3a1 = -rya1a3 / ra1a3norm + rya3a1mey / ra3a1menorm;
      double repxa4a1 = -rxa1a4 / ra1a4norm + rxa4a1mex / ra4a1menorm;
      double repya4a1 = -rya1a4 / ra1a4norm + rya4a1mey / ra4a1menorm;
      double repxa5a1 = -rxa1a5 / ra1a5norm + rxa5a1mex / ra5a1menorm;
      double repya5a1 = -rya1a5 / ra1a5norm + rya5a1mey / ra5a1menorm;
      double repxa6a1 = -rxa1a6 / ra1a6norm + rxa6a1mex / ra6a1menorm;
      double repya6a1 = -rya1a6 / ra1a6norm + rya6a1mey / ra6a1menorm;
      double repxa7a1 = -rxa1a7 / ra1a7norm + rxa7a1mex / ra7a1menorm;
      double repya7a1 = -rya1a7 / ra1a7norm + rya7a1mey / ra7a1menorm;

      double repxa2a3 = rxa2a3 / ra2a3norm + rxa2a3mex / ra2a3menorm;
      double repya2a3 = rya2a3 / ra2a3norm + rya2a3mey / ra2a3menorm;
      double repxa2a4 = rxa2a4 / ra2a4norm + rxa2a4mex / ra2a4menorm;
      double repya2a4 = rya2a4 / ra2a4norm + rya2a4mey / ra2a4menorm;
      double repxa2a5 = rxa2a5 / ra2a5norm + rxa2a5mex / ra2a5menorm;
      double repya2a5 = rya2a5 / ra2a5norm + rya2a5mey / ra2a5menorm;
      double repxa2a6 = rxa2a6 / ra2a6norm + rxa2a6mex / ra2a6menorm;
      double repya2a6 = rya2a6 / ra2a6norm + rya2a6mey / ra2a6menorm;
      double repxa2a7 = rxa2a7 / ra2a7norm + rxa2a7mex / ra2a7menorm;
      double repya2a7 = rya2a7 / ra2a7norm + rya2a7mey / ra2a7menorm;

      double repxa3a2 = -rxa2a3 / ra2a3norm + rxa3a2mex / ra3a2menorm;
      double repya3a2 = -rya2a3 / ra2a3norm + rya3a2mey / ra3a2menorm;
      double repxa4a2 = -rxa2a4 / ra2a4norm + rxa4a2mex / ra4a2menorm;
      double repya4a2 = -rya2a4 / ra2a4norm + rya4a2mey / ra4a2menorm;
      double repxa5a2 = -rxa2a5 / ra2a5norm + rxa5a2mex / ra5a2menorm;
      double repya5a2 = -rya2a5 / ra2a5norm + rya5a2mey / ra5a2menorm;
      double repxa6a2 = -rxa2a6 / ra2a6norm + rxa6a2mex / ra6a2menorm;
      double repya6a2 = -rya2a6 / ra2a6norm + rya6a2mey / ra6a2menorm;
      double repxa7a2 = -rxa2a7 / ra2a7norm + rxa7a2mex / ra7a2menorm;
      double repya7a2 = -rya2a7 / ra2a7norm + rya7a2mey / ra7a2menorm;

      double repxa3a4 = rxa3a4 / ra3a4norm + rxa3a4mex / ra3a4menorm;
      double repya3a4 = rya3a4 / ra3a4norm + rya3a4mey / ra3a4menorm;
      double repxa3a5 = rxa3a5 / ra3a5norm + rxa3a5mex / ra3a5menorm;
      double repya3a5 = rya3a5 / ra3a5norm + rya3a5mey / ra3a5menorm;
      double repxa3a6 = rxa3a6 / ra3a6norm + rxa3a6mex / ra3a6menorm;
      double repya3a6 = rya3a6 / ra3a6norm + rya3a6mey / ra3a6menorm;
      double repxa3a7 = rxa3a7 / ra3a7norm + rxa3a7mex / ra3a7menorm;
      double repya3a7 = rya3a7 / ra3a7norm + rya3a7mey / ra3a7menorm;

      double repxa4a3 = -rxa3a4 / ra3a4norm + rxa4a3mex / ra4a3menorm;
      double repya4a3 = -rya3a4 / ra3a4norm + rya4a3mey / ra4a3menorm;
      double repxa5a3 = -rxa3a5 / ra3a5norm + rxa5a3mex / ra5a3menorm;
      double repya5a3 = -rya3a5 / ra3a5norm + rya5a3mey / ra5a3menorm;
      double repxa6a3 = -rxa3a6 / ra3a6norm + rxa6a3mex / ra6a3menorm;
      double repya6a3 = -rya3a6 / ra3a6norm + rya6a3mey / ra6a3menorm;
      double repxa7a3 = -rxa3a7 / ra3a7norm + rxa7a3mex / ra7a3menorm;
      double repya7a3 = -rya3a7 / ra3a7norm + rya7a3mey / ra7a3menorm;

      double repxa4a5 = rxa4a5 / ra4a5norm + rxa4a5mex / ra4a5menorm;
      double repya4a5 = rya4a5 / ra4a5norm + rya4a5mey / ra4a5menorm;
      double repxa4a6 = rxa4a6 / ra4a6norm + rxa4a6mex / ra4a6menorm;
      double repya4a6 = rya4a6 / ra4a6norm + rya4a6mey / ra4a6menorm;
      double repxa4a7 = rxa4a7 / ra4a7norm + rxa4a7mex / ra4a7menorm;
      double repya4a7 = rya4a7 / ra4a7norm + rya4a7mey / ra4a7menorm;

      double repxa5a4 = -rxa4a5 / ra4a5norm + rxa5a4mex / ra5a4menorm;
      double repya5a4 = -rya4a5 / ra4a5norm + rya5a4mey / ra5a4menorm;
      double repxa6a4 = -rxa4a6 / ra4a6norm + rxa6a4mex / ra6a4menorm;
      double repya6a4 = -rya4a6 / ra4a6norm + rya6a4mey / ra6a4menorm;
      double repxa7a4 = -rxa4a7 / ra4a7norm + rxa7a4mex / ra7a4menorm;
      double repya7a4 = -rya4a7 / ra4a7norm + rya7a4mey / ra7a4menorm;

      double repxa5a6 = rxa5a6 / ra5a6norm + rxa5a6mex / ra5a6menorm;
      double repya5a6 = rya5a6 / ra5a6norm + rya5a6mey / ra5a6menorm;
      double repxa5a7 = rxa5a7 / ra5a7norm + rxa5a7mex / ra5a7menorm;
      double repya5a7 = rya5a7 / ra5a7norm + rya5a7mey / ra5a7menorm;

      double repxa6a5 = -rxa5a6 / ra5a6norm + rxa6a5mex / ra6a5menorm;
      double repya6a5 = -rya5a6 / ra5a6norm + rya6a5mey / ra6a5menorm;
      double repxa7a5 = -rxa5a7 / ra5a7norm + rxa7a5mex / ra7a5menorm;
      double repya7a5 = -rya5a7 / ra5a7norm + rya7a5mey / ra7a5menorm;

      double repxa6a7 = rxa6a7 / ra6a7norm + rxa6a7mex / ra6a7menorm;
      double repya6a7 = rya6a7 / ra6a7norm + rya6a7mey / ra6a7menorm;

      double repxa7a6 = -rxa6a7 / ra6a7norm + rxa7a6mex / ra7a6menorm;
      double repya7a6 = -rya6a7 / ra6a7norm + rya7a6mey / ra7a6menorm;

      double ba0a1 = normsuma0a1 * normsuma0a1 - da1 * da1;
      double ba0a2 = normsuma0a2 * normsuma0a2 - da2 * da2;
      double ba0a3 = normsuma0a3 * normsuma0a3 - da3 * da3;
      double ba0a4 = normsuma0a4 * normsuma0a4 - da4 * da4;
      double ba0a5 = normsuma0a5 * normsuma0a5 - da5 * da5;
      double ba0a6 = normsuma0a6 * normsuma0a6 - da6 * da6;
      double ba0a7 = normsuma0a7 * normsuma0a7 - da7 * da7;

      double ba1a0 = normsuma1a0 * normsuma1a0 - da0 * da0;
      double ba2a0 = normsuma2a0 * normsuma2a0 - da0 * da0;
      double ba3a0 = normsuma3a0 * normsuma3a0 - da0 * da0;
      double ba4a0 = normsuma4a0 * normsuma4a0 - da0 * da0;
      double ba5a0 = normsuma5a0 * normsuma5a0 - da0 * da0;
      double ba6a0 = normsuma6a0 * normsuma6a0 - da0 * da0;
      double ba7a0 = normsuma7a0 * normsuma7a0 - da0 * da0;

      double ba1a2 = normsuma1a2 * normsuma1a2 - da2 * da2;
      double ba1a3 = normsuma1a3 * normsuma1a3 - da3 * da3;
      double ba1a4 = normsuma1a4 * normsuma1a4 - da4 * da4;
      double ba1a5 = normsuma1a5 * normsuma1a5 - da5 * da5;
      double ba1a6 = normsuma1a6 * normsuma1a6 - da6 * da6;
      double ba1a7 = normsuma1a7 * normsuma1a7 - da7 * da7;

      double ba2a1 = normsuma2a1 * normsuma2a1 - da1 * da1;
      double ba3a1 = normsuma3a1 * normsuma3a1 - da1 * da1;
      double ba4a1 = normsuma4a1 * normsuma4a1 - da1 * da1;
      double ba5a1 = normsuma5a1 * normsuma5a1 - da1 * da1;
      double ba6a1 = normsuma6a1 * normsuma6a1 - da1 * da1;
      double ba7a1 = normsuma7a1 * normsuma7a1 - da1 * da1;

      double ba2a3 = normsuma2a3 * normsuma2a3 - da3 * da3;
      double ba2a4 = normsuma2a4 * normsuma2a4 - da4 * da4;
      double ba2a5 = normsuma2a5 * normsuma2a5 - da5 * da5;
      double ba2a6 = normsuma2a6 * normsuma2a6 - da6 * da6;
      double ba2a7 = normsuma2a7 * normsuma2a7 - da7 * da7;

      double ba3a2 = normsuma3a2 * normsuma3a2 - da2 * da2;
      double ba4a2 = normsuma4a2 * normsuma4a2 - da2 * da2;
      double ba5a2 = normsuma5a2 * normsuma5a2 - da2 * da2;
      double ba6a2 = normsuma6a2 * normsuma6a2 - da2 * da2;
      double ba7a2 = normsuma7a2 * normsuma7a2 - da2 * da2;

      double ba3a4 = normsuma3a4 * normsuma3a4 - da4 * da4;
      double ba3a5 = normsuma3a5 * normsuma3a5 - da5 * da5;
      double ba3a6 = normsuma3a6 * normsuma3a6 - da6 * da6;
      double ba3a7 = normsuma3a7 * normsuma3a7 - da7 * da7;

      double ba4a3 = normsuma4a3 * normsuma4a3 - da3 * da3;
      double ba5a3 = normsuma5a3 * normsuma5a3 - da3 * da3;
      double ba6a3 = normsuma6a3 * normsuma6a3 - da3 * da3;
      double ba7a3 = normsuma7a3 * normsuma7a3 - da3 * da3;

      double ba4a5 = normsuma4a5 * normsuma4a5 - da5 * da5;
      double ba4a6 = normsuma4a6 * normsuma4a6 - da6 * da6;
      double ba4a7 = normsuma4a7 * normsuma4a7 - da7 * da7;

      double ba5a4 = normsuma5a4 * normsuma5a4 - da4 * da4;
      double ba6a4 = normsuma6a4 * normsuma6a4 - da4 * da4;
      double ba7a4 = normsuma7a4 * normsuma7a4 - da4 * da4;

      double ba5a6 = normsuma5a6 * normsuma5a6 - da6 * da6;
      double ba5a7 = normsuma5a7 * normsuma5a7 - da7 * da7;

      double ba6a5 = normsuma6a5 * normsuma6a5 - da5 * da5;
      double ba7a5 = normsuma7a5 * normsuma7a5 - da5 * da5;

      double ba6a7 = normsuma6a7 * normsuma6a7 - da7 * da7;

      double ba7a6 = normsuma7a6 * normsuma7a6 - da6 * da6;

      double ba0a1norm = sqrt(ba0a1) * 0.5;
      double ba0a2norm = sqrt(ba0a2) * 0.5;
      double ba0a3norm = sqrt(ba0a3) * 0.5;
      double ba0a4norm = sqrt(ba0a4) * 0.5;
      double ba0a5norm = sqrt(ba0a5) * 0.5;
      double ba0a6norm = sqrt(ba0a6) * 0.5;
      double ba0a7norm = sqrt(ba0a7) * 0.5;

      double ba1a0norm = sqrt(ba1a0) * 0.5;
      double ba2a0norm = sqrt(ba2a0) * 0.5;
      double ba3a0norm = sqrt(ba3a0) * 0.5;
      double ba4a0norm = sqrt(ba4a0) * 0.5;
      double ba5a0norm = sqrt(ba5a0) * 0.5;
      double ba6a0norm = sqrt(ba6a0) * 0.5;
      double ba7a0norm = sqrt(ba7a0) * 0.5;

      double ba1a2norm = sqrt(ba1a2) * 0.5;
      double ba1a3norm = sqrt(ba1a3) * 0.5;
      double ba1a4norm = sqrt(ba1a4) * 0.5;
      double ba1a5norm = sqrt(ba1a5) * 0.5;
      double ba1a6norm = sqrt(ba1a6) * 0.5;
      double ba1a7norm = sqrt(ba1a7) * 0.5;

      double ba2a1norm = sqrt(ba2a1) * 0.5;
      double ba3a1norm = sqrt(ba3a1) * 0.5;
      double ba4a1norm = sqrt(ba4a1) * 0.5;
      double ba5a1norm = sqrt(ba5a1) * 0.5;
      double ba6a1norm = sqrt(ba6a1) * 0.5;
      double ba7a1norm = sqrt(ba7a1) * 0.5;

      double ba2a3norm = sqrt(ba2a3) * 0.5;
      double ba2a4norm = sqrt(ba2a4) * 0.5;
      double ba2a5norm = sqrt(ba2a5) * 0.5;
      double ba2a6norm = sqrt(ba2a6) * 0.5;
      double ba2a7norm = sqrt(ba2a7) * 0.5;

      double ba3a2norm = sqrt(ba3a2) * 0.5;
      double ba4a2norm = sqrt(ba4a2) * 0.5;
      double ba5a2norm = sqrt(ba5a2) * 0.5;
      double ba6a2norm = sqrt(ba6a2) * 0.5;
      double ba7a2norm = sqrt(ba7a2) * 0.5;

      double ba3a4norm = sqrt(ba3a4) * 0.5;
      double ba3a5norm = sqrt(ba3a5) * 0.5;
      double ba3a6norm = sqrt(ba3a6) * 0.5;
      double ba3a7norm = sqrt(ba3a7) * 0.5;

      double ba4a3norm = sqrt(ba4a3) * 0.5;
      double ba5a3norm = sqrt(ba5a3) * 0.5;
      double ba6a3norm = sqrt(ba6a3) * 0.5;
      double ba7a3norm = sqrt(ba7a3) * 0.5;

      double ba4a5norm = sqrt(ba4a5) * 0.5;
      double ba4a6norm = sqrt(ba4a6) * 0.5;
      double ba4a7norm = sqrt(ba4a7) * 0.5;

      double ba5a4norm = sqrt(ba5a4) * 0.5;
      double ba6a4norm = sqrt(ba6a4) * 0.5;
      double ba7a4norm = sqrt(ba7a4) * 0.5;

      double ba5a6norm = sqrt(ba5a6) * 0.5;
      double ba5a7norm = sqrt(ba5a7) * 0.5;

      double ba6a5norm = sqrt(ba6a5) * 0.5;
      double ba7a5norm = sqrt(ba7a5) * 0.5;

      double ba6a7norm = sqrt(ba6a7) * 0.5;

      double ba7a6norm = sqrt(ba7a6) * 0.5;

      double cfa0a1 = exp_fast(-ba0a1norm * inv_sigma) * normsuma0a1 * DIV_FACTOR / ba0a1norm;
      double cfa0a2 = exp_fast(-ba0a2norm * inv_sigma) * normsuma0a2 * DIV_FACTOR / ba0a2norm;
      double cfa0a3 = exp_fast(-ba0a3norm * inv_sigma) * normsuma0a3 * DIV_FACTOR / ba0a3norm;
      double cfa0a4 = exp_fast(-ba0a4norm * inv_sigma) * normsuma0a4 * DIV_FACTOR / ba0a4norm;
      double cfa0a5 = exp_fast(-ba0a5norm * inv_sigma) * normsuma0a5 * DIV_FACTOR / ba0a5norm;
      double cfa0a6 = exp_fast(-ba0a6norm * inv_sigma) * normsuma0a6 * DIV_FACTOR / ba0a6norm;
      double cfa0a7 = exp_fast(-ba0a7norm * inv_sigma) * normsuma0a7 * DIV_FACTOR / ba0a7norm;

      double cfa1a0 = exp_fast(-ba1a0norm * inv_sigma) * normsuma1a0 * DIV_FACTOR / ba1a0norm;
      double cfa2a0 = exp_fast(-ba2a0norm * inv_sigma) * normsuma2a0 * DIV_FACTOR / ba2a0norm;
      double cfa3a0 = exp_fast(-ba3a0norm * inv_sigma) * normsuma3a0 * DIV_FACTOR / ba3a0norm;
      double cfa4a0 = exp_fast(-ba4a0norm * inv_sigma) * normsuma4a0 * DIV_FACTOR / ba4a0norm;
      double cfa5a0 = exp_fast(-ba5a0norm * inv_sigma) * normsuma5a0 * DIV_FACTOR / ba5a0norm;
      double cfa6a0 = exp_fast(-ba6a0norm * inv_sigma) * normsuma6a0 * DIV_FACTOR / ba6a0norm;
      double cfa7a0 = exp_fast(-ba7a0norm * inv_sigma) * normsuma7a0 * DIV_FACTOR / ba7a0norm;

      double cfa1a2 = exp_fast(-ba1a2norm * inv_sigma) * normsuma1a2 * DIV_FACTOR / ba1a2norm;
      double cfa1a3 = exp_fast(-ba1a3norm * inv_sigma) * normsuma1a3 * DIV_FACTOR / ba1a3norm;
      double cfa1a4 = exp_fast(-ba1a4norm * inv_sigma) * normsuma1a4 * DIV_FACTOR / ba1a4norm;
      double cfa1a5 = exp_fast(-ba1a5norm * inv_sigma) * normsuma1a5 * DIV_FACTOR / ba1a5norm;
      double cfa1a6 = exp_fast(-ba1a6norm * inv_sigma) * normsuma1a6 * DIV_FACTOR / ba1a6norm;
      double cfa1a7 = exp_fast(-ba1a7norm * inv_sigma) * normsuma1a7 * DIV_FACTOR / ba1a7norm;

      double cfa2a1 = exp_fast(-ba2a1norm * inv_sigma) * normsuma2a1 * DIV_FACTOR / ba2a1norm;
      double cfa3a1 = exp_fast(-ba3a1norm * inv_sigma) * normsuma3a1 * DIV_FACTOR / ba3a1norm;
      double cfa4a1 = exp_fast(-ba4a1norm * inv_sigma) * normsuma4a1 * DIV_FACTOR / ba4a1norm;
      double cfa5a1 = exp_fast(-ba5a1norm * inv_sigma) * normsuma5a1 * DIV_FACTOR / ba5a1norm;
      double cfa6a1 = exp_fast(-ba6a1norm * inv_sigma) * normsuma6a1 * DIV_FACTOR / ba6a1norm;
      double cfa7a1 = exp_fast(-ba7a1norm * inv_sigma) * normsuma7a1 * DIV_FACTOR / ba7a1norm;

      double cfa2a3 = exp_fast(-ba2a3norm * inv_sigma) * normsuma2a3 * DIV_FACTOR / ba2a3norm;
      double cfa2a4 = exp_fast(-ba2a4norm * inv_sigma) * normsuma2a4 * DIV_FACTOR / ba2a4norm;
      double cfa2a5 = exp_fast(-ba2a5norm * inv_sigma) * normsuma2a5 * DIV_FACTOR / ba2a5norm;
      double cfa2a6 = exp_fast(-ba2a6norm * inv_sigma) * normsuma2a6 * DIV_FACTOR / ba2a6norm;
      double cfa2a7 = exp_fast(-ba2a7norm * inv_sigma) * normsuma2a7 * DIV_FACTOR / ba2a7norm;

      double cfa3a2 = exp_fast(-ba3a2norm * inv_sigma) * normsuma3a2 * DIV_FACTOR / ba3a2norm;
      double cfa4a2 = exp_fast(-ba4a2norm * inv_sigma) * normsuma4a2 * DIV_FACTOR / ba4a2norm;
      double cfa5a2 = exp_fast(-ba5a2norm * inv_sigma) * normsuma5a2 * DIV_FACTOR / ba5a2norm;
      double cfa6a2 = exp_fast(-ba6a2norm * inv_sigma) * normsuma6a2 * DIV_FACTOR / ba6a2norm;
      double cfa7a2 = exp_fast(-ba7a2norm * inv_sigma) * normsuma7a2 * DIV_FACTOR / ba7a2norm;

      double cfa3a4 = exp_fast(-ba3a4norm * inv_sigma) * normsuma3a4 * DIV_FACTOR / ba3a4norm;
      double cfa3a5 = exp_fast(-ba3a5norm * inv_sigma) * normsuma3a5 * DIV_FACTOR / ba3a5norm;
      double cfa3a6 = exp_fast(-ba3a6norm * inv_sigma) * normsuma3a6 * DIV_FACTOR / ba3a6norm;
      double cfa3a7 = exp_fast(-ba3a7norm * inv_sigma) * normsuma3a7 * DIV_FACTOR / ba3a7norm;

      double cfa4a3 = exp_fast(-ba4a3norm * inv_sigma) * normsuma4a3 * DIV_FACTOR / ba4a3norm;
      double cfa5a3 = exp_fast(-ba5a3norm * inv_sigma) * normsuma5a3 * DIV_FACTOR / ba5a3norm;
      double cfa6a3 = exp_fast(-ba6a3norm * inv_sigma) * normsuma6a3 * DIV_FACTOR / ba6a3norm;
      double cfa7a3 = exp_fast(-ba7a3norm * inv_sigma) * normsuma7a3 * DIV_FACTOR / ba7a3norm;

      double cfa4a5 = exp_fast(-ba4a5norm * inv_sigma) * normsuma4a5 * DIV_FACTOR / ba4a5norm;
      double cfa4a6 = exp_fast(-ba4a6norm * inv_sigma) * normsuma4a6 * DIV_FACTOR / ba4a6norm;
      double cfa4a7 = exp_fast(-ba4a7norm * inv_sigma) * normsuma4a7 * DIV_FACTOR / ba4a7norm;

      double cfa5a4 = exp_fast(-ba5a4norm * inv_sigma) * normsuma5a4 * DIV_FACTOR / ba5a4norm;
      double cfa6a4 = exp_fast(-ba6a4norm * inv_sigma) * normsuma6a4 * DIV_FACTOR / ba6a4norm;
      double cfa7a4 = exp_fast(-ba7a4norm * inv_sigma) * normsuma7a4 * DIV_FACTOR / ba7a4norm;

      double cfa5a6 = exp_fast(-ba5a6norm * inv_sigma) * normsuma5a6 * DIV_FACTOR / ba5a6norm;
      double cfa5a7 = exp_fast(-ba5a7norm * inv_sigma) * normsuma5a7 * DIV_FACTOR / ba5a7norm;

      double cfa6a5 = exp_fast(-ba6a5norm * inv_sigma) * normsuma6a5 * DIV_FACTOR / ba6a5norm;
      double cfa7a5 = exp_fast(-ba7a5norm * inv_sigma) * normsuma7a5 * DIV_FACTOR / ba7a5norm;

      double cfa6a7 = exp_fast(-ba6a7norm * inv_sigma) * normsuma6a7 * DIV_FACTOR / ba6a7norm;

      double cfa7a6 = exp_fast(-ba7a6norm * inv_sigma) * normsuma7a6 * DIV_FACTOR / ba7a6norm;

      double repxa0a11 = repxa0a1 * cfa0a1;
      double repya0a11 = repya0a1 * cfa0a1;
      double repxa0a21 = repxa0a2 * cfa0a2;
      double repya0a21 = repya0a2 * cfa0a2;
      double repxa0a31 = repxa0a3 * cfa0a3;
      double repya0a31 = repya0a3 * cfa0a3;
      double repxa0a41 = repxa0a4 * cfa0a4;
      double repya0a41 = repya0a4 * cfa0a4;
      double repxa0a51 = repxa0a5 * cfa0a5;
      double repya0a51 = repya0a5 * cfa0a5;
      double repxa0a61 = repxa0a6 * cfa0a6;
      double repya0a61 = repya0a6 * cfa0a6;
      double repxa0a71 = repxa0a7 * cfa0a7;
      double repya0a71 = repya0a7 * cfa0a7;

      double repxa1a01 = repxa1a0 * cfa1a0;
      double repya1a01 = repya1a0 * cfa1a0;
      double repxa2a01 = repxa2a0 * cfa2a0;
      double repya2a01 = repya2a0 * cfa2a0;
      double repxa3a01 = repxa3a0 * cfa3a0;
      double repya3a01 = repya3a0 * cfa3a0;
      double repxa4a01 = repxa4a0 * cfa4a0;
      double repya4a01 = repya4a0 * cfa4a0;
      double repxa5a01 = repxa5a0 * cfa5a0;
      double repya5a01 = repya5a0 * cfa5a0;
      double repxa6a01 = repxa6a0 * cfa6a0;
      double repya6a01 = repya6a0 * cfa6a0;
      double repxa7a01 = repxa7a0 * cfa7a0;
      double repya7a01 = repya7a0 * cfa7a0;

      double repxa1a21 = repxa1a2 * cfa1a2;
      double repya1a21 = repya1a2 * cfa1a2;
      double repxa1a31 = repxa1a3 * cfa1a3;
      double repya1a31 = repya1a3 * cfa1a3;
      double repxa1a41 = repxa1a4 * cfa1a4;
      double repya1a41 = repya1a4 * cfa1a4;
      double repxa1a51 = repxa1a5 * cfa1a5;
      double repya1a51 = repya1a5 * cfa1a5;
      double repxa1a61 = repxa1a6 * cfa1a6;
      double repya1a61 = repya1a6 * cfa1a6;
      double repxa1a71 = repxa1a7 * cfa1a7;
      double repya1a71 = repya1a7 * cfa1a7;

      double repxa2a11 = repxa2a1 * cfa2a1;
      double repya2a11 = repya2a1 * cfa2a1;
      double repxa3a11 = repxa3a1 * cfa3a1;
      double repya3a11 = repya3a1 * cfa3a1;
      double repxa4a11 = repxa4a1 * cfa4a1;
      double repya4a11 = repya4a1 * cfa4a1;
      double repxa5a11 = repxa5a1 * cfa5a1;
      double repya5a11 = repya5a1 * cfa5a1;
      double repxa6a11 = repxa6a1 * cfa6a1;
      double repya6a11 = repya6a1 * cfa6a1;
      double repxa7a11 = repxa7a1 * cfa7a1;
      double repya7a11 = repya7a1 * cfa7a1;

      double repxa2a31 = repxa2a3 * cfa2a3;
      double repya2a31 = repya2a3 * cfa2a3;
      double repxa2a41 = repxa2a4 * cfa2a4;
      double repya2a41 = repya2a4 * cfa2a4;
      double repxa2a51 = repxa2a5 * cfa2a5;
      double repya2a51 = repya2a5 * cfa2a5;
      double repxa2a61 = repxa2a6 * cfa2a6;
      double repya2a61 = repya2a6 * cfa2a6;
      double repxa2a71 = repxa2a7 * cfa2a7;
      double repya2a71 = repya2a7 * cfa2a7;

      double repxa3a21 = repxa3a2 * cfa3a2;
      double repya3a21 = repya3a2 * cfa3a2;
      double repxa4a21 = repxa4a2 * cfa4a2;
      double repya4a21 = repya4a2 * cfa4a2;
      double repxa5a21 = repxa5a2 * cfa5a2;
      double repya5a21 = repya5a2 * cfa5a2;
      double repxa6a21 = repxa6a2 * cfa6a2;
      double repya6a21 = repya6a2 * cfa6a2;
      double repxa7a21 = repxa7a2 * cfa7a2;
      double repya7a21 = repya7a2 * cfa7a2;

      double repxa3a41 = repxa3a4 * cfa3a4;
      double repya3a41 = repya3a4 * cfa3a4;
      double repxa3a51 = repxa3a5 * cfa3a5;
      double repya3a51 = repya3a5 * cfa3a5;
      double repxa3a61 = repxa3a6 * cfa3a6;
      double repya3a61 = repya3a6 * cfa3a6;
      double repxa3a71 = repxa3a7 * cfa3a7;
      double repya3a71 = repya3a7 * cfa3a7;

      double repxa4a31 = repxa4a3 * cfa4a3;
      double repya4a31 = repya4a3 * cfa4a3;
      double repxa5a31 = repxa5a3 * cfa5a3;
      double repya5a31 = repya5a3 * cfa5a3;
      double repxa6a31 = repxa6a3 * cfa6a3;
      double repya6a31 = repya6a3 * cfa6a3;
      double repxa7a31 = repxa7a3 * cfa7a3;
      double repya7a31 = repya7a3 * cfa7a3;

      double repxa4a51 = repxa4a5 * cfa4a5;
      double repya4a51 = repya4a5 * cfa4a5;
      double repxa4a61 = repxa4a6 * cfa4a6;
      double repya4a61 = repya4a6 * cfa4a6;
      double repxa4a71 = repxa4a7 * cfa4a7;
      double repya4a71 = repya4a7 * cfa4a7;

      double repxa5a41 = repxa5a4 * cfa5a4;
      double repya5a41 = repya5a4 * cfa5a4;
      double repxa6a41 = repxa6a4 * cfa6a4;
      double repya6a41 = repya6a4 * cfa6a4;
      double repxa7a41 = repxa7a4 * cfa7a4;
      double repya7a41 = repya7a4 * cfa7a4;

      double repxa5a61 = repxa5a6 * cfa5a6;
      double repya5a61 = repya5a6 * cfa5a6;
      double repxa5a71 = repxa5a7 * cfa5a7;
      double repya5a71 = repya5a7 * cfa5a7;

      double repxa6a51 = repxa6a5 * cfa6a5;
      double repya6a51 = repya6a5 * cfa6a5;
      double repxa7a51 = repxa7a5 * cfa7a5;
      double repya7a51 = repya7a5 * cfa7a5;

      double repxa6a71 = repxa6a7 * cfa6a7;
      double repya6a71 = repya6a7 * cfa6a7;

      double repxa7a61 = repxa7a6 * cfa7a6;
      double repya7a61 = repya7a6 * cfa7a6;

      double ca0a1 = exa0 * repxa0a11 + eya0 * repya0a11;
      double ca0a2 = exa0 * repxa0a21 + eya0 * repya0a21;
      double ca0a3 = exa0 * repxa0a31 + eya0 * repya0a31;
      double ca0a4 = exa0 * repxa0a41 + eya0 * repya0a41;
      double ca0a5 = exa0 * repxa0a51 + eya0 * repya0a51;
      double ca0a6 = exa0 * repxa0a61 + eya0 * repya0a61;
      double ca0a7 = exa0 * repxa0a71 + eya0 * repya0a71;

      double ca1a0 = exa1 * repxa1a01 + eya1 * repya1a01;
      double ca2a0 = exa2 * repxa2a01 + eya2 * repya2a01;
      double ca3a0 = exa3 * repxa3a01 + eya3 * repya3a01;
      double ca4a0 = exa4 * repxa4a01 + eya4 * repya4a01;
      double ca5a0 = exa5 * repxa5a01 + eya5 * repya5a01;
      double ca6a0 = exa6 * repxa6a01 + eya6 * repya6a01;
      double ca7a0 = exa7 * repxa7a01 + eya7 * repya7a01;

      double ca1a2 = exa1 * repxa1a21 + eya1 * repya1a21;
      double ca1a3 = exa1 * repxa1a31 + eya1 * repya1a31;
      double ca1a4 = exa1 * repxa1a41 + eya1 * repya1a41;
      double ca1a5 = exa1 * repxa1a51 + eya1 * repya1a51;
      double ca1a6 = exa1 * repxa1a61 + eya1 * repya1a61;
      double ca1a7 = exa1 * repxa1a71 + eya1 * repya1a71;

      double ca2a1 = exa2 * repxa2a11 + eya2 * repya2a11;
      double ca3a1 = exa3 * repxa3a11 + eya3 * repya3a11;
      double ca4a1 = exa4 * repxa4a11 + eya4 * repya4a11;
      double ca5a1 = exa5 * repxa5a11 + eya5 * repya5a11;
      double ca6a1 = exa6 * repxa6a11 + eya6 * repya6a11;
      double ca7a1 = exa7 * repxa7a11 + eya7 * repya7a11;

      double ca2a3 = exa2 * repxa2a31 + eya2 * repya2a31;
      double ca2a4 = exa2 * repxa2a41 + eya2 * repya2a41;
      double ca2a5 = exa2 * repxa2a51 + eya2 * repya2a51;
      double ca2a6 = exa2 * repxa2a61 + eya2 * repya2a61;
      double ca2a7 = exa2 * repxa2a71 + eya2 * repya2a71;

      double ca3a2 = exa3 * repxa3a21 + eya3 * repya3a21;
      double ca4a2 = exa4 * repxa4a21 + eya4 * repya4a21;
      double ca5a2 = exa5 * repxa5a21 + eya5 * repya5a21;
      double ca6a2 = exa6 * repxa6a21 + eya6 * repya6a21;
      double ca7a2 = exa7 * repxa7a21 + eya7 * repya7a21;

      double ca3a4 = exa3 * repxa3a41 + eya3 * repya3a41;
      double ca3a5 = exa3 * repxa3a51 + eya3 * repya3a51;
      double ca3a6 = exa3 * repxa3a61 + eya3 * repya3a61;
      double ca3a7 = exa3 * repxa3a71 + eya3 * repya3a71;

      double ca4a3 = exa4 * repxa4a31 + eya4 * repya4a31;
      double ca5a3 = exa5 * repxa5a31 + eya5 * repya5a31;
      double ca6a3 = exa6 * repxa6a31 + eya6 * repya6a31;
      double ca7a3 = exa7 * repxa7a31 + eya7 * repya7a31;

      double ca4a5 = exa4 * repxa4a51 + eya4 * repya4a51;
      double ca4a6 = exa4 * repxa4a61 + eya4 * repya4a61;
      double ca4a7 = exa4 * repxa4a71 + eya4 * repya4a71;

      double ca5a4 = exa5 * repxa5a41 + eya5 * repya5a41;
      double ca6a4 = exa6 * repxa6a41 + eya6 * repya6a41;
      double ca7a4 = exa7 * repxa7a41 + eya7 * repya7a41;

      double ca5a6 = exa5 * repxa5a61 + eya5 * repya5a61;
      double ca5a7 = exa5 * repxa5a71 + eya5 * repya5a71;

      double ca6a5 = exa6 * repxa6a51 + eya6 * repya6a51;
      double ca7a5 = exa7 * repxa7a51 + eya7 * repya7a51;

      double ca6a7 = exa6 * repxa6a71 + eya6 * repya6a71;

      double ca7a6 = exa7 * repxa7a61 + eya7 * repya7a61;

      double tha0a1 = sqrt(repxa0a11 * repxa0a11 + repya0a11 * repya0a11) * cospsi;
      double tha0a2 = sqrt(repxa0a21 * repxa0a21 + repya0a21 * repya0a21) * cospsi;
      double tha0a3 = sqrt(repxa0a31 * repxa0a31 + repya0a31 * repya0a31) * cospsi;
      double tha0a4 = sqrt(repxa0a41 * repxa0a41 + repya0a41 * repya0a41) * cospsi;
      double tha0a5 = sqrt(repxa0a51 * repxa0a51 + repya0a51 * repya0a51) * cospsi;
      double tha0a6 = sqrt(repxa0a61 * repxa0a61 + repya0a61 * repya0a61) * cospsi;
      double tha0a7 = sqrt(repxa0a71 * repxa0a71 + repya0a71 * repya0a71) * cospsi;

      double tha1a0 = sqrt(repxa1a01 * repxa1a01 + repya1a01 * repya1a01) * cospsi;
      double tha2a0 = sqrt(repxa2a01 * repxa2a01 + repya2a01 * repya2a01) * cospsi;
      double tha3a0 = sqrt(repxa3a01 * repxa3a01 + repya3a01 * repya3a01) * cospsi;
      double tha4a0 = sqrt(repxa4a01 * repxa4a01 + repya4a01 * repya4a01) * cospsi;
      double tha5a0 = sqrt(repxa5a01 * repxa5a01 + repya5a01 * repya5a01) * cospsi;
      double tha6a0 = sqrt(repxa6a01 * repxa6a01 + repya6a01 * repya6a01) * cospsi;
      double tha7a0 = sqrt(repxa7a01 * repxa7a01 + repya7a01 * repya7a01) * cospsi;

      double tha1a2 = sqrt(repxa1a21 * repxa1a21 + repya1a21 * repya1a21) * cospsi;
      double tha1a3 = sqrt(repxa1a31 * repxa1a31 + repya1a31 * repya1a31) * cospsi;
      double tha1a4 = sqrt(repxa1a41 * repxa1a41 + repya1a41 * repya1a41) * cospsi;
      double tha1a5 = sqrt(repxa1a51 * repxa1a51 + repya1a51 * repya1a51) * cospsi;
      double tha1a6 = sqrt(repxa1a61 * repxa1a61 + repya1a61 * repya1a61) * cospsi;
      double tha1a7 = sqrt(repxa1a71 * repxa1a71 + repya1a71 * repya1a71) * cospsi;

      double tha2a1 = sqrt(repxa2a11 * repxa2a11 + repya2a11 * repya2a11) * cospsi;
      double tha3a1 = sqrt(repxa3a11 * repxa3a11 + repya3a11 * repya3a11) * cospsi;
      double tha4a1 = sqrt(repxa4a11 * repxa4a11 + repya4a11 * repya4a11) * cospsi;
      double tha5a1 = sqrt(repxa5a11 * repxa5a11 + repya5a11 * repya5a11) * cospsi;
      double tha6a1 = sqrt(repxa6a11 * repxa6a11 + repya6a11 * repya6a11) * cospsi;
      double tha7a1 = sqrt(repxa7a11 * repxa7a11 + repya7a11 * repya7a11) * cospsi;

      double tha2a3 = sqrt(repxa2a31 * repxa2a31 + repya2a31 * repya2a31) * cospsi;
      double tha2a4 = sqrt(repxa2a41 * repxa2a41 + repya2a41 * repya2a41) * cospsi;
      double tha2a5 = sqrt(repxa2a51 * repxa2a51 + repya2a51 * repya2a51) * cospsi;
      double tha2a6 = sqrt(repxa2a61 * repxa2a61 + repya2a61 * repya2a61) * cospsi;
      double tha2a7 = sqrt(repxa2a71 * repxa2a71 + repya2a71 * repya2a71) * cospsi;

      double tha3a2 = sqrt(repxa3a21 * repxa3a21 + repya3a21 * repya3a21) * cospsi;
      double tha4a2 = sqrt(repxa4a21 * repxa4a21 + repya4a21 * repya4a21) * cospsi;
      double tha5a2 = sqrt(repxa5a21 * repxa5a21 + repya5a21 * repya5a21) * cospsi;
      double tha6a2 = sqrt(repxa6a21 * repxa6a21 + repya6a21 * repya6a21) * cospsi;
      double tha7a2 = sqrt(repxa7a21 * repxa7a21 + repya7a21 * repya7a21) * cospsi;

      double tha3a4 = sqrt(repxa3a41 * repxa3a41 + repya3a41 * repya3a41) * cospsi;
      double tha3a5 = sqrt(repxa3a51 * repxa3a51 + repya3a51 * repya3a51) * cospsi;
      double tha3a6 = sqrt(repxa3a61 * repxa3a61 + repya3a61 * repya3a61) * cospsi;
      double tha3a7 = sqrt(repxa3a71 * repxa3a71 + repya3a71 * repya3a71) * cospsi;

      double tha4a3 = sqrt(repxa4a31 * repxa4a31 + repya4a31 * repya4a31) * cospsi;
      double tha5a3 = sqrt(repxa5a31 * repxa5a31 + repya5a31 * repya5a31) * cospsi;
      double tha6a3 = sqrt(repxa6a31 * repxa6a31 + repya6a31 * repya6a31) * cospsi;
      double tha7a3 = sqrt(repxa7a31 * repxa7a31 + repya7a31 * repya7a31) * cospsi;

      double tha4a5 = sqrt(repxa4a51 * repxa4a51 + repya4a51 * repya4a51) * cospsi;
      double tha4a6 = sqrt(repxa4a61 * repxa4a61 + repya4a61 * repya4a61) * cospsi;
      double tha4a7 = sqrt(repxa4a71 * repxa4a71 + repya4a71 * repya4a71) * cospsi;

      double tha5a4 = sqrt(repxa5a41 * repxa5a41 + repya5a41 * repya5a41) * cospsi;
      double tha6a4 = sqrt(repxa6a41 * repxa6a41 + repya6a41 * repya6a41) * cospsi;
      double tha7a4 = sqrt(repxa7a41 * repxa7a41 + repya7a41 * repya7a41) * cospsi;

      double tha5a6 = sqrt(repxa5a61 * repxa5a61 + repya5a61 * repya5a61) * cospsi;
      double tha5a7 = sqrt(repxa5a71 * repxa5a71 + repya5a71 * repya5a71) * cospsi;

      double tha6a5 = sqrt(repxa6a51 * repxa6a51 + repya6a51 * repya6a51) * cospsi;
      double tha7a5 = sqrt(repxa7a51 * repxa7a51 + repya7a51 * repya7a51) * cospsi;

      double tha6a7 = sqrt(repxa6a71 * repxa6a71 + repya6a71 * repya6a71) * cospsi;

      double tha7a6 = sqrt(repxa7a61 * repxa7a61 + repya7a61 * repya7a61) * cospsi;

      double wa0a1 = -ca0a1 >= tha0a1 ? 1 : INFLUENCE;
      double wa0a2 = -ca0a2 >= tha0a2 ? 1 : INFLUENCE;
      double wa0a3 = -ca0a3 >= tha0a3 ? 1 : INFLUENCE;
      double wa0a4 = -ca0a4 >= tha0a4 ? 1 : INFLUENCE;
      double wa0a5 = -ca0a5 >= tha0a5 ? 1 : INFLUENCE;
      double wa0a6 = -ca0a6 >= tha0a6 ? 1 : INFLUENCE;
      double wa0a7 = -ca0a7 >= tha0a7 ? 1 : INFLUENCE;

      double wa1a0 = -ca1a0 >= tha1a0 ? 1 : INFLUENCE;
      double wa2a0 = -ca2a0 >= tha2a0 ? 1 : INFLUENCE;
      double wa3a0 = -ca3a0 >= tha3a0 ? 1 : INFLUENCE;
      double wa4a0 = -ca4a0 >= tha4a0 ? 1 : INFLUENCE;
      double wa5a0 = -ca5a0 >= tha5a0 ? 1 : INFLUENCE;
      double wa6a0 = -ca6a0 >= tha6a0 ? 1 : INFLUENCE;
      double wa7a0 = -ca7a0 >= tha7a0 ? 1 : INFLUENCE;

      double wa1a2 = -ca1a2 >= tha1a2 ? 1 : INFLUENCE;
      double wa1a3 = -ca1a3 >= tha1a3 ? 1 : INFLUENCE;
      double wa1a4 = -ca1a4 >= tha1a4 ? 1 : INFLUENCE;
      double wa1a5 = -ca1a5 >= tha1a5 ? 1 : INFLUENCE;
      double wa1a6 = -ca1a6 >= tha1a6 ? 1 : INFLUENCE;
      double wa1a7 = -ca1a7 >= tha1a7 ? 1 : INFLUENCE;

      double wa2a1 = -ca2a1 >= tha2a1 ? 1 : INFLUENCE;
      double wa3a1 = -ca3a1 >= tha3a1 ? 1 : INFLUENCE;
      double wa4a1 = -ca4a1 >= tha4a1 ? 1 : INFLUENCE;
      double wa5a1 = -ca5a1 >= tha5a1 ? 1 : INFLUENCE;
      double wa6a1 = -ca6a1 >= tha6a1 ? 1 : INFLUENCE;
      double wa7a1 = -ca7a1 >= tha7a1 ? 1 : INFLUENCE;

      double wa2a3 = -ca2a3 >= tha2a3 ? 1 : INFLUENCE;
      double wa2a4 = -ca2a4 >= tha2a4 ? 1 : INFLUENCE;
      double wa2a5 = -ca2a5 >= tha2a5 ? 1 : INFLUENCE;
      double wa2a6 = -ca2a6 >= tha2a6 ? 1 : INFLUENCE;
      double wa2a7 = -ca2a7 >= tha2a7 ? 1 : INFLUENCE;

      double wa3a2 = -ca3a2 >= tha3a2 ? 1 : INFLUENCE;
      double wa4a2 = -ca4a2 >= tha4a2 ? 1 : INFLUENCE;
      double wa5a2 = -ca5a2 >= tha5a2 ? 1 : INFLUENCE;
      double wa6a2 = -ca6a2 >= tha6a2 ? 1 : INFLUENCE;
      double wa7a2 = -ca7a2 >= tha7a2 ? 1 : INFLUENCE;

      double wa3a4 = -ca3a4 >= tha3a4 ? 1 : INFLUENCE;
      double wa3a5 = -ca3a5 >= tha3a5 ? 1 : INFLUENCE;
      double wa3a6 = -ca3a6 >= tha3a6 ? 1 : INFLUENCE;
      double wa3a7 = -ca3a7 >= tha3a7 ? 1 : INFLUENCE;

      double wa4a3 = -ca4a3 >= tha4a3 ? 1 : INFLUENCE;
      double wa5a3 = -ca5a3 >= tha5a3 ? 1 : INFLUENCE;
      double wa6a3 = -ca6a3 >= tha6a3 ? 1 : INFLUENCE;
      double wa7a3 = -ca7a3 >= tha7a3 ? 1 : INFLUENCE;

      double wa4a5 = -ca4a5 >= tha4a5 ? 1 : INFLUENCE;
      double wa4a6 = -ca4a6 >= tha4a6 ? 1 : INFLUENCE;
      double wa4a7 = -ca4a7 >= tha4a7 ? 1 : INFLUENCE;

      double wa5a4 = -ca5a4 >= tha5a4 ? 1 : INFLUENCE;
      double wa6a4 = -ca6a4 >= tha6a4 ? 1 : INFLUENCE;
      double wa7a4 = -ca7a4 >= tha7a4 ? 1 : INFLUENCE;

      double wa5a6 = -ca5a6 >= tha5a6 ? 1 : INFLUENCE;
      double wa5a7 = -ca5a7 >= tha5a7 ? 1 : INFLUENCE;

      double wa6a5 = -ca6a5 >= tha6a5 ? 1 : INFLUENCE;
      double wa7a5 = -ca7a5 >= tha7a5 ? 1 : INFLUENCE;

      double wa6a7 = -ca6a7 >= tha6a7 ? 1 : INFLUENCE;

      double wa7a6 = -ca7a6 >= tha7a6 ? 1 : INFLUENCE;

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

      //iterate over all people
      for (int j = i + 8; j < n; j++) // for (int j = i + 1; j < n; j++)
      {
        double rxb = position[IndexX(j)];
        double ryb = position[IndexY(j, n)];
        double exb0 = desired_direction[IndexX(j)];
        double eyb0 = desired_direction[IndexY(j, n)];
        double db0 = speed[j];

        double rxab0 = rxa0 - rxb;
        double ryab0 = rya0 - ryb;
        double rxab1 = rxa1 - rxb;
        double ryab1 = rya1 - ryb;
        double rxab2 = rxa2 - rxb;
        double ryab2 = rya2 - ryb;
        double rxab3 = rxa3 - rxb;
        double ryab3 = rya3 - ryb;
        double rxab4 = rxa4 - rxb;
        double ryab4 = rya4 - ryb;
        double rxab5 = rxa5 - rxb;
        double ryab5 = rya5 - ryb;
        double rxab6 = rxa6 - rxb;
        double ryab6 = rya6 - ryb;
        double rxab7 = rxa7 - rxb;
        double ryab7 = rya7 - ryb;

        double rab0 = rxab0 * rxab0 + ryab0 * ryab0;
        double rab1 = rxab1 * rxab1 + ryab1 * ryab1;
        double rab2 = rxab2 * rxab2 + ryab2 * ryab2;
        double rab3 = rxab3 * rxab3 + ryab3 * ryab3;
        double rab4 = rxab4 * rxab4 + ryab4 * ryab4;
        double rab5 = rxab5 * rxab5 + ryab5 * ryab5;
        double rab6 = rxab6 * rxab6 + ryab6 * ryab6;
        double rab7 = rxab7 * rxab7 + ryab7 * ryab7;

        double rabnorm0 = sqrt(rab0);
        double rabnorm1 = sqrt(rab1);
        double rabnorm2 = sqrt(rab2);
        double rabnorm3 = sqrt(rab3);
        double rabnorm4 = sqrt(rab4);
        double rabnorm5 = sqrt(rab5);
        double rabnorm6 = sqrt(rab6);
        double rabnorm7 = sqrt(rab7);

        //double everything

        double rxabmex0 = rxab0 - db0 * exb0;
        double ryabmey0 = ryab0 - db0 * eyb0;
        double rxabmex1 = rxab1 - db0 * exb0;
        double ryabmey1 = ryab1 - db0 * eyb0;
        double rxabmex2 = rxab2 - db0 * exb0;
        double ryabmey2 = ryab2 - db0 * eyb0;
        double rxabmex3 = rxab3 - db0 * exb0;
        double ryabmey3 = ryab3 - db0 * eyb0;
        double rxabmex4 = rxab4 - db0 * exb0;
        double ryabmey4 = ryab4 - db0 * eyb0;
        double rxabmex5 = rxab5 - db0 * exb0;
        double ryabmey5 = ryab5 - db0 * eyb0;
        double rxabmex6 = rxab6 - db0 * exb0;
        double ryabmey6 = ryab6 - db0 * eyb0;
        double rxabmex7 = rxab7 - db0 * exb0;
        double ryabmey7 = ryab7 - db0 * eyb0;

        double rxbamex0 = -rxab0 - da0 * exa0;
        double rybamey0 = -ryab0 - da0 * eya0;
        double rxbamex1 = -rxab1 - da1 * exa1;
        double rybamey1 = -ryab1 - da1 * eya1;
        double rxbamex2 = -rxab2 - da2 * exa2;
        double rybamey2 = -ryab2 - da2 * eya2;
        double rxbamex3 = -rxab3 - da3 * exa3;
        double rybamey3 = -ryab3 - da3 * eya3;
        double rxbamex4 = -rxab4 - da4 * exa4;
        double rybamey4 = -ryab4 - da4 * eya4;
        double rxbamex5 = -rxab5 - da5 * exa5;
        double rybamey5 = -ryab5 - da5 * eya5;
        double rxbamex6 = -rxab6 - da6 * exa6;
        double rybamey6 = -ryab6 - da6 * eya6;
        double rxbamex7 = -rxab7 - da7 * exa7;
        double rybamey7 = -ryab7 - da7 * eya7;

        double rabme0 = rxabmex0 * rxabmex0 + ryabmey0 * ryabmey0;
        double rabme1 = rxabmex1 * rxabmex1 + ryabmey1 * ryabmey1;
        double rabme2 = rxabmex2 * rxabmex2 + ryabmey2 * ryabmey2;
        double rabme3 = rxabmex3 * rxabmex3 + ryabmey3 * ryabmey3;
        double rabme4 = rxabmex4 * rxabmex4 + ryabmey4 * ryabmey4;
        double rabme5 = rxabmex5 * rxabmex5 + ryabmey5 * ryabmey5;
        double rabme6 = rxabmex6 * rxabmex6 + ryabmey6 * ryabmey6;
        double rabme7 = rxabmex7 * rxabmex7 + ryabmey7 * ryabmey7;

        double rbame0 = rxbamex0 * rxbamex0 + rybamey0 * rybamey0;
        double rbame1 = rxbamex1 * rxbamex1 + rybamey1 * rybamey1;
        double rbame2 = rxbamex2 * rxbamex2 + rybamey2 * rybamey2;
        double rbame3 = rxbamex3 * rxbamex3 + rybamey3 * rybamey3;
        double rbame4 = rxbamex4 * rxbamex4 + rybamey4 * rybamey4;
        double rbame5 = rxbamex5 * rxbamex5 + rybamey5 * rybamey5;
        double rbame6 = rxbamex6 * rxbamex6 + rybamey6 * rybamey6;
        double rbame7 = rxbamex7 * rxbamex7 + rybamey7 * rybamey7;

        double rabmenorm0 = sqrt(rabme0);
        double rabmenorm1 = sqrt(rabme1);
        double rabmenorm2 = sqrt(rabme2);
        double rabmenorm3 = sqrt(rabme3);
        double rabmenorm4 = sqrt(rabme4);
        double rabmenorm5 = sqrt(rabme5);
        double rabmenorm6 = sqrt(rabme6);
        double rabmenorm7 = sqrt(rabme7);

        double rbamenorm0 = sqrt(rbame0);
        double rbamenorm1 = sqrt(rbame1);
        double rbamenorm2 = sqrt(rbame2);
        double rbamenorm3 = sqrt(rbame3);
        double rbamenorm4 = sqrt(rbame4);
        double rbamenorm5 = sqrt(rbame5);
        double rbamenorm6 = sqrt(rbame6);
        double rbamenorm7 = sqrt(rbame7);

        double normsumab0 = rabnorm0 + rabmenorm0;
        double normsumab1 = rabnorm1 + rabmenorm1;
        double normsumab2 = rabnorm2 + rabmenorm2;
        double normsumab3 = rabnorm3 + rabmenorm3;
        double normsumab4 = rabnorm4 + rabmenorm4;
        double normsumab5 = rabnorm5 + rabmenorm5;
        double normsumab6 = rabnorm6 + rabmenorm6;
        double normsumab7 = rabnorm7 + rabmenorm7;

        double normsumba0 = rabnorm0 + rbamenorm0;
        double normsumba1 = rabnorm1 + rbamenorm1;
        double normsumba2 = rabnorm2 + rbamenorm2;
        double normsumba3 = rabnorm3 + rbamenorm3;
        double normsumba4 = rabnorm4 + rbamenorm4;
        double normsumba5 = rabnorm5 + rbamenorm5;
        double normsumba6 = rabnorm6 + rbamenorm6;
        double normsumba7 = rabnorm7 + rbamenorm7;

        double repxab00 = rxab0 / rabnorm0 + rxabmex0 / rabmenorm0;
        double repyab00 = ryab0 / rabnorm0 + ryabmey0 / rabmenorm0;
        double repxab01 = rxab1 / rabnorm1 + rxabmex1 / rabmenorm1;
        double repyab01 = ryab1 / rabnorm1 + ryabmey1 / rabmenorm1;
        double repxab02 = rxab2 / rabnorm2 + rxabmex2 / rabmenorm2;
        double repyab02 = ryab2 / rabnorm2 + ryabmey2 / rabmenorm2;
        double repxab03 = rxab3 / rabnorm3 + rxabmex3 / rabmenorm3;
        double repyab03 = ryab3 / rabnorm3 + ryabmey3 / rabmenorm3;
        double repxab04 = rxab4 / rabnorm4 + rxabmex4 / rabmenorm4;
        double repyab04 = ryab4 / rabnorm4 + ryabmey4 / rabmenorm4;
        double repxab05 = rxab5 / rabnorm5 + rxabmex5 / rabmenorm5;
        double repyab05 = ryab5 / rabnorm5 + ryabmey5 / rabmenorm5;
        double repxab06 = rxab6 / rabnorm6 + rxabmex6 / rabmenorm6;
        double repyab06 = ryab6 / rabnorm6 + ryabmey6 / rabmenorm6;
        double repxab07 = rxab7 / rabnorm7 + rxabmex7 / rabmenorm7;
        double repyab07 = ryab7 / rabnorm7 + ryabmey7 / rabmenorm7;

        double repxba00 = -rxab0 / rabnorm0 + rxbamex0 / rbamenorm0;
        double repyba00 = -ryab0 / rabnorm0 + rybamey0 / rbamenorm0;
        double repxba01 = -rxab1 / rabnorm1 + rxbamex1 / rbamenorm1;
        double repyba01 = -ryab1 / rabnorm1 + rybamey1 / rbamenorm1;
        double repxba02 = -rxab2 / rabnorm2 + rxbamex2 / rbamenorm2;
        double repyba02 = -ryab2 / rabnorm2 + rybamey2 / rbamenorm2;
        double repxba03 = -rxab3 / rabnorm3 + rxbamex3 / rbamenorm3;
        double repyba03 = -ryab3 / rabnorm3 + rybamey3 / rbamenorm3;
        double repxba04 = -rxab4 / rabnorm4 + rxbamex4 / rbamenorm4;
        double repyba04 = -ryab4 / rabnorm4 + rybamey4 / rbamenorm4;
        double repxba05 = -rxab5 / rabnorm5 + rxbamex5 / rbamenorm5;
        double repyba05 = -ryab5 / rabnorm5 + rybamey5 / rbamenorm5;
        double repxba06 = -rxab6 / rabnorm6 + rxbamex6 / rbamenorm6;
        double repyba06 = -ryab6 / rabnorm6 + rybamey6 / rbamenorm6;
        double repxba07 = -rxab7 / rabnorm7 + rxbamex7 / rbamenorm7;
        double repyba07 = -ryab7 / rabnorm7 + rybamey7 / rbamenorm7;

        double bab0 = normsumab0 * normsumab0 - db0 * db0;
        double bab1 = normsumab1 * normsumab1 - db0 * db0;
        double bab2 = normsumab2 * normsumab2 - db0 * db0;
        double bab3 = normsumab3 * normsumab3 - db0 * db0;
        double bab4 = normsumab4 * normsumab4 - db0 * db0;
        double bab5 = normsumab5 * normsumab5 - db0 * db0;
        double bab6 = normsumab6 * normsumab6 - db0 * db0;
        double bab7 = normsumab7 * normsumab7 - db0 * db0;

        double bba0 = normsumba0 * normsumba0 - da0 * da0;
        double bba1 = normsumba1 * normsumba1 - da1 * da1;
        double bba2 = normsumba2 * normsumba2 - da2 * da2;
        double bba3 = normsumba3 * normsumba3 - da3 * da3;
        double bba4 = normsumba4 * normsumba4 - da4 * da4;
        double bba5 = normsumba5 * normsumba5 - da5 * da5;
        double bba6 = normsumba6 * normsumba6 - da6 * da6;
        double bba7 = normsumba7 * normsumba7 - da7 * da7;

        double babnorm0 = sqrt(bab0) * 0.5;
        double babnorm1 = sqrt(bab1) * 0.5;
        double babnorm2 = sqrt(bab2) * 0.5;
        double babnorm3 = sqrt(bab3) * 0.5;
        double babnorm4 = sqrt(bab4) * 0.5;
        double babnorm5 = sqrt(bab5) * 0.5;
        double babnorm6 = sqrt(bab6) * 0.5;
        double babnorm7 = sqrt(bab7) * 0.5;

        double bbanorm0 = sqrt(bba0) * 0.5;
        double bbanorm1 = sqrt(bba1) * 0.5;
        double bbanorm2 = sqrt(bba2) * 0.5;
        double bbanorm3 = sqrt(bba3) * 0.5;
        double bbanorm4 = sqrt(bba4) * 0.5;
        double bbanorm5 = sqrt(bba5) * 0.5;
        double bbanorm6 = sqrt(bba6) * 0.5;
        double bbanorm7 = sqrt(bba7) * 0.5;

        double cfab0 = exp_fast(-babnorm0 * inv_sigma) * normsumab0 * DIV_FACTOR / babnorm0;
        double cfab1 = exp_fast(-babnorm1 * inv_sigma) * normsumab1 * DIV_FACTOR / babnorm1;
        double cfab2 = exp_fast(-babnorm2 * inv_sigma) * normsumab2 * DIV_FACTOR / babnorm2;
        double cfab3 = exp_fast(-babnorm3 * inv_sigma) * normsumab3 * DIV_FACTOR / babnorm3;
        double cfab4 = exp_fast(-babnorm4 * inv_sigma) * normsumab4 * DIV_FACTOR / babnorm4;
        double cfab5 = exp_fast(-babnorm5 * inv_sigma) * normsumab5 * DIV_FACTOR / babnorm5;
        double cfab6 = exp_fast(-babnorm6 * inv_sigma) * normsumab6 * DIV_FACTOR / babnorm6;
        double cfab7 = exp_fast(-babnorm7 * inv_sigma) * normsumab7 * DIV_FACTOR / babnorm7;

        double cfba0 = exp_fast(-bbanorm0 * inv_sigma) * normsumba0 * DIV_FACTOR / bbanorm0; //2 mult, 2 div, 1 exp
        double cfba1 = exp_fast(-bbanorm1 * inv_sigma) * normsumba1 * DIV_FACTOR / bbanorm1; //2 mult, 2 div, 1 exp
        double cfba2 = exp_fast(-bbanorm2 * inv_sigma) * normsumba2 * DIV_FACTOR / bbanorm2; //2 mult, 2 div, 1 exp
        double cfba3 = exp_fast(-bbanorm3 * inv_sigma) * normsumba3 * DIV_FACTOR / bbanorm3; //2 mult, 2 div, 1 exp
        double cfba4 = exp_fast(-bbanorm4 * inv_sigma) * normsumba4 * DIV_FACTOR / bbanorm4; //2 mult, 2 div, 1 exp
        double cfba5 = exp_fast(-bbanorm5 * inv_sigma) * normsumba5 * DIV_FACTOR / bbanorm5; //2 mult, 2 div, 1 exp
        double cfba6 = exp_fast(-bbanorm6 * inv_sigma) * normsumba6 * DIV_FACTOR / bbanorm6; //2 mult, 2 div, 1 exp
        double cfba7 = exp_fast(-bbanorm7 * inv_sigma) * normsumba7 * DIV_FACTOR / bbanorm7; //2 mult, 2 div, 1 exp

        double repxab10 = repxab00 * cfab0; //1 mult
        double repyab10 = repyab00 * cfab0; //1 mult
        double repxab11 = repxab01 * cfab1; //1 mult
        double repyab11 = repyab01 * cfab1; //1 mult
        double repxab12 = repxab02 * cfab2; //1 mult
        double repyab12 = repyab02 * cfab2; //1 mult
        double repxab13 = repxab03 * cfab3; //1 mult
        double repyab13 = repyab03 * cfab3; //1 mult
        double repxab14 = repxab04 * cfab4; //1 mult
        double repyab14 = repyab04 * cfab4; //1 mult
        double repxab15 = repxab05 * cfab5; //1 mult
        double repyab15 = repyab05 * cfab5; //1 mult
        double repxab16 = repxab06 * cfab6; //1 mult
        double repyab16 = repyab06 * cfab6; //1 mult
        double repxab17 = repxab07 * cfab7; //1 mult
        double repyab17 = repyab07 * cfab7; //1 mult

        double repxba10 = repxba00 * cfba0; //1 mult
        double repyba10 = repyba00 * cfba0; //1 mult
        double repxba11 = repxba01 * cfba1; //1 mult
        double repyba11 = repyba01 * cfba1; //1 mult
        double repxba12 = repxba02 * cfba2; //1 mult
        double repyba12 = repyba02 * cfba2; //1 mult
        double repxba13 = repxba03 * cfba3; //1 mult
        double repyba13 = repyba03 * cfba3; //1 mult
        double repxba14 = repxba04 * cfba4; //1 mult
        double repyba14 = repyba04 * cfba4; //1 mult
        double repxba15 = repxba05 * cfba5; //1 mult
        double repyba15 = repyba05 * cfba5; //1 mult
        double repxba16 = repxba06 * cfba6; //1 mult
        double repyba16 = repyba06 * cfba6; //1 mult
        double repxba17 = repxba07 * cfba7; //1 mult
        double repyba17 = repyba07 * cfba7; //1 mult

        double cab0 = exa0 * repxab10 + eya0 * repyab10; //1 add, 2 mult
        double cab1 = exa1 * repxab11 + eya1 * repyab11; //1 add, 2 mult
        double cab2 = exa2 * repxab12 + eya2 * repyab12; //1 add, 2 mult
        double cab3 = exa3 * repxab13 + eya3 * repyab13; //1 add, 2 mult
        double cab4 = exa4 * repxab14 + eya4 * repyab14; //1 add, 2 mult
        double cab5 = exa5 * repxab15 + eya5 * repyab15; //1 add, 2 mult
        double cab6 = exa6 * repxab16 + eya6 * repyab16; //1 add, 2 mult
        double cab7 = exa7 * repxab17 + eya7 * repyab17; //1 add, 2 mult

        double cba0 = exb0 * repxba10 + eyb0 * repyba10; //1 add, 2 mult
        double cba1 = exb0 * repxba11 + eyb0 * repyba11; //1 add, 2 mult
        double cba2 = exb0 * repxba12 + eyb0 * repyba12; //1 add, 2 mult
        double cba3 = exb0 * repxba13 + eyb0 * repyba13; //1 add, 2 mult
        double cba4 = exb0 * repxba14 + eyb0 * repyba14; //1 add, 2 mult
        double cba5 = exb0 * repxba15 + eyb0 * repyba15; //1 add, 2 mult
        double cba6 = exb0 * repxba16 + eyb0 * repyba16; //1 add, 2 mult
        double cba7 = exb0 * repxba17 + eyb0 * repyba17; //1 add, 2 mult

        double thab0 = sqrt(repxab10 * repxab10 + repyab10 * repyab10) * cospsi;
        double thab1 = sqrt(repxab11 * repxab11 + repyab11 * repyab11) * cospsi;
        double thab2 = sqrt(repxab12 * repxab12 + repyab12 * repyab12) * cospsi;
        double thab3 = sqrt(repxab13 * repxab13 + repyab13 * repyab13) * cospsi;
        double thab4 = sqrt(repxab14 * repxab14 + repyab14 * repyab14) * cospsi;
        double thab5 = sqrt(repxab15 * repxab15 + repyab15 * repyab15) * cospsi;
        double thab6 = sqrt(repxab16 * repxab16 + repyab16 * repyab16) * cospsi;
        double thab7 = sqrt(repxab17 * repxab17 + repyab17 * repyab17) * cospsi;

        double thba0 = sqrt(repxba10 * repxba10 + repyba10 * repyba10) * cospsi;
        double thba1 = sqrt(repxba11 * repxba11 + repyba11 * repyba11) * cospsi;
        double thba2 = sqrt(repxba12 * repxba12 + repyba12 * repyba12) * cospsi;
        double thba3 = sqrt(repxba13 * repxba13 + repyba13 * repyba13) * cospsi;
        double thba4 = sqrt(repxba14 * repxba14 + repyba14 * repyba14) * cospsi;
        double thba5 = sqrt(repxba15 * repxba15 + repyba15 * repyba15) * cospsi;
        double thba6 = sqrt(repxba16 * repxba16 + repyba16 * repyba16) * cospsi;
        double thba7 = sqrt(repxba17 * repxba17 + repyba17 * repyba17) * cospsi;

        double wab0 = -cab0 >= thab0 ? 1 : INFLUENCE;
        double wab1 = -cab1 >= thab1 ? 1 : INFLUENCE;
        double wab2 = -cab2 >= thab2 ? 1 : INFLUENCE;
        double wab3 = -cab3 >= thab3 ? 1 : INFLUENCE;
        double wab4 = -cab4 >= thab4 ? 1 : INFLUENCE;
        double wab5 = -cab5 >= thab5 ? 1 : INFLUENCE;
        double wab6 = -cab6 >= thab6 ? 1 : INFLUENCE;
        double wab7 = -cab7 >= thab7 ? 1 : INFLUENCE;

        double wba0 = -cba0 >= thba0 ? 1 : INFLUENCE;
        double wba1 = -cba1 >= thba1 ? 1 : INFLUENCE;
        double wba2 = -cba2 >= thba2 ? 1 : INFLUENCE;
        double wba3 = -cba3 >= thba3 ? 1 : INFLUENCE;
        double wba4 = -cba4 >= thba4 ? 1 : INFLUENCE;
        double wba5 = -cba5 >= thba5 ? 1 : INFLUENCE;
        double wba6 = -cba6 >= thba6 ? 1 : INFLUENCE;
        double wba7 = -cba7 >= thba7 ? 1 : INFLUENCE;

        sfx0 += repxab10 * wab0;
        sfy0 += repyab10 * wab0;
        sfx1 += repxab11 * wab1;
        sfy1 += repyab11 * wab1;
        sfx2 += repxab12 * wab2;
        sfy2 += repyab12 * wab2;
        sfx3 += repxab13 * wab3;
        sfy3 += repyab13 * wab3;
        sfx4 += repxab14 * wab4;
        sfy4 += repyab14 * wab4;
        sfx5 += repxab15 * wab5;
        sfy5 += repyab15 * wab5;
        sfx6 += repxab16 * wab6;
        sfy6 += repyab16 * wab6;
        sfx7 += repxab17 * wab7;
        sfy7 += repyab17 * wab7;

        double sfxbeta0 = (repxba10 * wba0) + (repxba11 * wba1);
        double sfybeta0 = (repyba10 * wba0) + (repyba11 * wba1);
        double sfxbeta1 = (repxba12 * wba2) + (repxba13 * wba3);
        double sfybeta1 = (repyba12 * wba2) + (repyba13 * wba3);
        double sfxbeta2 = (repxba14 * wba4) + (repxba15 * wba5);
        double sfybeta2 = (repyba14 * wba4) + (repyba15 * wba5);
        double sfxbeta3 = (repxba16 * wba6) + (repxba17 * wba7);
        double sfybeta3 = (repyba16 * wba6) + (repyba17 * wba7);

        social_force[IndexX(j)] += (sfxbeta0 + sfxbeta1) + (sfxbeta2 + sfxbeta3);
        social_force[IndexY(j, n)] += (sfybeta0 + sfybeta1) + (sfybeta2 + sfybeta3);
      } // n-1 * (12 adds, 18 mults, 6 divs, 1 exp, 4 sqrts)

      /************************************************/
      //UPDATE ACCELERATION TERM
      /************************************************/
      // get actual velocity, desired direction, desired speed

      // compute velocity difference
      double vdx00 = dsv0 * exa0; // 1 mul, 1 flop
      double vdy00 = dsv0 * eya0; // 1 mul, 1 flop
      double vdx01 = dsv1 * exa1; // 1 mul, 1 flop
      double vdy01 = dsv1 * eya1; // 1 mul, 1 flop
      double vdx02 = dsv2 * exa2; // 1 mul, 1 flop
      double vdy02 = dsv2 * eya2; // 1 mul, 1 flop
      double vdx03 = dsv3 * exa3; // 1 mul, 1 flop
      double vdy03 = dsv3 * eya3; // 1 mul, 1 flop
      double vdx04 = dsv4 * exa4; // 1 mul, 1 flop
      double vdy04 = dsv4 * eya4; // 1 mul, 1 flop
      double vdx05 = dsv5 * exa5; // 1 mul, 1 flop
      double vdy05 = dsv5 * eya5; // 1 mul, 1 flop
      double vdx06 = dsv6 * exa6; // 1 mul, 1 flop
      double vdy06 = dsv6 * eya6; // 1 mul, 1 flop
      double vdx07 = dsv7 * exa7; // 1 mul, 1 flop
      double vdy07 = dsv7 * eya7; // 1 mul, 1 flop

      double vdx10 = vdx00 - avx0; // 1 add, 1 flop
      double vdy10 = vdy00 - avy0; // 1 add, 1 flop
      double vdx11 = vdx01 - avx1; // 1 add, 1 flop
      double vdy11 = vdy01 - avy1; // 1 add, 1 flop
      double vdx12 = vdx02 - avx2; // 1 add, 1 flop
      double vdy12 = vdy02 - avy2; // 1 add, 1 flop
      double vdx13 = vdx03 - avx3; // 1 add, 1 flop
      double vdy13 = vdy03 - avy3; // 1 add, 1 flop
      double vdx14 = vdx04 - avx4; // 1 add, 1 flop
      double vdy14 = vdy04 - avy4; // 1 add, 1 flop
      double vdx15 = vdx05 - avx5; // 1 add, 1 flop
      double vdy15 = vdy05 - avy5; // 1 add, 1 flop
      double vdx16 = vdx06 - avx6; // 1 add, 1 flop
      double vdy16 = vdy06 - avy6; // 1 add, 1 flop
      double vdx17 = vdx07 - avx7; // 1 add, 1 flop
      double vdy17 = vdy07 - avy7; // 1 add, 1 flop

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
      //n/8 * 
      // ( 689 add, 1106 mult, 280 div, 140 sqrt, 56 exp_fast +
      //   (n-(i+8))*(200 add, 320 mult, 80 div, 56 sqrt, 16 exp_fast) + 
      //   32 add, 32 mult 
      // )
    }
    for (int i = 0; i < n - 3; i += 8)
    {
      /************************************************/
      // LOADS
      /************************************************/
      double cx0 = position[IndexX(i)];
      double cy0 = position[IndexY(i, n)];
      double cx1 = position[IndexX(i + 1)];
      double cy1 = position[IndexY(i + 1, n)];
      double cx2 = position[IndexX(i + 2)];
      double cy2 = position[IndexY(i + 2, n)];
      double cx3 = position[IndexX(i + 3)];
      double cy3 = position[IndexY(i + 3, n)];
      double cx4 = position[IndexX(i + 4)];
      double cy4 = position[IndexY(i + 4, n)];
      double cx5 = position[IndexX(i + 5)];
      double cy5 = position[IndexY(i + 5, n)];
      double cx6 = position[IndexX(i + 6)];
      double cy6 = position[IndexY(i + 6, n)];
      double cx7 = position[IndexX(i + 7)];
      double cy7 = position[IndexY(i + 7, n)];

      double max0 = desired_max_speed[i];
      double max1 = desired_max_speed[i + 1];
      double max2 = desired_max_speed[i + 2];
      double max3 = desired_max_speed[i + 3];
      double max4 = desired_max_speed[i + 4];
      double max5 = desired_max_speed[i + 5];
      double max6 = desired_max_speed[i + 6];
      double max7 = desired_max_speed[i + 7];

      //compute prefered velocity by integrating over the social force for the timestep, assuming the social force is constant over \delta t
      double pvx0 = actual_velocity[IndexX(i)] + social_force[IndexX(i)] * TIMESTEP;               // 1 add, 1 mult => 2 flops
      double pvy0 = actual_velocity[IndexY(i, n)] + social_force[IndexY(i, n)] * TIMESTEP;         // 1 add, 1 mult => 2 flops
      double pvx1 = actual_velocity[IndexX(i + 1)] + social_force[IndexX(i + 1)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      double pvy1 = actual_velocity[IndexY(i + 1, n)] + social_force[IndexY(i + 1, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      double pvx2 = actual_velocity[IndexX(i + 2)] + social_force[IndexX(i + 2)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      double pvy2 = actual_velocity[IndexY(i + 2, n)] + social_force[IndexY(i + 2, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      double pvx3 = actual_velocity[IndexX(i + 3)] + social_force[IndexX(i + 3)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      double pvy3 = actual_velocity[IndexY(i + 3, n)] + social_force[IndexY(i + 3, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      double pvx4 = actual_velocity[IndexX(i + 4)] + social_force[IndexX(i + 4)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      double pvy4 = actual_velocity[IndexY(i + 4, n)] + social_force[IndexY(i + 4, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      double pvx5 = actual_velocity[IndexX(i + 5)] + social_force[IndexX(i + 5)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      double pvy5 = actual_velocity[IndexY(i + 5, n)] + social_force[IndexY(i + 5, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      double pvx6 = actual_velocity[IndexX(i + 6)] + social_force[IndexX(i + 6)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      double pvy6 = actual_velocity[IndexY(i + 6, n)] + social_force[IndexY(i + 6, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops
      double pvx7 = actual_velocity[IndexX(i + 7)] + social_force[IndexX(i + 7)] * TIMESTEP;       // 1 add, 1 mult => 2 flops
      double pvy7 = actual_velocity[IndexY(i + 7, n)] + social_force[IndexY(i + 7, n)] * TIMESTEP; // 1 add, 1 mult => 2 flops

      /************************************************/
      // UPDATE POSITION
      /************************************************/
      //compute the norm of the preferd velocity
      double xysq00 = pvx0 * pvx0; // 1 mult => 1 flop
      double xysq10 = pvx1 * pvx1; // 1 mult => 1 flop
      double xysq20 = pvx2 * pvx2; // 1 mult => 1 flop
      double xysq30 = pvx3 * pvx3; // 1 mult => 1 flop
      double xysq40 = pvx4 * pvx4; // 1 mult => 1 flop
      double xysq50 = pvx5 * pvx5; // 1 mult => 1 flop
      double xysq60 = pvx6 * pvx6; // 1 mult => 1 flop
      double xysq70 = pvx7 * pvx7; // 1 mult => 1 flop

      double xysq01 = xysq00 + (pvy0 * pvy0); // 1 add, 1 mult1 => 2 flops
      double xysq11 = xysq10 + (pvy1 * pvy1); // 1 add, 1 mult1 => 2 flops
      double xysq21 = xysq20 + (pvy2 * pvy2); // 1 add, 1 mult1 => 2 flops
      double xysq31 = xysq30 + (pvy3 * pvy3); // 1 add, 1 mult1 => 2 flops
      double xysq41 = xysq40 + (pvy4 * pvy4); // 1 add, 1 mult1 => 2 flops
      double xysq51 = xysq50 + (pvy5 * pvy5); // 1 add, 1 mult1 => 2 flops
      double xysq61 = xysq60 + (pvy6 * pvy6); // 1 add, 1 mult1 => 2 flops
      double xysq71 = xysq70 + (pvy7 * pvy7); // 1 add, 1 mult1 => 2 flops

      double nv0 = sqrt(xysq01); // 1 sqrt => 1 flops
      double nv1 = sqrt(xysq11); // 1 sqrt => 1 flops
      double nv2 = sqrt(xysq21); // 1 sqrt => 1 flops
      double nv3 = sqrt(xysq31); // 1 sqrt => 1 flops
      double nv4 = sqrt(xysq41); // 1 sqrt => 1 flops
      double nv5 = sqrt(xysq51); // 1 sqrt => 1 flops
      double nv6 = sqrt(xysq61); // 1 sqrt => 1 flops
      double nv7 = sqrt(xysq71); // 1 sqrt => 1 flops

      //formula 12 in the paper --> compute control_value according to norm
      double cv0 = nv0 > max0 ? (max0 / nv0) : 1.0; // 1 div => 1 flops
      double cv1 = nv1 > max1 ? (max1 / nv1) : 1.0; // 1 div => 1 flops
      double cv2 = nv2 > max2 ? (max2 / nv2) : 1.0; // 1 div => 1 flops
      double cv3 = nv3 > max3 ? (max3 / nv3) : 1.0; // 1 div => 1 flops
      double cv4 = nv4 > max4 ? (max4 / nv4) : 1.0; // 1 div => 1 flops
      double cv5 = nv5 > max5 ? (max5 / nv5) : 1.0; // 1 div => 1 flops
      double cv6 = nv6 > max6 ? (max6 / nv6) : 1.0; // 1 div => 1 flops
      double cv7 = nv7 > max7 ? (max7 / nv7) : 1.0; // 1 div => 1 flops

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
      double dx0 = final_destination[IndexX(i)] - cx0;        // 1 add => 1 flop
      double dy0 = final_destination[IndexY(i, n)] - cy0;     // 1 add => 1 flop
      double dx1 = final_destination[IndexX(i + 1)] - cx1;    // 1 add => 1 flop
      double dy1 = final_destination[IndexY(i + 1, n)] - cy1; // 1 add => 1 flop
      double dx2 = final_destination[IndexX(i + 2)] - cx2;    // 1 add => 1 flop
      double dy2 = final_destination[IndexY(i + 2, n)] - cy2; // 1 add => 1 flop
      double dx3 = final_destination[IndexX(i + 3)] - cx3;    // 1 add => 1 flop
      double dy3 = final_destination[IndexY(i + 3, n)] - cy3; // 1 add => 1 flop
      double dx4 = final_destination[IndexX(i + 4)] - cx4;    // 1 add => 1 flop
      double dy4 = final_destination[IndexY(i + 4, n)] - cy4; // 1 add => 1 flop
      double dx5 = final_destination[IndexX(i + 5)] - cx5;    // 1 add => 1 flop
      double dy5 = final_destination[IndexY(i + 5, n)] - cy5; // 1 add => 1 flop
      double dx6 = final_destination[IndexX(i + 6)] - cx6;    // 1 add => 1 flop
      double dy6 = final_destination[IndexY(i + 6, n)] - cy6; // 1 add => 1 flop
      double dx7 = final_destination[IndexX(i + 7)] - cx7;    // 1 add => 1 flop
      double dy7 = final_destination[IndexY(i + 7, n)] - cy7; // 1 add => 1 flop

      double d0_0 = dx0 * dx0; // 1 add, 2 mult => 3 flops
      double d1_0 = dx1 * dx1; // 1 add, 2 mult => 3 flops
      double d2_0 = dx2 * dx2; // 1 add, 2 mult => 3 flops
      double d3_0 = dx3 * dx3; // 1 add, 2 mult => 3 flops
      double d4_0 = dx4 * dx4; // 1 add, 2 mult => 3 flops
      double d5_0 = dx5 * dx5; // 1 add, 2 mult => 3 flops
      double d6_0 = dx6 * dx6; // 1 add, 2 mult => 3 flops
      double d7_0 = dx7 * dx7; // 1 add, 2 mult => 3 flops

      double d0_1 = d0_0 + dy0 * dy0; // 1 add, 2 mult => 3 flops
      double d1_1 = d1_0 + dy1 * dy1; // 1 add, 2 mult => 3 flops
      double d2_1 = d2_0 + dy2 * dy2; // 1 add, 2 mult => 3 flops
      double d3_1 = d3_0 + dy3 * dy3; // 1 add, 2 mult => 3 flops
      double d4_1 = d4_0 + dy4 * dy4; // 1 add, 2 mult => 3 flops
      double d5_1 = d5_0 + dy5 * dy5; // 1 add, 2 mult => 3 flops
      double d6_1 = d6_0 + dy6 * dy6; // 1 add, 2 mult => 3 flops
      double d7_1 = d7_0 + dy7 * dy7; // 1 add, 2 mult => 3 flops

      double n0 = sqrt(d0_1); // 1 sqrt => 1 flop
      double n1 = sqrt(d1_1); // 1 sqrt => 1 flop
      double n2 = sqrt(d2_1); // 1 sqrt => 1 flop
      double n3 = sqrt(d3_1); // 1 sqrt => 1 flop
      double n4 = sqrt(d4_1); // 1 sqrt => 1 flop
      double n5 = sqrt(d5_1); // 1 sqrt => 1 flop
      double n6 = sqrt(d6_1); // 1 sqrt => 1 flop
      double n7 = sqrt(d7_1); // 1 sqrt => 1 flop

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
      desired_direction[IndexX(i)] = dx0 / n0;        
      desired_direction[IndexY(i, n)] = dy0 / n0;     
      desired_direction[IndexX(i + 1)] = dx1 / n1;    
      desired_direction[IndexY(i + 1, n)] = dy1 / n1; 
      desired_direction[IndexX(i + 2)] = dx2 / n2;    
      desired_direction[IndexY(i + 2, n)] = dy2 / n2; 
      desired_direction[IndexX(i + 3)] = dx3 / n3;    
      desired_direction[IndexY(i + 3, n)] = dy3 / n3; 
      desired_direction[IndexX(i + 4)] = dx4 / n4;    
      desired_direction[IndexY(i + 4, n)] = dy4 / n4; 
      desired_direction[IndexX(i + 5)] = dx5 / n5;    
      desired_direction[IndexY(i + 5, n)] = dy5 / n5; 
      desired_direction[IndexX(i + 6)] = dx6 / n6;    
      desired_direction[IndexY(i + 6, n)] = dy6 / n6; 
      desired_direction[IndexX(i + 7)] = dx7 / n7;    
      desired_direction[IndexY(i + 7, n)] = dy7 / n7; 
      
      //n/8 * (64 add, 88 mult, 16 div, 16 sqrt)
    }
  }
}