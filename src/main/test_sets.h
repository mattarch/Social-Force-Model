#ifndef TEST_SET_H
#define TEST_SET_H

// acceleration from start straight
double acceleration_direction0[] = {1, 0, -1, 0};
double acceleration_expected0[]  = {2.68, 0, -2.68, 0};
double acceleration_vel0[]       = {0, 0, 0, 0};
int acceleration_n0           = 2;
// acceleration from start 45°
double acceleration_direction1[] = {1, 1, -1, 1};
double acceleration_expected1[]  = {2.68, 2.68, -2.68, 2.68};
double acceleration_vel1[]       = {0, 0, 0, 0};
int acceleration_n1           = 2;
// acceleration with velocity straight
double acceleration_direction2[] = {1, 0, -1, 0};
double acceleration_expected2[]  = {1.68, 0, -1.68, 0};
double acceleration_vel2[]       = {0.5, 0, -0.5, 0};
int acceleration_n2           = 2;
// acceleration with velocity 45°
double acceleration_direction3[] = {1, 1, -1, 1};
double acceleration_expected3[]  = {1.68, 2.68, -1.68, 2.68};
double acceleration_vel3[]       = {0.5, 0, -0.5, 0};
int acceleration_n3           = 2;
// acceleration to zero
double acceleration_direction4[] = {1, 0, -1, 0};
double acceleration_expected4[]  = {0, 0, 0, 0};
double acceleration_vel4[]       = {1.34, 0, -1.34, 0};
int acceleration_n4           = 2;
// deacceleration straight
double acceleration_direction5[] = {1, 0, -1, 0};
double acceleration_expected5[]  = {-0.8, 0, 0.8, 0};
double acceleration_vel5[]       = {1.74, 0, -1.74, 0};
int acceleration_n5              = 2;
// deacceleration 45°
double acceleration_direction6[] = {1, 1, -1, 1};
double acceleration_expected6[]  = {-0.8, 2.68, 0.8, 2.68};
double acceleration_vel6[]       = {1.74, 0, -1.74, 0};
int acceleration_n6              = 2;
/*
double repulsion_position0[] = {5, 3, 8, 2, 9, 3};
double repulsion_speed0[]    = {0.7, 0.5, 1.3};
double repulsion_direction0[]= {0, 1, 1, 0, 4, 3};
int repulsion_n0 = 3;
double repulsion_expected_result0[] = { 0, 0, -0.000020266568837, 0.000005901638773, -0.000000000000016, -0.000000000000004, 
                                        0.000033293428727, -0.000018112519707, 0, 0, -0.000000000027150, -0.000000000023879,
                                        0.000004591945972, -0.000000780381604, 0.034905007547125, 0.084268142615003, 0, 0};

double border_pos0[] = {0, 0, 10 , 10};
double border_borders0[] = {0, 10};
double border_expected0[] = {1, 0, 1, 1, 0, 0, 1, 0};
int border_n0 = 2;
int border_nb0 = 2;

double social_acc0[] = {3, 2, -5, 4};
double social_prep0[] = {0, 0, 1.5, -0.5, -1.5, 0.5, 0, 0};
double social_brep0[] = {1, 2, -3, 4, 5, 2, -2, -2};
double social_expected0[] = {2.5, 7.5, -3.5, 4.5};
int social_n0 = 2;
int social_nb0 = 2;
*/
double position_pos0[] = {0, 0, 10, 10};
double position_dir0[] = {1, 1, 0, 0};
double position_speed0[] = {1, 2};
double position_force0[] = {1, -1, 0, 5};
double position_vel0[] = {5, 5, 1, 1};
double position_expected0[] = {1, 1, 0, 0};
int position_n0 = 2;
#endif