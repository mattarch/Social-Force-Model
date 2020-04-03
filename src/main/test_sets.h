#ifndef TEST_SET_H
#define TEST_SET_H

double direction_position0[] = {0, 0, 10, 10};
double direction_fdest0[] = {10, 10, 0, 0};
int direction_n0 = 2;
double direction_expected0[] = {0.7071, 0.7071, -0.7071, -0.7071};
double direction_expected1[] = {0, 0.7071, -20, -0.7071};

double acceleration_direction0[] = {1, 1, 1, 1};
double acceleration_expected0[]  = {0, 0, 0, 0};
double acceleration_vel0[]       = {1, 0, 1, 0};
int acceleration_n0           = 2;

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

double social_acc0[] = {1, 0, 0, 0};
double social_prep0[] = {1, 1, 0, 0, 1, 1, 0, 0};
double social_brep0[] = {1, 1, 0, 0, 1, 1, 0, 0};
double social_expected0[] = {0, 1, 8, 3};
int social_n0 = 2;
int social_nb0 = 2;

double position_pos0[] = {0, 0, 10, 10};
double position_dir0[] = {1, 1, 0, 0};
double position_speed0[] = {1, 2};
double position_force0[] = {1, -1, 0, 5};
double position_vel0[] = {5, 5, 1, 1};
double position_expected0[] = {1, 1, 0, 0};
int position_n0 = 2;
#endif