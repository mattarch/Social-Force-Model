#ifndef TEST_SET_H
#define TEST_SET_H

// basic direction test case for 4 persons
double direction_position0[] = {0,0,10,10,0,0,10,10};
double direction_fdest0[] = {10,10,0,0,10,10,0,0};
double direction_expected0[] = {0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071};
int direction_n0 = 4;

double acceleration_desired_speed[] = {1.34, 1.34, 1.34, 1.34};
// acceleration from start straight and 45°
double acceleration_direction0[] = {1, 0, -1, 0, 1, 1, -1, 1};
double acceleration_expected0[] = {2.68, 0, -2.68, 0, 2.68, 2.68, -2.68, 2.68};
double acceleration_vel0[] = {0, 0, 0, 0, 0, 0, 0, 0};
int acceleration_n0 = 4;
// acceleration with velocity straight and 45°
double acceleration_direction1[] = {1, 0, -1, 0, 1, 1, -1, 1};
double acceleration_expected1[] = {1.68, 0, -1.68, 0, 1.68, 2.68, -1.68, 2.68};
double acceleration_vel1[] = {0.5, 0, -0.5, 0, 0.5, 0, -0.5, 0};
int acceleration_n1 = 4;
// acceleration to zero straight and 45°
double acceleration_direction2[] = {1, 0, -1, 0, 1, 1, -1, 1};
double acceleration_expected2[] = {0, 0, 0, 0, 0, 0, 0, 0};
double acceleration_vel2[] = {1.34, 0, -1.34, 0, 1.34, 1.34, -1.34, 1.34};
int acceleration_n2 = 4;
// deacceleration straight and 45°
double acceleration_direction3[] = {1, 0, -1, 0 ,1, 1, -1, 1};
double acceleration_expected3[] = {-0.8, 0, 0.8, 0, -0.8, 2.68, 0.8, 2.68};
double acceleration_vel3[] = {1.74, 0, -1.74, 0, 1.74, 0, -1.74, 0};
int acceleration_n3 = 4;

double repulsion_position0[] = {5, 3, 8, 2, 9, 3};
double repulsion_speed0[] = {0.7, 0.5, 1.3};
double repulsion_direction0[] = {0, 1, 1, 0, 4, 3};
int repulsion_n0 = 3;
double repulsion_expected_result0[] = {0, 0, -0.000020266568837, 0.000005901638773, -0.000000000000016, -0.000000000000004,
                                       0.000033293428727, -0.000018112519707, 0, 0, -0.000000000027150, -0.000000000023879,
                                       0.000004591945972, -0.000000780381604, 0.034905007547125, 0.084268142615003, 0, 0};

double border_pos0[] = {0, 0, 10, 10};
double border_borders0[] = {0, 10};
double border_expected0[] = {1, 0, 1, 1, 0, 0, 1, 0};
int border_n0 = 2;
int border_nb0 = 2;

// basic social force test for 4 persons and 2 walls
double social_acc0[] = {1, 2, -3, 5, 2, 2, -4, -2};
double social_prep0[] = {0, 0, 1, 1, 2, 3, -1, 1,
                        -1, -1, 0, 0, 2, 2, -2, 1,
                        -2, -3, -2, -2, 0, 0, 3, -1,
                        1, -1, 2, -1, -3, 1, 0, 0};
double social_brep0[] = {0, 1, 1, 1,
                        -1, 1, 0, 1,
                        0, 0, 1, 2,
                        2, 1, -1, 2};
double social_expected0[] = {4,9,-5,9,2,-2,-3,0};
int social_n0 = 4;
int social_nb0 = 2;

double position_desired_speed[] = {1.34, 1.34, 1.34, 1.34};
double position_dir[] = {0, 0, 0, 0, 0, 0, 0, 0}; // not needed for computation, only written to
double position_speed[] = {0, 0, 0, 0}; // not needed for computation, only written to
double position_pos[] = {0, 0, 0, 0, 0, 0, 0, 0}; // starting at origin
// Test0, 4 unit directions, 0 starting speed, unit social force
double position_force0[] = {1, 0, 0, 1, -1, 0, 0, -1}; // unit social force
double position_vel0[] = {0, 0, 0, 0, 0, 0, 0, 0};
double position_expected0[] = {0.04, 0, 0, 0.04, -0.04, 0, 0, -0.04};
int position_n0 = 4;
// Test1, 4 unit directions, 1.742 starting speed, unit social force
double position_force1[] = {1, 0, 0, 1, -1, 0, 0, -1}; // unit social force
double position_vel1[] = {1.742, 0, 0, 1.742, -1.742, 0, 0, -1.742};
double position_expected1[] = {0.3484, 0, 0, 0.3484, -0.3484, 0, 0, -0.3484};
int position_n1 = 4;
#endif