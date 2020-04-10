/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow
*/

#ifndef SOCIAL_FORCE_H_ /* Include guard */
#define SOCIAL_FORCE_H_

// default values
#define AVG_SPEED 1.34 // in basic scenario this is the desired speed of every person
#define RELAX_TIME 0.5 // in [seconds]; used for smooth acceleration
#define N_BORDERS 2    // number of biarders used in the scenario
#define TIMESTEP 0.2   // in [seconds]; timestep for simulation

// parameters model PAGE 8
#define V_ALPHA_BETA 2.1 // in m^{2}s^{-2}
#define SIGMA 0.3        // in m
#define U_ALPHA_B 10.0   // in m^{2}s^{-2}
#define R 0.2            // in m
#define PSI 1.75         // in radians
#define INFLUENCE 0.5    // pure (?)

#ifdef DEBUG
#define CONSOLE_PRINT(x) printf x
#else
#define CONSOLE_PRINT(x) \
    do                   \
    {                    \
    } while (0)
#endif

void initialize_people(double *position, double *desired_direction, double *final_destination, double *desired_speed, int n);
void initialize_borders(double *borders, int n_borders);

#endif