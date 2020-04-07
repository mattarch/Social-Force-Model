#ifndef PARSE_ARGS_H_   /* Include guard */
#define PARSE_ARGS_H_
// inspired by: https://makework.blog/blog/2018/10/5/argument-parsing-in-c

#include <argp.h>
#include <stdbool.h>
#include <stdlib.h> 

char doc[] = "Social Force Model Basic Version";

enum options_with_no_short_option {
    OPTION_N_PEOPLE = 0x100,
    OPTION_N_TIMESTEPS = 0x1023,
};

static struct argp_option options[] = {
    { "n_people", OPTION_N_PEOPLE, "Int", 0, "Number of people used in simulation, default: 300." },
    { "n_timesteps", OPTION_N_TIMESTEPS, "Int", 0, "Number of timesteps, default: 300." },
    { "test", 't', 0, 0, "Run tests." },
    { 0 }
};

struct arguments {
    int n_people;
    int n_timesteps;
    bool test;
};

static error_t parse_option( int key, char *arg, struct argp_state *state ) {
    // state->input points to our user-defined argument structure
    struct arguments *arguments = state->input;
    char *p;

    switch ( key ) {
    case 't':
        arguments->test = true;
        break;
    case OPTION_N_PEOPLE:
        arguments->n_people = strtol(arg, &p, 10);
        break;
    case OPTION_N_TIMESTEPS:
        arguments->n_timesteps = strtol(arg, &p, 10);
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct argp argp = { options, parse_option, 0, doc };

#endif

