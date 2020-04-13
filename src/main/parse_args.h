/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow

  This file provides the functionality to process command line input for a C application.
  code inspired by: https://makework.blog/blog/2018/10/5/argument-parsing-in-c
*/

#ifndef PARSE_ARGS_H_ /* Include guard */
#define PARSE_ARGS_H_

#include <stdbool.h>

struct arguments
{
    int n_people;
    int n_timesteps;
    double walkway_width;
    double walkway_length;
    char* benchmark;
    int benchmark_max;
    bool test;
};

struct arg_options
{
    char option_string[24];
    int option_key;
    bool takes_arguments;
    char default_arg[20];
    char option_description[80];
};

enum option_keys
{
    OPTION_N_PEOPLE = 0x100,
    OPTION_N_TIMESTEPS = 0x101,
    OPTION_WALKWAY_WIDTH = 0x102,
    OPTION_WALKWAY_LENGTH = 0x103,
    OPTION_BENCHMARK = 0x104,
    OPTION_TEST = 0x146,
    OPTION_HELP = 0x147,
};

/* function declarations */
char *custom_strsep(char **stringp, const char *delimiter);
static void display_help();
static void process_option(int key, char *arg, struct arguments *arguments);
static size_t find_matching_option(char *parsed_string);
void parse_args(int argc, char *argv[], struct arguments *arguments);

#endif
