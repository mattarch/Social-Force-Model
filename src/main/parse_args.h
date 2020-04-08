#ifndef PARSE_ARGS_H_ /* Include guard */
#define PARSE_ARGS_H_

// inspired by: https://makework.blog/blog/2018/10/5/argument-parsing-in-c

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct arguments
{
    int n_people;
    int n_timesteps;
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
    OPTION_TEST = 0x102,
    OPTION_HELP = 0x42,
};

static void parse_option(int key, char *arg, struct arguments *arguments);
size_t find_matching_option(char *parsed_string);
void parse_args(int argc, char *argv[], struct arguments *arguments);

char doc[] = "HELP: Social Force Model Basic Version";

static struct arg_options options[] = {
    {"--n_people", OPTION_N_PEOPLE, 1, "300", "Number of people used in simulation."},
    {"--n_timesteps", OPTION_N_TIMESTEPS, 1, "300", "Number of timesteps."},
    {"--test", OPTION_TEST, 0, "false", "Run tests."},
    {"--help", OPTION_HELP, 0, "-", "Displays this help."},
    {0}};

static void display_help()
{
    printf("%s\n\n",doc);
    printf("%-20s%-20s%-100s\n", "option string", "default arg", "description");

    size_t i = 0;
    while (strcmp("", options[i].option_string))
    {
        char *option_string = options[i].option_string;
        if (options[i].takes_arguments)
        {
            strcat(option_string, "=arg");
        }
        printf("%-20s%-20s%-100s\n", option_string, options[i].default_arg, options[i].option_description);

        i++;
    }

    exit(1);
}

static void parse_option(int key, char *arg, struct arguments *arguments)
{
    char *p;

    switch (key)
    {
    case OPTION_N_PEOPLE:
        arguments->n_people = strtol(arg, &p, 10); // TODO: error handling
        break;
    case OPTION_N_TIMESTEPS:
        arguments->n_timesteps = strtol(arg, &p, 10); // TODO: error handling
        break;
    case OPTION_TEST:
        arguments->test = true;
        break;
    case OPTION_HELP:
        display_help();
        break;
    default:
        return ;
    }
    return;
}

size_t find_matching_option(char *parsed_string)
{
    size_t i = 0;
    while (strcmp(parsed_string, options[i].option_string) && strcmp("", options[i].option_string))
    {
        i++;
    }
    if (!strcmp("", options[i].option_string))
    {
        return -1;
    }
    else
    {
        return i;
    }
}

void parse_args(int argc, char *argv[], struct arguments *arguments)
{
    int *option_key = (int *)calloc(argc, sizeof(int *));
    char **pointer_to_args = (char **)calloc(argc, sizeof(char **));

    bool parse_error = false;
    bool missing_arguments = false;

    for (int i = 1; i < argc; i++)
    {
        char *option = argv[i];
        pointer_to_args[i] = argv[i];
        option = strsep(&pointer_to_args[i], "=");

        int option_index = find_matching_option(option);

        if (option_index < 0)
        {
            /* no command line option found */
            parse_error = true;
            break;
        }
        else if (options[option_index].takes_arguments && !pointer_to_args[i])
        {
            /* parsing of an option succeeded, but no argument specified*/
            missing_arguments = true;
        }
        else if (options[option_index].takes_arguments && pointer_to_args[i])
        {
            /* parsing of an option with argument succeeded */
            option_key[i] = options[option_index].option_key;
        }
        else
        {
            /* parsing of an option without argument succeeded */
            /* no argument needed*/
            option_key[i] = options[option_index].option_key;
        }
    }

    if (parse_error)
    {
        char *msg = argv[0];
        strcat(msg, " --help: no command line option found");
        puts(msg);
        exit(1);
    }
    else if (missing_arguments)
    {
        char *msg = argv[0];
        strcat(msg, " --help: missing some arguments");
        puts(msg);
        exit(1);
    }

    for (int i = 1; i < argc; i++)
    {
        parse_option(option_key[i], pointer_to_args[i], arguments);
    }

    free(option_key);
    free(pointer_to_args);
}
#endif
