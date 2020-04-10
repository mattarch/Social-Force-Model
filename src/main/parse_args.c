/*
  Code of the ASL project.
  Topic: social force model
  Group: Simon Zurfluh, Matthias Matti, Tommaso Pegolotti, Lino Telschow

  This file provides the functionality to process command line input for a C application.
  code inspired by: https://makework.blog/blog/2018/10/5/argument-parsing-in-c
*/

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "parse_args.h"

/* string displayed when the help option is called.*/
char doc[] = "HELP: Social Force Model Basic Version";

/*
   Supported options
   - to add a new option do the following:
    (1) add a new option_key to the option_keys ENUM that is distinct from the other options
    (2) add option in this option array
    (3) add state to the arguments struct if necessary 
    (4) add a case for the switch statement using the defined option_key (from the ENUM), 
        in the function process_option that processes the option
*/
static struct arg_options options[] = {
    /* add options here */
    {"--n_people", OPTION_N_PEOPLE, 1, "300", "Number of people used in simulation."},
    {"--n_timesteps", OPTION_N_TIMESTEPS, 1, "300", "Number of timesteps."},
    {"--test", OPTION_TEST, 0, "false", "Run tests."},
    {"--width", OPTION_WALKWAY_WIDTH, 1, "4", "Width of the walkway."},
    {"--length", OPTION_WALKWAY_LENGTH, 1, "50", "Length of the walkway."},

    /* don't change last two entries*/
    {"--help", OPTION_HELP, 0, "-", "Displays this help."},
    {0}}; // do not delete the last entry, severs as flag for the end of the array

/*
  Same functionality as strsep in "string.h".
  Assumptions: none
  Parameters : none
*/
char *custom_strsep(char **stringp, const char *delimiter)
{
    char *start = *stringp;
    char *p;

    if (start != NULL)
    {
        p = strpbrk(start, delimiter);
    }
    else
    {
        p = NULL;
    }

    if (p == NULL)
    {
        *stringp = NULL;
    }
    else
    {
        *p = '\0';
        *stringp = p + 1;
    }

    return start;
}

/*
  This function prints the available command line options to the console upon the option --help.

  Assumptions: none
  Parameters : none
*/
static void display_help()
{
    printf("%s\n\n", doc);
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

/*
  This function parses the arguments input via command line and updates the state of the arguments.
  This function is called only from parse_args.

  Assumptions: option_key and struct arguments match
  Parameters :
            option_key : option_key defined in the ENUM option_keys
                   arg : contains the arguments if the option has any, otherwise the pointer is NULL
              arguments: struct of arguments
*/
static void process_option(int option_key, char *arg, struct arguments *arguments)
{
    char *p;

    switch (option_key)
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
    case OPTION_WALKWAY_WIDTH:
        arguments->walkway_width = strtol(arg, &p, 10); // TODO: error handling
        break;
    case OPTION_WALKWAY_LENGTH:
        arguments->walkway_length = strtol(arg, &p, 10); // TODO: error handling
        break;
    case OPTION_HELP:
        display_help();
        break;
    default:
        return;
    }
    return;
}

/*
  This function checks whether the parsed_string matches any option.
  If the parsed_string matches any option it returns the index of the option in the options[] array.
  Otherwise the return value is negative, indicating that the parsed_string does not match any option.
  This function is called only from parse_args.

  Assumptions: The last entry of struct arg_options options[] array is {0}.
  Parameters :
         parsed_string : parsed option string
*/
static size_t find_matching_option(char *parsed_string)
{
    size_t i = 0;
    while (strcmp(parsed_string, options[i].option_string) && strcmp("", options[i].option_string))
    {
        i++;
    }
    if (!strcmp("", options[i].option_string))
    {
        //
        return -1;
    }
    else
    {
        return i;
    }
}

/*
  This function parses the command line input, 
  checks if the inputs match any option and if all arguments are specified.
  Function succeeds if command line input is valid.

  Assumptions: - options which require an argument use exactly one '=' as a delimiter for the option string and the argument
               - options without delimiters do not use '=' in their string
  Parameters :
         argc : argument count passed from the main function
         argv : argument value[] array passed from the main function, argv[0] contains the program name
    arguments : struct of arguments
*/
void parse_args(int argc, char *argv[], struct arguments *arguments)
{
    int *option_key = (int *)calloc(argc, sizeof(int *));
    char **pointer_to_args = (char **)calloc(argc, sizeof(char **));

    bool parse_error = false;
    bool missing_arguments = false;

    /* parse command line input, check if inputs match options and if all arguments are specified */
    for (int i = 1; i < argc; i++)
    {
        char *option = argv[i];
        pointer_to_args[i] = argv[i];
        option = custom_strsep(&pointer_to_args[i], "=");

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
        /* command line input is not valid, some input does not match any option*/
        char *msg = argv[0];
        strcat(msg, " --help: no command line option found");
        puts(msg);
        exit(1);
    }
    else if (missing_arguments)
    {
        /* command line input is valid, but some arguments are missing for a specified option*/
        char *msg = argv[0];
        strcat(msg, " --help: missing some arguments for a specified option");
        puts(msg);
        exit(1);
    }

    /* command line input is valid, parse the arguments*/
    for (int i = 1; i < argc; i++)
    {
        process_option(option_key[i], pointer_to_args[i], arguments);
    }

    free(option_key);
    free(pointer_to_args);
}
