/** @file histo_cmdline.hpp
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.4
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef HISTO_CMDLINE_H
#define HISTO_CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef HISTO_CMDLINE_PACKAGE
/** @brief the program name (used for printing errors) */
#define HISTO_CMDLINE_PACKAGE "jellyfish histo"
#endif

#ifndef HISTO_CMDLINE_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define HISTO_CMDLINE_PACKAGE_NAME "jellyfish histo"
#endif

#ifndef HISTO_CMDLINE_VERSION
/** @brief the program version */
#define HISTO_CMDLINE_VERSION VERSION
#endif

/** @brief Where the command line options are stored */
struct histo_args
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  long buffer_size_arg;	/**< @brief Length in bytes of input buffer (default='10000000').  */
  char * buffer_size_orig;	/**< @brief Length in bytes of input buffer original value given at command line.  */
  const char *buffer_size_help; /**< @brief Length in bytes of input buffer help description.  */
  long low_arg;	/**< @brief Low count value of histogram (default='1').  */
  char * low_orig;	/**< @brief Low count value of histogram original value given at command line.  */
  const char *low_help; /**< @brief Low count value of histogram help description.  */
  long high_arg;	/**< @brief High count value of histogram (default='10000').  */
  char * high_orig;	/**< @brief High count value of histogram original value given at command line.  */
  const char *high_help; /**< @brief High count value of histogram help description.  */
  long increment_arg;	/**< @brief Increment value for buckets (default='1').  */
  char * increment_orig;	/**< @brief Increment value for buckets original value given at command line.  */
  const char *increment_help; /**< @brief Increment value for buckets help description.  */
  int threads_arg;	/**< @brief Number of threads (default='1').  */
  char * threads_orig;	/**< @brief Number of threads original value given at command line.  */
  const char *threads_help; /**< @brief Number of threads help description.  */
  char * output_arg;	/**< @brief Output file (default='/dev/fd/1').  */
  char * output_orig;	/**< @brief Output file original value given at command line.  */
  const char *output_help; /**< @brief Output file help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int buffer_size_given ;	/**< @brief Whether buffer-size was given.  */
  unsigned int low_given ;	/**< @brief Whether low was given.  */
  unsigned int high_given ;	/**< @brief Whether high was given.  */
  unsigned int increment_given ;	/**< @brief Whether increment was given.  */
  unsigned int threads_given ;	/**< @brief Whether threads was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct histo_cmdline_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure histo_args (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure histo_args (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *histo_args_purpose;
/** @brief the usage string of the program */
extern const char *histo_args_usage;
/** @brief all the lines making the help output */
extern const char *histo_args_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int histo_cmdline (int argc, char **argv,
  struct histo_args *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use histo_cmdline_ext() instead
 */
int histo_cmdline2 (int argc, char **argv,
  struct histo_args *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int histo_cmdline_ext (int argc, char **argv,
  struct histo_args *args_info,
  struct histo_cmdline_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int histo_cmdline_dump(FILE *outfile,
  struct histo_args *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int histo_cmdline_file_save(const char *filename,
  struct histo_args *args_info);

/**
 * Print the help
 */
void histo_cmdline_print_help(void);
/**
 * Print the version
 */
void histo_cmdline_print_version(void);

/**
 * Initializes all the fields a histo_cmdline_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void histo_cmdline_params_init(struct histo_cmdline_params *params);

/**
 * Allocates dynamically a histo_cmdline_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized histo_cmdline_params structure
 */
struct histo_cmdline_params *histo_cmdline_params_create(void);

/**
 * Initializes the passed histo_args structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void histo_cmdline_init (struct histo_args *args_info);
/**
 * Deallocates the string fields of the histo_args structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void histo_cmdline_free (struct histo_args *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int histo_cmdline_required (struct histo_args *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* HISTO_CMDLINE_H */