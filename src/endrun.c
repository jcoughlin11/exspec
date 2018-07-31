/*****************************************************************************************
Title: endrun.c
Purpose: Contains functions related to checking error codes and exiting as gracefully as
            I can.
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/********************************************
              check_exit_code
********************************************/
void check_exit_code(void)
{
    // This function just does an MPI_Reduce with multiplication on all of the
    // exit_code_local variables. If 0, then it frees all global variables that need to
    // be freed and then calls MPI_Finalize() and exits. The variable write_error is
    // used because sometimes there's an error where writing it to the log file is
    // not possible (such as if the log file itself failed to open, or memory for the
    // error_message_local array could not be allocated).

    int write_err_global = 1; 

    MPI_Allreduce(&exit_code_local, &exit_code_global, 1, MPI_INT, MPI_PROD, \
        MPI_COMM_WORLD);

    MPI_Allreduce(&write_error_to_log, &write_err_global, 1, MPI_INT, MPI_PROD, \
        MPI_COMM_WORLD);

    // If there's an error, we need to gracefully exit
    if(exit_code_global == 0)
    {
        // Check to see if the error is writable to the log file
        if(write_err_global)
        {
            // Write error to log file
            log_write_errors();
        }

        // Free all of the globals (all locals that need to be freed are freed in their
        // respective functions)
        free_globals();

        // Clean up mpi
        MPI_Finalize();

        // Exit
        exit(EXIT_FAILURE);
    }
}



/********************************************
                free_globals
********************************************/
void free_globals(void)
{
    // All of the global variables that are dynamic arrays are initialized to NULL in
    // allvars.c. This just goes over all of them and, if not NULL, then I know it has
    // been allocated and therefore needs to be freed.

    // Random number generator
    if(rng != NULL)
    {
        gsl_rng_free(rng);
    }

    // Particle data
    if(P != NULL)
    {
        free(P);
    }

    // Pixels
    if(pixels != NULL)
    {
        free(pixels);
    }

    // Log file
    if(log_file != NULL)
    {
        fclose(log_file);
    }

    // Error message string
    if(error_message_local != NULL)
    {
        free(error_message_local);
    }
}
