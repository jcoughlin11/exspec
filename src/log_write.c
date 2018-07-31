/*****************************************************************************************
Title: log_write.c
Purpose: Contains functions for writing information to the log file
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/********************************************
         log_write_compile_options
********************************************/
void log_write_compile_options(void)
{
    // This function writes the compile time options to the log file

    // Write the section heading
    fprintf(log_file, "****************************************\n");
    fprintf(log_file, "         Compile Time Options\n");
    fprintf(log_file, "****************************************\n");

    // Write the options

    // Debugging
    fprintf(log_file, "Debugging: ");
    #ifdef DEBUGGING
        fprintf(log_file, "on\n");
    #else
        fprintf(log_file, "off\n");
    #endif

    // MPI Debugging
    fprintf(log_file, "MPI Debugging: ");
    #ifdef MPI_DEBUGGING
        fprintf(log_file, "on\n");
    #else
        fprintf(log_file, "off\n");
    #endif

    // Irreg len
    fprintf(log_file, "Irregular Length LOS: ");
    #ifdef IRREG_LEN
        fprintf(log_file, "on\n");
    #else
        fprintf(log_file, "off\n");
    #endif

    // General LOS
    fprintf(log_file, "General LOS: ");
    #ifdef GENERAL_LOS
        fprintf(log_file, "on\n");
    #else
        fprintf(log_file, "off\n");
    #endif

    // Time Subroutines
    fprintf(log_file, "Time Subroutines: ");
    #ifdef TIME_SUBROUTINES
        fprintf(log_file, "on\n");
    #else
        fprintf(log_file, "off\n");
    #endif

    // LOS Redshift
    fprintf(log_file, "LOS Redshift Evolution: ");
    #ifdef LOS_REDSHIFT
        fprintf(log_file, "on\n");
    #else
        fprintf(log_file, "off\n");
    #endif

    // Frac pixels
    fprintf(log_file, "Use Custom Pixel Fraction: ");
    #ifdef FRAC_PIXELS
        fprintf(log_file, "on\n");
    #else
        fprintf(log_file, "off\n");
    #endif

    // Add some space between sections
    fprintf(log_file, "\n\n");
}



/********************************************
           log_write_parameters
********************************************/
void log_write_parameters(void)
{
    // This function writes the parameters passed in the parameter file to the log file

    // Write the section heading
    fprintf(log_file, "****************************************\n");
    fprintf(log_file, "              Parameters\n");
    fprintf(log_file, "****************************************\n");

    // Snapshot
    fprintf(log_file, "Input Snapshot: %s\n", snap_file);

    // HM file
    fprintf(log_file, "Haardt-Madau UV Background File: %s\n", hm_file);

    // Dark param
    fprintf(log_file, "Dark Param: %d\n", dark_param);

    // Number of LOS
    fprintf(log_file, "Number of LOS: %d\n", nlos);

    // Number of Pixels per LOS
    fprintf(log_file, "Pixels per LOS: %d\n", npixels);

    // Baryon Fraction
    fprintf(log_file, "Baryon Fraction: %.02g\n", baryon_fraction);

    // GUL
    fprintf(log_file, "Gadget Unit of Length in comoving cm/h: %.04g\n", GUL_in_cm);

    // GUM
    fprintf(log_file, "Gadget Unit of Mass in g/h: %.04g\n", GUM_in_g);

    // Normalization Tolerance
    fprintf(log_file, "Spectral Normalization Tolerance: %.02g\n", tolerance);

    // Max allowed iterations in normalization loop
    fprintf(log_file, "Max allowed iterations in normalization loop: %d\n", max_iters);

    // W0
    fprintf(log_file, "w0: %.02lf\n", w0);

    // Wa
    fprintf(log_file, "wa: %.02lf\n", wa);

    // Hydrogen mass frac
    fprintf(log_file, "Hydrogren Mass Fraction: %.02lf\n", hydrogen_mass_frac);

    // Helium mass frac
    fprintf(log_file, "Helium Mass Frac: %.02lf\n", helium_mass_frac);

    // Mean Molecular Weight
    fprintf(log_file, "Mean Molecular Weight: %.02lf\n", molweight);

    // Frac pixels
    #ifdef FRAC_PIXELS
        fprintf(log_file, "Minimum Fraction of Pixels Used: %.02lf\n", min_pixels_frac);
    #else
        fprintf(log_file, "Minimum Fraction of Pixels Used: NA\n");
    #endif

    // Number of processors
    fprintf(log_file, "Number of Processors Used: %d\n", ntasks);
}



/********************************************
             log_write_timings
********************************************/
void log_write_timings(void)
{
    // This function writes the subroutine timings to the log file

    if(!(timer.is_instantiated))
    {
        // Write the section heading
        fprintf(log_file, "****************************************\n");
        fprintf(log_file, "              Timings\n");
        fprintf(log_file, "****************************************\n");
    }

    // Instantiate so we don't re-draw the section heading (that's the only purpose of
    // this variable...)
    timer.is_instantiated = 1;

    // Write timings

    // Time to initialize
    if(timer.init_time > 0.0)
    {
        fprintf(log_file, "Time to initialize run: %lf seconds\n", timer.init_time);

        // Reset so this isn't rewritten every time
        timer.init_time = -1.0;
    }

    // Time to read snapshot
    if(timer.read_snapshot_time > 0.0)
    {
        fprintf(log_file, "Time to read snapshot run: %lf seconds\n", \
            timer.read_snapshot_time);

        // Reset so this isn't rewritten every time
        timer.read_snapshot_time = -1.0;
    }

    // Time get gamma
    if(timer.gamma_time > 0.0)
    {
        fprintf(log_file, "Time to calculate photoionization rate: %lf seconds\n", \
            timer.gamma_time);

        // Reset so this isn't rewritten every time
        timer.gamma_time = -1.0;
    }

    // Time to generate spectra
    if(timer.spectra_time > 0.0)
    {
        fprintf(log_file, "Time to generate spectra: %lf seconds\n", timer.spectra_time);

        // Reset so this isn't rewritten every time
        timer.spectra_time = -1.0;
    }

    // Total time to run
    if(timer.run_time > 0.0)
    {
        fprintf(log_file, "Total time to run: %lf seconds\n", timer.run_time);

        // Reset so this isn't rewritten every time
        timer.run_time = -1.0;
    }

    // Add some space between sections
    fprintf(log_file, "\n\n");
}



/********************************************
            log_write_errors
********************************************/
void log_write_errors(void)
{
    // This function collects any error messages that each processor might have and
    // writes them to the log file.

    char *error_message        = NULL;
    int  *exit_codes           = NULL;
    int  i                     = 0;
    int  root_exit_code        = 1;
    int  global_root_exit_code = 1;

    // Allocate memory for exit codes
    if(thistask == 0)
    {
        // Write the section heading
        fprintf(log_file, "\n");
        fprintf(log_file, "****************************************\n");
        fprintf(log_file, "              Errors\n");
        fprintf(log_file, "****************************************\n");

        if(!(exit_codes = calloc(ntasks, sizeof(int))))
        {
            fprintf(log_file, "Error, could not allocate memory for exit_codes!\n");
            root_exit_code = 0;
        }

        if(!(error_message = calloc(err_msg_size, sizeof(char))))
        {
            fprintf(log_file, "Error, could not allocate memory for error_message!\n");
            free(exit_codes);
            root_exit_code = 0;
        }
    }
    
    // Check for errors (can't call check_exit_code because that results in an infinite
    // loop. Also can't use exit_code_local and global because that product will give 0.
    // After all, that's how we got into this function in the first place.
    MPI_Allreduce(&global_root_exit_code, &root_exit_code, 1, MPI_INT, MPI_PROD, \
        MPI_COMM_WORLD);

    if(global_root_exit_code == 0)
    {
        return;
    }

    // Figure out which processors have error messages to send
    MPI_Gather(&exit_code_local, 1, MPI_INT, exit_codes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Now get the error messages from those processors that have them (if root)
    if(thistask == 0)
    {
        for(i = 0; i < ntasks; i++)
        {
            // If the root processor has an error, we don't want to call a Recv b/c
            // there won't be a matching send, so we just write the error to the file
            if(i == 0)
            {
                if(exit_codes[i] == 0)
                {
                    fprintf(log_file, "%s\n\n", error_message_local);
                }
            }

            else
            {
                if(exit_codes[i] == 0)
                {
                    // Wait for a message from this processor
                    MPI_Recv(error_message, err_msg_size, MPI_CHAR, i, MPI_ANY_TAG, \
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    // Write the message to the log file
                    fprintf(log_file, "%s\n\n", error_message);

                    // Reset the error message string
                    memset(error_message, 0, err_msg_size * sizeof(char));
                }
            }
        }

        // Free memory
        free(exit_codes);
        free(error_message);
    }

    // Have other processors send their error message to root
    else
    {
        if(exit_code_local == 0)
        {
            MPI_Send(error_message_local, err_msg_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
    }
}
