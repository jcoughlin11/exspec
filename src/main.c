/*****************************************************************************************
Title:   exspec
Author:  Jared Coughlin
Date:    9/29/16
Purpose: Extract artificial Lyman alpha spectra from GADGET-2 snapshots
Notes:   1.) Snapshots must contain densities, temperatures, and smoothing lengths
         2.) Snapshot header must be the GADGET-2 standard 256 bytes (unless you want to
            change the code)
         3.) References:
            A.) Reading snapshot: read_snapshot.c packaged with vanilla GADGET-2
                (https://wwwmpa.mpa-garching.mpg.de/gadget/)

            B.) Photoionization rate:
                i.) Haardt-Madau 2005 UV background file (obtained from Cloudy10): 
                    https://arxiv.org/abs/astro-ph/9509093
                ii.) Equations: Osterbrock 1989 (Astrophysics of Gaseous Nebulae and 
                    Active Galactic Nucleii)

            C.) Dark Energy: 
                i.) Hubble equation: Amendola and Tsujikawa (Dark Energy: Theory
                    and Observations)
                ii.) w0-wa parameterization: Linder 2002 (Exploring the Expansion
                    History of the Universe, https://arxiv.org/abs/astro-ph/0208512)

            D.) LOS: This is all based on code originally written by S. Bertone and
                kindly provided via private communication (this includes the optimization
                of projecting the smoothing sphere onto the LOS in order to loop over
                only the possible range of overlapped pixels instead of the entire LOS)
                    i.) Trident paper: https://arxiv.org/abs/1612.03935 (this is used for
                        the LOS_REDSHIFT option)

            E.) Optical Depth (including normalization):
                i.) Theuns et al 1998 (P^3M-SPH Simulations of the Lyman Alpha Forest,
                    https://arxiv.org/abs/astro-ph/9805119)
                ii.) Bertone and White 2005 (How Do Galactic Winds Affect the Lyman
                    Alpha Forest?, https://arxiv.org/abs/astro-ph/0511028)
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



int main(int argc, char **argv)
{
    // Set up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &thistask);

    int    i            = 0;
    double start        = 0.0;
    double end          = 0.0;
    double sub_start    = 0.0;
    double sub_end      = 0.0;
    char   command[256];

    // Get start time
    start = clock();

    // Args check and read param file
    if(thistask == 0)
    {
        printf("Initializing...\n");
    }

    sub_start = clock();

    init(argc, argv);

    // Check for errors
    check_exit_code();

    sub_end = clock();

    // Save timing
    if(thistask == 0)
    {
        timer.init_time = (sub_end - sub_start) / (double)CLOCKS_PER_SEC;
        printf("Time to initialize: %lf s\n", timer.init_time);
    }

    // Read in header and snapshot data
    if(thistask == 0)
    {
        printf("Reading snapshot...\n");

        sub_start = clock();

        load_snapshot();

        sub_end = clock();

        // Save timing
        timer.read_snapshot_time = (sub_end - sub_start) / (double)CLOCKS_PER_SEC;
        printf("Time to read snapshot: %lf s\n", timer.read_snapshot_time); 
    }

    // Check for errors
    check_exit_code();

    // Distribute number of particles and header to other processors
    if(ntasks > 1)
    {
        MPI_Bcast(&numpart, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&header, 1, mpi_header_type, 0, MPI_COMM_WORLD);
    }

    // Get the photoionization rate from HM table and H
    if(thistask == 0)
    {
        printf("Calculating HI photoionization rate...\n");

        sub_start = clock();

        photoionization_rate();

        sub_end = clock();

        // Save timing
        timer.gamma_time = (sub_end - sub_start) / (double)CLOCKS_PER_SEC;
        printf("Time to get photoionization rate: %lf s\n", timer.gamma_time);

        // Get the Hubble parameter for snapshot redshift
        printf("Calculating Hubble parameter...\n");
        H_sim = get_hubble((1.0 / header.time) - 1.0);
    }

    // Check for errors
    check_exit_code();

    // Broadcast the photoionization rate and H to other processors
    if(ntasks > 1)
    {
        MPI_Bcast(&GammaXHI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&H_sim, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Set pixel width in COMOVING GUL. This is only true if
    // the LOS is parallel to an axis
    #ifndef GENERAL_LOS
        delta_pix = header.BoxSize / (float)npixels;
    #endif

    // Write delta_pix to logging file
    if(thistask == 0)
    {
        fprintf(log_file, "Pixel width: %.02lf\n\n", delta_pix); 
    }

    // Now we can extract spectra
    if(thistask == 0)
    {
        printf("Generating spectra...\n");
    }

    sub_start = clock();

    // Loop over the desired number of spectra to be extracted
    // Use start and end index for compatibility with mpi
    // If running in serial, these are just 0 and nlos
    for(i = start_index; i < end_index; i++)
    {
        // Set up the los
        los();

        // Get optical depths for pixels. This is really the heart
        // of the code.
        if(exit_code_local)
        {
            optical_depth();
        }

        // Normalize the spectrum
        if(exit_code_local)
        {
            normalize_spectrum();
        }

        // Write the spectrum to a file
        if(exit_code_local)
        {
            write_spectrum(i);
        }

        // Check for errors
        check_exit_code();

        // If using GENERAL_LOS, npixels changes from los to los, so the memory
        // needs to be freed before it is reallocated
        #ifdef GENERAL_LOS
            free(pixels);
            pixels = NULL;
        #endif
    }

    sub_end = clock();

    // Save timing
    if(thistask == 0)
    {
        timer.spectra_time = (sub_end - sub_start) / (double)CLOCKS_PER_SEC;
        printf("Time to create spectra: %lf s\n", timer.spectra_time); 
    }

    // Get end time and total run time
    end = clock();

    if(thistask == 0)
    {
        timer.run_time = (end - start) / (double)CLOCKS_PER_SEC;
        printf("\nTotal run time: %lfs\n", timer.run_time);
        log_write_timings();
    }

    // Free global variables
    free_globals();

    if(thistask == 0)
    {
        // Move the log file to the output directory
        sprintf(command, "mv exspec_logging_file.txt %s", parent_dir);
        system(command); 

        // Bid the user adieu :)
        printf("\nFinished. Have a nice day!\n");
    }
   
    MPI_Finalize();

    return 0;
}
