/*****************************************************************************************
Title:   write.c
Purpose: Contains functions related to, well, writing the spectrum to a file for plotting
             and analysis.
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"



/********************************************
               write_spectrum
********************************************/
void write_spectrum(int file_num)
{
    // Driver routine for writing the spectrum data to a file

    FILE *fd = NULL;
    int  i   = 0;
    char spec_file[256];
    char command[256];
    int  npix_to_use = 0;
    int  last        = 0;

    // Open file for writing
    if(nlos <= 999)
    {
        sprintf(spec_file, "%s-spectrum_%03d.txt", snap_file, file_num);
    }

    else if((nlos > 999) && (nlos < 10000))
    {
        sprintf(spec_file, "%s-spectrum_%04d.txt", snap_file, file_num);
    }

    else if((nlos >= 10000) && (nlos < 100000))
    {
        sprintf(spec_file, "%s-spectrum_%05d.txt", snap_file, file_num);
    }

    if(!(fd = fopen(spec_file, "w")))
    {
        sprintf(error_message_local, "Error, could not open file for writing spectrum "\
            "on task %d!\n", thistask);
        exit_code_local = 0;
        return;
    }

    #ifdef IRREG_LEN
        // gsl_rng_uniform_int generates integers between 0 and n - 1
        // inclusive with equal probability. I use npixels + 1 so that
        // npixels is a probability.
        // I need to enforce the number of pixels such that we have LOS that
        // are long enough to be useful. I have no idea what that number is, so
        // I'm going to guess and then mess around with it, I think.
        while(npix_to_use < (npixels * min_pixels_frac))
        {
            npix_to_use = gsl_rng_uniform_int(rng, npixels + 1);
        } 
    #else 
        npix_to_use = npixels;
    #endif

    // Write the file header
    fprintf(fd, "# Npixels: %d\n", npix_to_use);
    fprintf(fd, "# First pixel in kpc (x,y,z): %f \t %f \t %f\n", \
          pixels[0].Pos[0], pixels[0].Pos[1], pixels[0].Pos[2]);
    fprintf(fd, "# w0 = %.1f wa = %.1f\n", w0, wa);
    fprintf(fd, "# z = %.2f\n", header.redshift);
    fprintf(fd, "# Velocity(km/s) \t Wavelength(cm) \t nHI(cm_p^-3) \t T(K) \t "\
        "vHI(km/s) \t tau \t rho_tot(g/cm_p^3) \t flux\n");

    // Write the data
    // Prior to changing where the observer is (now they're at point B), starting at
    // pixel 0 here was fine. However, now that the observer is at point B and the
    // quasar at A (and pixel 0 is basically at point A), starting at point A and writing
    // npix_to_use will only save the portion of the LOS that starts at the quasar, not
    // the observer. npix_to_use should effectively switch where the quasar is, not where
    // the observer is, and since my v_H are zeroed out at the observer, I need to keep
    // those points

    // LOS:
    // A[pix0]-------------------[pix n -1]B

    // if using IRREG_LEN, npix_to_use != npixels, so, in general, we have:
    // A[pix0]-------------[npix_to_use]-----[pix n - 1]B.

    // By looping from pix 0 -> pix npix_to_use - 1, I'm keeping the pixels from A to
    // npix_to_use - 1, when it should be the other way around in order to include the
    // observer and zero point for v_H.

    // If IRREG_LEN is not on, then this still includes all of the pixels, but instead of
    // being written from 0 -> npix -1, they're written from npix-1 -> 0.

    last = npixels - 1;

    for(i = 0; i < npix_to_use; i++)
    {
        fprintf(fd, "%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e\n",\
        pixels[last - i].z_vel, pixels[last - i].lambda, pixels[last - i].nHI, \
        pixels[last - i].Temp_HI, pixels[last - i].Vel_HI, pixels[last - i].tau, \
        pixels[last - i].RhoTot, exp(-1.0 * pixels[last - i].tau)); 
    }

    // Close the file
    fclose(fd);

    // Now move the spectrum file to the spec dir
    sprintf(command, "mv %s %s", spec_file, spec_dir);
    system(command);
}
