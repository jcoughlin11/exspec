/*****************************************************************************************
Title:   normal.c
Purpose: Normalizes the spectrum so the mean flux matches the mean observed flux at the 
            current redshift.
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"



/********************************************
            normalize_spectrum
********************************************/
void normalize_spectrum(void)
{
    // Normalizes the spectrum so that the mean simulated optical
    // depth matches the mean observed optical depth at the given
    // redshift. Uses method described in section 3.3 of Kim02
    // (there is no equation number for it)

    int   i      = 0;
    float tauobs = 0.0;
    float tausim = 0.0;
    float diff   = 0.0;
    float fac    = 0.0;
    int   iter   = 0;

    // Calculate the mean observed optical depth
    tauobs = 0.0032 * pow(1.0 + header.redshift, 3.37);

    // Get the mean simulated optical depth
    tausim = get_tausim();

    // Get difference between the values
    diff = tauobs - tausim;

    if(diff < 0.0)
    {
        diff = -1.0 * diff;
    }

    // Iterate until we're within the desired tolerance
    while(diff > tolerance)
    {
        iter++;
      
        if(iter > max_iters)
        {
            sprintf(error_message_local, "Error, max iterations for normalization "\
                "exceeded on task %d!\nlos_start = %.02lf \t %.02lf \t %.02lf\n"\
                "los_end = %.02lf \t %.02lf \t %.02lf\ndelta_pix = %.02lf\n", thistask,\
                los_start[0], los_start[1], los_start[2], los_end[0], los_end[1], \
                los_end[2], delta_pix);
            exit_code_local = 0;
            return;
        }

        // Get normalizing factor
        fac = tauobs / tausim;

        // Apply the normalizing factor
        for(i = 0; i < npixels; i++)
        {
            pixels[i].tau = pixels[i].tau * fac;

            // Check for nans and infs
            if(isnan(pixels[i].tau) == 1)
            {
                sprintf(error_message_local, "Error, pixels[%d].tau is a nan in normal.c"\
                    " on task %d!\npixels[0].Pos[0] = %lf\npixels[0].Pos[1] = %lf\n"\
                    "pixels[0].Pos[2] = %lf\nnpixels = %d\nlos_len = %lf\ndelta = %lf\n"\
                    "los_end[0] = %lf\nlos_end[1] = %lf\nlos_end[2] = %lf\niter = %d\n",\
                    i, thistask, pixels[0].Pos[0], pixels[0].Pos[1], pixels[0].Pos[2],\
                    npixels, los_len, delta_pix, los_end[0], los_end[1], los_end[2], \
                    iter);
                exit_code_local = 0;
                return;
            }

            if(isfinite(pixels[i].tau) == 0)
            {
                sprintf(error_message_local, "Error, pixels[%d].tau is infinite in "\
                    "normal.c on task %d!\npix[0].Pos[0] = %lf\npix[0].Pos[1] = %lf\n"\
                    "pixels[0].Pos[2] = %lf\nnpixels = %d\nlos_len = %lf\ndelta = %lf\n"\
                    "los_end[0] = %lf\nlos_end[1] = %lf\nlos_end[2] = %lf\niter = %d\n",\
                    i, thistask, pixels[0].Pos[0], pixels[0].Pos[1], pixels[0].Pos[2],\
                    npixels, los_len, delta_pix, los_end[0], los_end[1], los_end[2], \
                    iter);
                exit_code_local = 0;
                return;
            }
        }

        // Recalculate tausim
        tausim = get_tausim();

        // Find new diff
        diff = tauobs - tausim;

        if(diff < 0.0)
        {
            diff = -1.0 * diff;
        }
    }
}



/********************************************
                 get_tausim
********************************************/
float get_tausim(void)
{
    // Calculates the average optical depth in the simulated spectrum

    float sum = 0.0;
    float tausim;
    int i;

    for(i = 0; i < npixels; i++)
    {
        sum += exp(-1.0 * pixels[i].tau);
    }

    tausim = -1.0 * log(sum / (float)npixels);

    return tausim;
}
