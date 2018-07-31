/*****************************************************************************************
Title:   de.c
Purpose: Contains functions related to calculating the dark factor and Hubble parameter
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/********************************************
                 get_hubble
********************************************/
double get_hubble(double z)
{
    // This function calculated H(z) in s^-1

    double H0;
    double H;

    // Calculate H0 in km/s/Mpc
    H0 = 100.0 * header.HubbleParam;

    // Convert to s^-1
    H0 = H0 * 3.24e-20;

    H = expansion_factor(z); 

    H = H * H0;

    return H;
}



/********************************************
             expansion_factor
********************************************/
double expansion_factor(double z)
{
    // This function gets H(z) / H0. It ignores the contribution from radiation
    // since I'm not going back to high enough redshifts for that to matter.
    
    double E;

    E = ((1.0 - (header.Omega0 + header.OmegaLambda)) * pow(1.0 + z, 2.0)) +
        (header.Omega0 * pow(1.0 + z, 3.0)) + 
        (header.OmegaLambda * get_dark_factor(z));

    E = sqrt(E);

    return E; 
}



/********************************************
              get_dark_factor
********************************************/
double get_dark_factor(double z)
{
    // Gets the dark factor for Linder02 parameterization
    // See handwritten notes for derivation

    double dark_factor = 1.0;

    if(dark_param == 1)
    {
        dark_factor = pow(1. + z, 3.0 * (1.0 + w0 + wa)) * exp(-3.0 * wa * \
            (1.0 - (1.0 / (1.0 + z))));
    }

    return dark_factor;
}
