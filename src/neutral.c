/*****************************************************************************************
Title:   neutral.c
Purpose: Contains functions related to getting the HI neutral fraction for a particle
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"



/********************************************
             get_neutral_frac
********************************************/
float get_neutral_frac(PARTICLE_DATA pd)
{
    // This function calculates the neutral hydrogen fraction for the 
    // particle. Follows eqs 6.2-6.7 in Bertone's thesis.

    float ne       = 0.0;
    float recrate  = 0.0;
    float collrate = 0.0;
    float XHI      = 0.0;

    // Get electron number density in proper cm^-3
    ne = get_electron_number_density(pd);

    // Get HII recombination rate in s^-1
    recrate = get_recombination_rate(pd);

    // Get HI collisional ionization rate proper cm^3 s^-1
    collrate = get_collisional_rate(pd);

    // Get HI neutral fraction dimensionless
    XHI = (recrate * ne) / (GammaXHI + (collrate * ne));

    // Check for nans and infs and negs
    if((isnan(XHI) == 1) || (isfinite(XHI) == 0) || (XHI < 0.0))
    {
        sprintf(error_message_local, "Error, XHI is either a nan, infinite, "\
            "or negative on task %d!\n", thistask);
        exit_code_local = 0;
    }

    return XHI;
}



/********************************************
       get_electron_number_density
********************************************/
float get_electron_number_density(PARTICLE_DATA pd)
{
    // This function evaluates Bertone thesis eq. 6.3.

    float rho = 0.0;
    float nH  = 0.0;
    float ne  = 0.0;

    // Convert from comoving GUL to proper GUL
    rho = pd.Rho / (header.time * header.time * header.time);

    // Take care of little h
    rho = rho * header.HubbleParam * header.HubbleParam;

    // Convert from proper GUL to proper cgs
    rho = rho * GUD_in_cgs;

    // Multiply by baryon frac if using dm only
    if(header.npartTot[0] == 0)
    {
        rho = rho * baryon_fraction;
    }

    // Get hydrogen number density (in proper cm^-3)
    nH = hydrogen_mass_frac * rho / hydrogen_mass;

    // Get electron number density
    ne = (2.0 - helium_mass_frac) / (2.0 * (1.0 - helium_mass_frac)) * nH;

   return ne;
}



/********************************************
          get_recombination_rate
********************************************/
float get_recombination_rate(PARTICLE_DATA pd)
{
    // Evaluates Bertone05 eq. A2

    float recomb = 0.0;
    float T4     = 0.0;

    T4 = pd.temperature / 1.0e4;

    if(pd.temperature <= 1.0e4)
    {
        recomb = 1.58e-13 * pow(T4, -0.51);
    }

    else
    {
        recomb = 1.58e-13 * pow(T4, -0.51 - (0.1 * log10(T4)));
    }

    // In s^-1
    return recomb;
}



/********************************************
           get_collisional_rate
********************************************/
float get_collisional_rate(PARTICLE_DATA pd)
{
    // Evaluates Bertone05 eq. A4

    float colrate = 0.0;
    float T5      = 0.0;

    T5 = pd.temperature / 1.0e5;

    colrate = 1.17e-10 * sqrt(pd.temperature) * 
                exp(-157809.1 / pd.temperature) / (1.0 + sqrt(T5));

    // Proper cm^3 / s
    return colrate;
}
