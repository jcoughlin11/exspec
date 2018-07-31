/*****************************************************************************************
Title:   los_redshift.c
Purpose: Contains functions used to calculate the redshift range spanned by the los
Notes:   1.) This is based on the method described in the Trident paper (Hummels17):
                https://arxiv.org/abs/1612.03935

         2.) Uses the gsl to do root finding (as seen in this example): 
                https://www.gnu.org/software/gsl/doc/html/roots.html

         3.) Uses the gsl to do integration
https://www.gnu.org/software/gsl/manual/html_node/Numerical-integration-examples.html

         4.) This is so much easier in python
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "allvars.h"
#include "proto.h"



/********************************************
              redshift_range
********************************************/
float redshift_range(float L)
{
    // This is the main driver function for calculating the redshift range spanned by
    // the los. The basic idea is this: we have L = D * \int_b^a f(z) dz. We know
    // f(z), L, a, and D. We want b. In order to do this we set up:
    // g(x) = D * \int_x^a f(z)dz - L and then solve for x such that g(x) = 0.
    // We root find, in other words. The gsl is used to do both the root finding and
    // the integration of the expansion factor 1 / E(z). All of the cosmological info
    // is taken from the snapshot header, which is a global variable. We also assume,
    // following Trident, that z_A = z_snapshot. The quasar is at z_A and the observer
    // is at z_B. As such, z_A > z_B. Note, though that my limits of integration are
    // opposite those of trident. Trident is set up such that L is negative. I believe the
    // gsl requires lower_lim < upper_lim, so I've swapped them such that L > 0. The
    // negative sign is properly handled when getting z_i later on in los.c

    // Parameters:
    // -----------
    //      L : float
    //          This is the length of the LOS in comoving GUL

    // Returns:
    // --------
    //      z_B : float
    //          The redshift corresponding to the endpoint of the los (the observer)

    int iters                      = 0;
    int max_iters                  = 10000;
    int status                     = 0;
    double lower_lim               = 0.0;
    double upper_lim               = 0.0;
    double z_lo                    = 0.0;
    double z_hi                    = 0.0;
    double hubble_dist             = 0.0;
    const gsl_root_fsolver_type *T = NULL;
    gsl_root_fsolver *s            = NULL;
    gsl_function F;
    gsl_integration_workspace *w   = NULL;
    Integration_Params params;

    // Get the Hubble distance. This is defined as D_H = c / H0. speed_of_light is in
    // km/s, and H0 = 100 * h is in km/s/Mpc, so c / H0 has units of proper Mpc. I'm
    // going to convert L from comoving GUL to proper Mpc.

    // Convert L from comoving GUL to proper GUL
    L = header.time * L;

    // Convert from GUL to Mpc. This is just dividing by 1000 and h.
    L = L / (1000.0 * header.HubbleParam);

    // Get the Hubble distance
    hubble_dist = speed_of_light / (100.0 * header.HubbleParam);

    // Set up the initial search interval. Trident has the integration limits reversed
    // from what they are in Hogg99, which results in Trident's L being negative. This
    // means that their lower_limit of integration is actually z_A. I'm reversing it
    // here in order to get a positive answer because, I think, the gsl requires the
    // lower limit to actually be lower than the upper limit. The negative sign is 
    // correctly handled when actually getting dz_i, though. As such, my lower_limit 
    // of integration is what I'm trying to find (since that's z_B).
    upper_lim = header.redshift;
    z_lo = upper_lim  * 0.01;
    z_hi = upper_lim + 0.1;

    // Set up integrator
    w = gsl_integration_workspace_alloc(1000);

    // Set up params to be passed to my_func
    params.upper_lim = upper_lim;
    params.answer = L;
    params.w = w;
    params.hubble_dist = hubble_dist;

    // Set up the root solver. Use Brent's method. For the guess, the gsl is a little
    // different than python. You have to give it bounds to work within.
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    F.function = &redshift_roots;
    F.params = &params;
    gsl_root_fsolver_set(s, &F, z_lo, z_hi);

    // Loop until we find the right upper limit
    do
    {
        // Update the state of the root solver
        status = gsl_root_fsolver_iterate(s);

        // Perform the root finding algorithm
        lower_lim = gsl_root_fsolver_root(s);

        // Update the interval? Not sure what this does, but it
        // doesn't work unless this is done
        status = gsl_root_test_interval(z_lo, z_hi, 0, 0.001);

        // Update iterations to prevent infinite loop
        iters++;
    
    } while((status == GSL_CONTINUE) && (iters < max_iters));

    // Clean up
    gsl_root_fsolver_free(s);
    gsl_integration_workspace_free(w);

    // Check to make sure it worked
    if(lower_lim > upper_lim)
    {
        sprintf(error_message_local, "Error, lower_lim > upper_lim "\
            "(task %d)!\n", thistask);
        exit_code_local = 0;
        return (lower_lim = -1.0);
    }

    if(isnan(lower_lim) == 1)
    {
        sprintf(error_message_local, "Error, lower_lim is a nan "\
            "(task %d)!\n", thistask);
        exit_code_local = 0;
        return (lower_lim = -1.0);
    }

    if(isfinite(lower_lim) == 0)
    {
        sprintf(error_message_local, "Error, lower_lim is infinite "\
            "(task %d)!\n", thistask);
        exit_code_local = 0;
        return (lower_lim = -1.0);
    }

    return lower_lim;
}



/********************************************
             redshift_roots
********************************************/
double redshift_roots(double lower_lim, void *params)
{
    // As described in the redshift_range doc string, what we want is the root of
    // g(x) = \int_low^upper f(z)dz - L. This function is g(x) that the gsl will
    // use to do the actual root finding. This function relies on the gsl to evaluate the
    // integral of f(z)dz.

    // Parameters:
    // -----------
    //      lower_lim : double
    //          This is z_B, which is the redshift of the observer, and is the lowre limit
    //          of integration.

    //      params : void *
    //          This is a pointer to an array holding all of the extra parameters needed.
    //          This in includes the upper_limit of integration and the value of L.

    // Returns:
    // --------
    //      g(lower_limit) : double
    //          This is the value of the function g(x) defined above evaluated at the
    //          currently passed value of lower_limit.

    gsl_function F;
    gsl_integration_workspace *w = NULL;
    double result                = 0.0;
    double err                   = 0.0;
    double upper_lim             = 0.0;
    double answer                = 0.0;
    double hubble_dist           = 0.0;
    double g                     = 0.0;
    Integration_Params *p        = NULL;

    // Set up
    p = (Integration_Params *)params;
    upper_lim = p->upper_lim;
    answer = p->answer;
    w = p->w;
    hubble_dist = p->hubble_dist;
    F.function = &integrand;

    // Do the integration
    gsl_integration_qags(&F, lower_lim, upper_lim, 0, 1e-7, 1000, w, &result, &err);

    // Get the value of g(x) for this lower limit. Both L and hubble_dist should have
    // units of proper Mpc.
    g = (hubble_dist * result) - answer;

    return g;
}



/********************************************
                 integrand
********************************************/
double integrand(double z, void *params)
{
    // The function whose root we are trying to find g(x) = \int_x^upper f(z)dz - L has
    // an integral in it. This function is f(z). In particular, f(z) = 1/E(z), where
    // E(z) is the expansion factor = H(z) / H0. This function is integrated by the gsl
    // in redshift_roots().

    // Parameters:
    // -----------
    //      z : double
    //          The value of the redshift that the expansion factor should be evaluated at

    //      params : void *
    //          Required by the gsl. Contains extra parameters to the function. Not used.

    // Returns:
    // --------
    //      1 / E(z) : double
    //          The value of H0 / H(z), aka the inverse of the expansion factor. This 
    //          depends on the dark energy model being used.

    double result;

    result = 1.0 / expansion_factor(z);

    return result;
}
