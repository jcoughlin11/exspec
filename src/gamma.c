/*****************************************************************************************
Title:   gamma.c
Purpose: Contains functions related to calculating
         the HI photoionization rate from the hm file
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "allvars.h"
#include "proto.h"



/********************************************
           photoionization_rate
********************************************/
void photoionization_rate(void)
{
    // This function calculates gamma, the photoionization rate of neutral hydrogen
    // due to the cosmic UV background that is due to galaxies and quasars. The
    // core component here is the HM (Haardt-Madau) file, which contains J(nu, z),
    // which is the UV background spectrum. This goes into an integral given in
    // Osterbrock 1989 equation 2.1 (left-hand side)

    // Driver routine for getting gamma
    double *zarray           = NULL;
    double *lambda_array     = NULL;
    double **J               = NULL;
    double **gamma_integrand = NULL;
    double *gamma_array      = NULL;
    int num_lambdas          = 0;
    int nionlambdas          = 0;
    int num_redshifts        = 1;
    long fmarker             = 0;
    int i;

    // Get the HM redshifts
    if(!(zarray = get_HM_redshifts(&fmarker, &num_redshifts)))
    {
        return;
    }

    // Read in the wavelengths. Returns them in cm
    if(!(lambda_array = get_HM_lambdas(&fmarker, &num_lambdas, &nionlambdas)))
    {
        free(zarray);
        return;
    }

    // Read in J
    if(!(J = get_HM_J(&fmarker, num_lambdas, num_redshifts)))
    {
        free(zarray);
        free(lambda_array);
        return;
    }

    // Tabulate the integrand of Bertone05 eq. A3
    // using lambda instead of nu
    if(!(gamma_integrand = get_gamma_integrand(J, lambda_array, nionlambdas, \
        num_redshifts)))
    {
        free(zarray);
        free(lambda_array);

        for(i = 0; i < num_lambdas; i++)
        {
            free(J[i]);
        }
        free(J);
        return;
    }

    // Do the integral in Bertone05 eq. A3
    if(!(gamma_array = integrate_gamma(gamma_integrand, lambda_array, nionlambdas, \
        num_redshifts)))
    {
        free(zarray);
        free(lambda_array);

        for(i = 0; i < num_lambdas; i++)
        {
            free(J[i]);
        }
        free(J);

        for(i = 0; i < nionlambdas; i++)
        {
            free(gamma_integrand[i]);
        }
        free(gamma_integrand);
        return;
    }

    // Now interpolate to get the right value of gamma at the redshift of the snapshot
    if((GammaXHI = interp_gamma(gamma_array, zarray, num_redshifts)) < 0.0)
    {
        free(zarray);
        free(lambda_array);

        for(i = 0; i < num_lambdas; i++)
        {
            free(J[i]);
        }
        free(J);

        for(i = 0; i < nionlambdas; i++)
        {
            free(gamma_integrand[i]);
        }
        free(gamma_integrand);
        return; 
    }

    if(thistask == 0)
    {
        printf("Gamma = %g s^-1\n", GammaXHI);
    }
 
    // Free memory
    free(zarray);
    free(lambda_array);

    for(i = 0; i < num_lambdas; i++)
    {
        free(J[i]);
    }

    free(J);

    for(i = 0; i < nionlambdas; i++)
    {
        free(gamma_integrand[i]);
    }

    free(gamma_integrand);
    free(gamma_array);
}



/********************************************
              get_HM_redshifts
********************************************/
double *get_HM_redshifts(long *fmarker, int *nz)
{
    // Reads the redshifts from the HM table. num_z starts at 1 b/c I check for z's 
    // based on =, and if that's where we start, then we miss the first one. nz is 
    // the number of redshifts contained in the file and fmarker is used to save 
    // the location in the file where the data actually starts.

    FILE   *fd   = NULL;
    double *z    = NULL;
    int    i     = 0;
    int    j     = 0;
    int    r     = 0;
    int    num_z = *nz;

    // Open the file
    if(!(fd = fopen(hm_file, "r")))
    {
        sprintf(error_message_local, "Error, could not open HM table for getting "\
            "redshifts!\n");
        exit_code_local = 0;
        return z;
    }

    // The redshifts start with the first '=' in the file, we read until then
    while((i = getc(fd)) != '=')
    {
        continue;
    }

    // Mark our location in the file so we can return
    if((*fmarker = ftell(fd)) == -1)
    {
        sprintf(error_message_local, "Error getting file position for redshifts!\n");
        fclose(fd);
        exit_code_local = 0;
        return z;
    }

    // Now get the number of redshifts
    while((i = getc(fd)) != '\n')
    {
        if(i == '=')
        {
            num_z++;
        }
    }

    // Allocate memory for z
    if(!(z = calloc(num_z, sizeof(double))))
    {
        sprintf(error_message_local, "Error, could not allocate memory for redshifts!\n");
        fclose(fd);
        exit_code_local = 0;
        return z;
    }

    // Go back to the first redshift. The minus puts us before 
    // the = instead of on it so we don't miss the first z
    fseek(fd, *fmarker - 1, SEEK_SET);

    // Read the redshifts
    while((i = getc(fd)) != '\n')
    {
        if(i == '=')
        {
            r = fscanf(fd, "%lf", &z[j]);

            if(r == EOF)
            {
                sprintf(error_message_local, "File read error for redshifts!\n");
                free(z);
                z = NULL;
                fclose(fd);
                exit_code_local = 0;
                return z;
            }

            j++;
        }
    }

    // Now mark the position of the end of the redshift line
    if((*fmarker = ftell(fd)) == -1)
    {
        sprintf(error_message_local, "Error getting file pos for end of z line!\n");
        fclose(fd);
        free(z);
        z = NULL;
        exit_code_local = 0;
        return z;
    }

    // Close the file
    fclose(fd);

    // Update num_redshifts
    *nz = num_z;

    return z;
}



/********************************************
               get_HM_lambdas
********************************************/
double *get_HM_lambdas(long *fmarker, int *num_lambdas, int *nionlambdas)
{
    // Reads the wavelengths from the HM table. fmarker is used to indicate the position
    // in the file where the data actually starts. num_lambdas is the number of
    // wavelengths contained in the file, and nionlambdas is the number of wavelengths
    // that have the appropriate energy to actually ionized neutral hydrogen

    FILE *fd = NULL;
    int i          = 0;
    int j          = 0;
    int r          = 0;
    int nl         = *num_lambdas;
    int nil        = *nionlambdas;
    double *lambda = NULL;

    // Open the file
    if(!(fd = fopen(hm_file, "r")))
    {
        sprintf(error_message_local, "Error, couldn't open HM table for reading "\
            "lambdas!\n");
        exit_code_local = 0;
        return lambda;
    }

    // Jump to the end of the redshifts line. The + puts us after the
    // \n at end of redshift line
    fseek(fd, *fmarker + 1, SEEK_SET);

    // Get number of wavelengths
    while((i = getc(fd)) != EOF)
    {
        if(i == '\n')
        {
            nl++;
        }
    }

    // Allocate memory for lambdas
    if(!(lambda = calloc(nl, sizeof(double))))
    {
        sprintf(error_message_local, "Error, could not allocate memory for HM "\
            "lambdas!\n");
        exit_code_local = 0;
        fclose(fd);
        return lambda;
    }

    // Go back to where lambdas start (no + so the first thing we read is
    // the first wavelength)
    fseek(fd, *fmarker, SEEK_SET);

    // Read in the wavelengths (in Angstroms). The file is laid out as:
    // wavelength J1 J2 ... JN\n
    for(j = 0; j < nl; j++)
    {
        r = fscanf(fd, "%lf", &lambda[j]);

        if(r == EOF)
        {
            sprintf(error_message_local, "File read error for wavelengths!\n");
            exit_code_local = 0;
            fclose(fd);
            free(lambda);
            lambda = NULL;
            return lambda;
        }

        // Read the rest of the line
        while((i = getc(fd)) != '\n')
        {
            continue;
        }
    }

    // Close the file
    fclose(fd);

    // Update num_lambdas
    *num_lambdas = nl;

    // Update nionlambdas
    for(i = 0; i < nl; i++)
    {
        // Convert from Angstroms to cm
        lambda[i] = lambda[i] * 1e-8;

        if(lambda[i] <= LYMANA)
        {
            nil++;
        }
    }

    *nionlambdas = nil;

    return lambda;
}



/********************************************
                  get_HM_J
********************************************/
double **get_HM_J(long *fmarker, int nl, int nz)
{
    // Actually reads J from the HM table (uv and J are the same)

    FILE *fd     = NULL;
    double **uv  = NULL;
    double dummy = 0.0;
    int i        = 0;
    int j        = 0;
    int r        = 0;

    // Open file
    if(!(fd = fopen(hm_file, "r")))
    {
        sprintf(error_message_local, "Error, could not open file for reading J!\n");
        exit_code_local = 0;
        return uv;
    }

    // Allocate memory for uv
    if(!(uv = calloc(nl, sizeof(double *))))
    {
        sprintf(error_message_local, "Error, could not allocate memory for J!\n");
        exit_code_local = 0;
        fclose(fd);
        return uv;
    }

    for(i = 0; i < nl; i++)
    {
        if(!(uv[i] = calloc(nz, sizeof(double))))
        {
            sprintf(error_message_local, "Error, could not allocate memeory for J!\n");
            exit_code_local = 0;
            fclose(fd);
            for(j = 0; j < i; j++)
            {
                free(uv[j]);
            }
            free(uv);
            uv = NULL;
            return uv;
        }
    }

    // Jump to the location in file where J data starts
    fseek(fd, *fmarker, SEEK_SET);

    // Read in the uv data: line layout is: wavelength J1 J2 ... JN\n
    for(i = 0; i < nl; i++)
    {
        // Read wavelength
        r = fscanf(fd, "%lf", &dummy);

        if(r == EOF)
        {
            sprintf(error_message_local, "File read error for dummy wavelength!\n");
            exit_code_local = 0;
            fclose(fd);
            for(j = 0; j < nl; j++)
            {
                free(uv[j]);
            }
            free(uv);
            uv = NULL;
            return uv;
        }
      
        // Read rest of line
        for(j = 0; j < nz; j++)
        {
            r = fscanf(fd, "%lf", &uv[i][j]);

            if(r == EOF)
            {
                sprintf(error_message_local, "File read error for rest of line!\n");
                exit_code_local = 0;
                fclose(fd);
                for(j = 0; j < nl; j++)
                {
                    free(uv[j]);
                }
                free(uv);
                uv = NULL;
                return uv;
            }
        }
    }

    // Close the file
    fclose(fd);

    return uv;
}



/********************************************
            get_gamma_integrand
********************************************/
double **get_gamma_integrand(double **J, double *lambda, int nionlambdas, int nredshifts)
{
    // Evaluates the integrand of Bertone05 eq. A3, but using lambda instead of nu

    int i          = 0;
    int j          = 0;
    double Zcharge = 1.0;       // Atomic number. 1 for H.
    double A0      = 6.3e-18;   // cm^2. A0 from Osterbrock89 eq. 2.4
    double epsilon = 0.0;       // The epsilon param in Osterbrock89 eq. 2.4
    double alpha[nionlambdas];
    double **gamma_integrand = NULL;

    // Allocate memory
    if(!(gamma_integrand = calloc(nionlambdas, sizeof(double *))))
    {
        sprintf(error_message_local, "Error, could not allocate memory for gamma "\
            "integrand!\n");
        exit_code_local = 0;
        return gamma_integrand;
    }

    for(i = 0; i < nionlambdas; i++)
    {
        if(!(gamma_integrand[i] = calloc(nredshifts, sizeof(double))))
        {
            sprintf(error_message_local, "Error, could not allocate memory for gamma "\
                "integrand!\n");
            exit_code_local = 0;

            for(j = 0; j < i; j++)
            {
                free(gamma_integrand[i]);
            }
            free(gamma_integrand);
            gamma_integrand = NULL;
            return gamma_integrand;
        }
    }

    // Calculate the values of alpha for those wavelengths <= LAMBDA0. Alpha is the cross
    // section of photoionization. This is equation 2.4 in Osterbrock's 1989 book
    // astrophysics of gaseous nebulae and active galactic nuclei.  
    for(i = 0; i < nionlambdas; i++)
    {
        epsilon = sqrt((LYMANA / lambda[i]) - 1.0);
        alpha[i] = (A0 / pow(Zcharge, 2.0)) * pow(lambda[i] / LYMANA, 4.0) * 
                    exp(4.0 - (4.0 * atan(epsilon)/epsilon)) / 
                    (1.0 - exp(-2.0 * M_PI / epsilon));
    }

    // Get integrand for every frequency at every redshift
    for(i = 0; i < nionlambdas; i++)
    {
        for(j = 0; j < nredshifts; j++)
        {
            gamma_integrand[i][j] = 4.0 * M_PI * J[i][j] * alpha[i] / 
                                    (PLANCK * lambda[i]);

            // Check for nans and infs
            if(isfinite(gamma_integrand[i][j]) == 0)
            {
                sprintf(error_message_local, "Error, got infinity for gamma "\
                    "integrand!\n");
                exit_code_local = 0;

                for(j = 0; j < nionlambdas; j++)
                {
                    free(gamma_integrand[j]);
                }

                free(gamma_integrand);

                gamma_integrand = NULL;

                return gamma_integrand;
            }

            if(isnan(gamma_integrand[i][j]) == 1)
            {
                sprintf(error_message_local, "Error, obtained a NaN for gamma "\
                    "integrand!\n");
                exit_code_local = 0;

                for(j = 0; j < nionlambdas; j++)
                {
                    free(gamma_integrand[j]);
                }

                free(gamma_integrand);

                gamma_integrand = NULL;

                return gamma_integrand;
            }
        }
    }

    return gamma_integrand;
}



/********************************************
              integrate_gamma
*********************************************/
double *integrate_gamma(double **gamma_integrand, double *lambda, int nionlambdas, int nz)
{
    // This function uses the gsl to do an interpolation integral on the
    // gamma_integrand to evaluate Bertone05 eq A3.

    int i               = 0;
    int k               = 0;
    double *gamma_array = NULL;
    double x[nionlambdas];
    double y[nionlambdas];
    gsl_interp_accel *accel = NULL;
    gsl_spline *spline      = NULL;

    // Allocate memory
    if(!(gamma_array = calloc(nz, sizeof(double))))
    {
        sprintf(error_message_local, "Error, could not allocate memory for gamma "\
            "array!\n");
        exit_code_local = 0;
        return gamma_array;
    }

    // Do the integral for every redshift in the HM table
    for(k = 0; k < nz; k++)
    {
        // Set the x and y data. The x data is lambda, and y is the integrand at
        // the redshift given by index k 
        for(i = 0; i < nionlambdas; i++)
        {
            y[i] = gamma_integrand[i][k];
            x[i] = lambda[i];
        }

        // Initialize the spline
        spline = gsl_spline_alloc(gsl_interp_cspline, nionlambdas);
        accel  = gsl_interp_accel_alloc();
        gsl_spline_init(spline, x, y, nionlambdas);

        // Actually perform the integral. Is from 0 to LAMBDA0
        gamma_array[k] = gsl_spline_eval_integ(spline, x[0], x[nionlambdas - 1], accel);

        // Check result for nans and infs
        if(isnan(gamma_array[k]) == 1)
        {
            sprintf(error_message_local, "Error, detected a nan in gamma array!\n");
            exit_code_local = 0;
            free(gamma_array);
            gamma_array = NULL;
            return gamma_array;
        }

        if(isfinite(gamma_array[k]) == 0)
        {
            sprintf(error_message_local, "Error, detected an inf in gamma array!\n");
            exit_code_local = 0;
            free(gamma_array);
            gamma_array = NULL;
            return gamma_array;
        }

        // Free the gsl resources for the next redshift
        gsl_spline_free(spline);
        gsl_interp_accel_free(accel);
    }   

    return gamma_array;
}



/********************************************
               interp_gamma
********************************************/
double interp_gamma(double *gamma_array, double *z, int nz)
{
    // This function uses the gsl to interpolate gamma_array to the value of gamma at
    // the redshift of the snapshot

    gsl_spline *spline      = NULL;
    gsl_interp_accel *accel = NULL;
    double photo_rate       = -1.0;

    // Set up the spline
    spline = gsl_spline_alloc(gsl_interp_cspline, nz);
    accel = gsl_interp_accel_alloc();
    gsl_spline_init(spline, z, gamma_array, nz);

    // Do the interpolation
    photo_rate = gsl_spline_eval(spline, header.redshift, accel);

    // Check for nans and infs
    if(isnan(photo_rate) == 1)
    {
        sprintf(error_message_local, "Error, gamma is a nan on task!\n");
        gsl_spline_free(spline);
        gsl_interp_accel_free(accel);
        exit_code_local = 0;
        return photo_rate;
    }

    if(isfinite(photo_rate) == 0)
    {
        sprintf(error_message_local, "Error, gamma is inf on task!\n");
        gsl_spline_free(spline);
        gsl_interp_accel_free(accel);
        exit_code_local = 0;
        return photo_rate;
    }

    // Free gsl resources
    gsl_spline_free(spline);
    gsl_interp_accel_free(accel);

    return photo_rate;
}
