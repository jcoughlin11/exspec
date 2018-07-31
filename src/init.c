/*****************************************************************************************
Title:   init.c
Purpose: Does an args check, reads in the parameter file, and sets up the rng, mpi data
            type, log file, units, number of los per processor, and the output directory  
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/********************************************
                    init
********************************************/
void init(int nargs, char **args)
{
    // This function does an args check and
    // reads in the parameter file

    FILE   *fd           = NULL;
    int    elements      = 0;
    size_t seed          = 0.0;
    char   *error        = NULL;
    char   buffer[256];
    char   buffer2[256];
    char   command[256];
    char   snap_base[256];

    // Open the logging file
    if(thistask == 0)
    {
        if(!(log_file = fopen("exspec_logging_file.txt", "w")))
        {
            printf("Error, could not open logging file for writing on task %d!\n", \
                thistask);
            exit_code_local = 0;
            write_error_to_log = 0;
            return;
        }

        // Write compile time options to log file
        log_write_compile_options();
    }

    // Allocate memory for the local processor's error message
    if(!(error_message_local = calloc(err_msg_size, sizeof(char))))
    {
        printf("Error, could not allocate memory for error_message_local on ");
        printf("task %d!\n", thistask);
        exit_code_local = 0;
        write_error_to_log = 0;
        return;
    }
    

    // Args check
    if(nargs != 2)
    {
        // This error will trigger for every processor, so I only need the error
        // message once
        if(thistask == 0)
        {
            sprintf(error_message_local, "Usage: exspec param_file.param\n");
            exit_code_local = 0;
        }
        return;
    }

    // Open the parameter file for reading
    if(!(fd = fopen(args[1], "r")))
    {
        sprintf(error_message_local, "Error, could not open parameter file for reading "\
            "on task %d!\n", thistask);
        exit_code_local = 0;
        return;
    }

    // Read the parameter file

    // Snapshot
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "Snapshot%s", snap_file);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }

    // Snapshot base (just the snapshot name, not including the path). This is used
    // for setting up the output directories because string manipulation in c is
    // awful and this is easy
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "Snapshot Base%s", snap_base);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }

    // HM file
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "HM_File%s", hm_file);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }

    // Dark param (if 1, using w0, wa. If 0, using Lambda)
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "Dark_Param%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    dark_param = atoi(buffer2);

    // Number of LOS
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "Nlos%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    nlos = atoi(buffer2);

    // Number of pixels per LOS
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "Npixels%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    npixels = atoi(buffer2);

    // Baryon Fraction (only matters in dm only sim)
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "Baryon_Frac%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    baryon_fraction = atof(buffer2);

    // GUL (Gadget Unit of Length in cm/h)
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "GUL%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    GUL_in_cm = atof(buffer2);

    // GUM (Gadget Unit of Mass in g/h)
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "GUM%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    GUM_in_g = atof(buffer2);

    // Normalization tolerance
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "Norm_tolerance%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    tolerance = atof(buffer2);
    
    // Maximum allowed iterations in normalization loop
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "max_norm_iters%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    max_iters = atof(buffer2);

    // W0 (only matters if de_flag == 1)
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "W0%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    w0 = atof(buffer2);

    // WA (only matters if de_flag == 1)
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "Wa%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    wa = atof(buffer2);

    // Hydrogen Mass Frac
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "H_Mass_Frac%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    hydrogen_mass_frac = atof(buffer2);

    // Helium Mass Frac
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "He_Mass_Frac%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    helium_mass_frac = atof(buffer2);

    // Mean molecular weight (this is actually just the coefficient. The mean
    // molecular weight mu is normally defined as mu = F * m_p, where m_p is
    // the mass of the proton. Here this is just F and m_p is explicitly added
    // where needed)
    error = fgets(buffer, sizeof(buffer), fd);
    elements = sscanf(buffer, "Molweight%s", buffer2);
    if(!(checkforerror(error, elements, &buffer[0])))
    {
        fclose(fd);
        return;
    }
    molweight = atof(buffer2);

    // Pixel width (comoving Mpc/h)
    #ifdef GENERAL_LOS
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "Delta_pix%s", buffer2);
        if(!(checkforerror(error, elements, &buffer[0])))
        {
            fclose(fd);
            return;
        }
        delta_pix = atof(buffer2);

        // Convert to comoving GUL. Since 1 GUL = 1 kpc/h,
        // I just need to multiply by 1000
        delta_pix *= 1000.0;
    #endif

    // Minimum number of pixels to use when truncating los. All that's needed here
    // is the denominator. That is, if I want to use a minimum of npixels / 3, I just
    // pass 3 for this argument. The default value is to use 1 / 3 of pixels.
    #ifdef FRAC_PIXELS
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "Min Pixel Frac for Truncating%s", buffer2);
        if(!(checkforerror(error, elements, &buffer[0])))
        {
            fclose(fd);
            return;
        }
        min_pixels_frac = 1.0 / atof(buffer2);

        if(min_pixels_frac > 1.0)
        {
            sprintf(error_message_local, "Error, min_pixels_frac just needs to be the "\
                "fraction denominator (i.e. 3 for 1 / 3 of pixels) (task %d)\n", \
                thistask);
            fclose(fd);
            exit_code_local = 0;
            return;
        }
    #endif
   
    // Close file
    fclose(fd);

    // Write parameter file options to log file (except for delta_pix, as it hasn't been
    // set yet if GENERAL_LOS isn't turned on)
    if(thistask == 0)
    {
        log_write_parameters();
    }

    // Check number of LOS. I currently don't support > 100000
    if(nlos >= 100000)
    {
        sprintf(error_message_local, "Error, no support for more than 100000 LOS! "\
            "(task %d)\n", thistask);
        exit_code_local = 0;
        return;
    }

    // Allocate memory for pixel struct only if GENERAL_LOS is not defined. This is
    // because, if GENERAL_LOS is turned on, the number of pixels varies from LOS to
    // LOS and so that memory allocation is done in the los() routine
    #ifndef GENERAL_LOS
        if(!(pixels = calloc(npixels, sizeof(PIXEL))))
        {
            sprintf(error_message_local, "Error, could not allocate memory for pixels! "\
                "on task %d\n", thistask);
            exit_code_local = 0;
            return;
        }
    #endif

    // Set the unit density for converting GUM / GUL^3 -> g / cm^3
    GUD_in_cgs = GUM_in_g / (GUL_in_cm * GUL_in_cm);
    GUD_in_cgs = GUD_in_cgs / GUL_in_cm;

    // Set up the rng for los selection so every processor uses a different seed
    seed = time(0) + thistask; 

    rng = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(rng, seed);

    // Set number of LOS to use to 1 if we're debugging
    #ifdef DEBUGGING
        nlos = 1;
    #endif

    #ifdef MPI_DEBUGGING
        nlos = 2;
    #endif

    start_index = 0;
    end_index = nlos;

    // Get the number of lines of sight for each processor. This is a hack. I'm just 
    // going to set the nlos to be the nearest number that is divisible by the number 
    // of processors that I'm using. Basically +/- a few los is not a big deal because
    // I don't feel like messing with a variable number of LOS per processor
    if(ntasks > 1)
    {
        if((nlos % ntasks) > 0)
        {
            nlos = nlos + (ntasks - (nlos % ntasks));
        }

        nlos_local = nlos / ntasks;

        // Now find the start and end indices
        start_index = thistask * nlos_local;
        end_index = start_index + nlos_local;

        // Create a custom MPI data type so we can easily send
        // the particle data to the other processors
        make_custom_mpi_type();
    }

    // Set up the timer
    if(thistask == 0)
    {
        timer.init_time = -1.0;
        timer.read_snapshot_time = -1.0;
        timer.gamma_time = -1.0;
        timer.spectra_time = -1.0;
        timer.run_time = -1.0;
        timer.is_instantiated = 0;
    }

    // Now create a directory for holding all of the data. The layout of this is:
    // exspec-snapshot_xxx/log_file and exspec-snapshot_xxx/specs/
    // Only need to do this once across all processors. The reason snap_base is used
    // instead of snap_file is because of the path. Stripping that out is a pain in c
    // and leaving it in leads to some hilariously named directory structures that,
    // while funny, are a pain to actually work with
    sprintf(parent_dir, "exspec-%s", snap_base);
    sprintf(spec_dir, "%s/specs", parent_dir);

    if(thistask == 0)
    {
        sprintf(command, "mkdir -p %s/specs", parent_dir);
        system(command);
    }
}



/********************************************
               checkforerror
********************************************/
int checkforerror(char *error, int elements, char *buf)
{
    // This function makes sure the proper number of elements,
    // as well as the proper element, was read from the param file

    char buffer2[256];

    if(error == NULL)
    {
        sprintf(error_message_local, "Error, could not read from parameter file "\
            "on task %d.\n", thistask);
        exit_code_local = 0;
        return 0;
    }

    if(elements == 0)
    {
        strncpy(buffer2, buf, strlen(buf));
      
        // Add null terminating
        buffer2[strlen(buf)-1]='\0';
        sprintf(error_message_local, "Could not read '%s' on task %d.\n", \
            buffer2, thistask);
        exit_code_local = 0;
        return 0;
    }

    return 1;
}



/********************************************
           make_custom_mpi_type
********************************************/
void make_custom_mpi_type(void)
{
    // See https://www.rc.colorado.edu/sites/default/files/Datatypes.pdf 
    // I need to make two: one for the header and one for the particle data
    // This makes sending the header and particle structure between processors
    // much, much easier.

    // Header type variables
    // h_blocks contains the number of elements that each header member contains
    // and h_types contains the data type of each header member
    int h_blocks[14] = {6,6,1,1,1,1,6,1,1,1,1,1,1,96};
    MPI_Datatype h_types[14] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT,\
                               MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE,\
                               MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR};
    MPI_Aint h_disp[14];
    MPI_Aint lb;

    // Particle type variables
    int p_blocks[11] = {3,3,1,1,1,1,1,1,1,1,1};
    MPI_Datatype p_types[11] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,\
                                MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT, MPI_INT, \
                                MPI_INT};
    MPI_Aint p_disp[11];
    MPI_Aint charex, intex, floatex, doublex;

    // Get the extents
    MPI_Type_get_extent(MPI_CHAR, &lb, &charex);
    MPI_Type_get_extent(MPI_INT, &lb, &intex);
    MPI_Type_get_extent(MPI_FLOAT, &lb, &floatex);
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &doublex);

    // Now fill the header displacements
    h_disp[0]  = (MPI_Aint)0;
    h_disp[1]  = 6 * intex;
    h_disp[2]  = h_disp[1]  + (6 * doublex);
    h_disp[3]  = h_disp[2]  + doublex;
    h_disp[4]  = h_disp[3]  + doublex;
    h_disp[5]  = h_disp[4]  + intex;
    h_disp[6]  = h_disp[5]  + intex;
    h_disp[7]  = h_disp[6]  + (6 * intex);
    h_disp[8]  = h_disp[7]  + intex;
    h_disp[9]  = h_disp[8]  + intex;
    h_disp[10] = h_disp[9]  + doublex;
    h_disp[11] = h_disp[10] + doublex;
    h_disp[12] = h_disp[11] + doublex;
    h_disp[13] = h_disp[12] + doublex;

    // Make MPI header structure
    MPI_Type_create_struct(14, h_blocks, h_disp, h_types, &mpi_header_type);
    MPI_Type_commit(&mpi_header_type);

    // Fill the particle displacements
    p_disp[0]  = (MPI_Aint)0;
    p_disp[1]  = 3 * floatex;
    p_disp[2]  = p_disp[1] + (3 * floatex);
    p_disp[3]  = p_disp[2] + floatex;
    p_disp[4]  = p_disp[3] + floatex;
    p_disp[5]  = p_disp[4] + floatex;
    p_disp[6]  = p_disp[5] + floatex;
    p_disp[7]  = p_disp[6] + floatex;
    p_disp[8]  = p_disp[7] + floatex;
    p_disp[9]  = p_disp[8] + intex;
    p_disp[10] = p_disp[9] + intex;

    // Make the MPI particle type
    MPI_Type_create_struct(11, p_blocks, p_disp, p_types, &mpi_particle_type);
    MPI_Type_commit(&mpi_particle_type);
}
