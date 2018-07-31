/*****************************************************************************************
Title:   read_snapshot
Purpose: Does as the name says...
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"



/********************************************
                load_snapshot
********************************************/
void load_snapshot(void)
{
    // Well, this does as the name says

    FILE *fd         = NULL;
    int  blksize1    = 0;
    int  blksize2    = 0;
    int  i           = 0;
    int  j           = 0;
    int  k           = 0;
    int  l           = 0;
    int  n           = 0;
    int  npart_file  = 0;
    int  n_with_mass = 0;
    int  npartTot    = 0;
    enum fields blocknr;
    char snap_file_mult[256];
    PARTICLE_DATA *Temp_P = NULL;

    // Try opening the file
    if(!(fd = fopen(snap_file, "rb")))
    {
        // We're here either because the file doesn't exist
        // or because there are multiple files (i.e. snap_file.0)
        sprintf(snap_file_mult, "%s.0", snap_file);

        if(!(fd = fopen(snap_file_mult, "rb")))
        {
            sprintf(error_message_local, "Error, could not open snapshot for reading "\
                "on task %d!\n", thistask);
            exit_code_local = 0;
            return;
        }
    }

    // Read the header from the file
    if(!(my_fread(&blksize1, sizeof(int), 1, fd)))
    {
        sprintf(error_message_local, "Error reading blksize1 before header!\n");
        exit_code_local = 0;
        fclose(fd);
        return;    
    }

    if(!(my_fread(&header, sizeof(SNAP_HEADER), 1, fd)))
    {
        sprintf(error_message_local, "Error reading header!\n");
        exit_code_local = 0;
        fclose(fd);
        return;
    }

    if(!(my_fread(&blksize2, sizeof(int), 1, fd)))
    {
        sprintf(error_message_local, "Error reading blksize2 after header!\n");
        exit_code_local = 0;
        fclose(fd);
        return;
    }

    blocknr = HEADER;
    if(!(block_check(blocknr, blksize1, blksize2, header)))
    {
        sprintf(error_message_local, "Error block checking header!\n");
        exit_code_local = 0;
        fclose(fd);
        return;
    }

    // Get total number of particles in sim
    for(i = 0, npartTot = 0; i < 6; i++)
    {
        npartTot += header.npartTot[i];
    }

    // Get the number of "gas" particles in the sim
    if(header.npartTot[0] > 0)
    {
        numpart = header.npartTot[0];
    } 

    else
    {
        numpart = header.npartTot[1];
    }

    // Allocate memory for master particle data container
    if(!(P = calloc(numpart, sizeof(PARTICLE_DATA))))
    {
        sprintf(error_message_local, "Error, could not allocate memory for P\n");
        fclose(fd);
        exit_code_local = 0;
        return;
    }

    // Loop over every file in the snapshot
    for(i = 0; i < header.num_files; i++)
    {
        // If applicable, open the next file and read its header
        if(i > 0)
        {
            sprintf(snap_file_mult, "%s.%d", snap_file, i);

            if(!(fd = fopen(snap_file_mult, "rb")))
            {
                sprintf(error_message_local, "Error, could not open snapshot %d\n", i);
                exit_code_local = 0;
                return;
            }

            if(!(my_fread(&blksize1, sizeof(int), 1, fd)))
            {
                sprintf(error_message_local, "Error reading blksize1 before header!\n");
                exit_code_local = 0;
                fclose(fd);
                return;
            }

            if(!(my_fread(&header, sizeof(SNAP_HEADER), 1, fd)))
            {
                sprintf(error_message_local, "Error reading header!\n");
                exit_code_local = 0;
                fclose(fd);
                return;
            }

            if(!(my_fread(&blksize2, sizeof(int), 1, fd)))
            {
                sprintf(error_message_local, "Error reading blksize2 after header!\n");
                exit_code_local = 0;
                fclose(fd);
                return;
            }

            blocknr = HEADER;
            if(!(block_check(blocknr, blksize1, blksize2, header)))
            {
                sprintf(error_message_local, "Error block checking header!\n");
                exit_code_local = 0;
                fclose(fd);
                return;
            }
        }

        // Get the total number of particles in the current file
        for(j = 0, npart_file = 0; j < 6; j++)
        {
            npart_file += header.npart[j];
        }

        // Allocate memory for current file's particles
        if(!(Temp_P = calloc(npart_file, sizeof(PARTICLE_DATA))))
        {
            sprintf(error_message_local, "Error, could not allocate memory for "\
                "Temp_P in file %d\n", i);
            fclose(fd);
            exit_code_local = 0;
            return;
        }

        // Read Positions
        if(!(my_fread(&blksize1, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize1 before pos "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        for(j = 0; j < npart_file; j++)
        {
            if(!(my_fread(&Temp_P[j].Pos[0], sizeof(float), 3, fd)))
            {
                sprintf(error_message_local, "Error reading pos "\
                    "(file %d, part %d)!\n", i, j);
                exit_code_local = 0;
                fclose(fd);
                free(Temp_P);
                return;
            }
        }

        if(!(my_fread(&blksize2, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize2 after pos "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        blocknr = POS;
        if(!(block_check(blocknr, blksize1, blksize2, header)))
        {
            sprintf(error_message_local, "Error block checking pos (file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        // Read Velocities
        if(!(my_fread(&blksize1, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize1 before vel "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        for(j = 0; j < npart_file; j++)
        {
            if(!(my_fread(&Temp_P[j].Vel[0], sizeof(float), 3, fd)))
            {
                sprintf(error_message_local, "Error reading vel "\
                    "(file %d, part %d)!\n", i, j); 
                exit_code_local = 0;
                fclose(fd);
                free(Temp_P);
                return;
            }
        } 

        if(!(my_fread(&blksize2, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize2 after vel "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        blocknr = VEL;
        if(!(block_check(blocknr, blksize1, blksize2, header)))
        {
            sprintf(error_message_local, "Error block checking vel (file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        // Read Ids
        if(!(my_fread(&blksize1, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize1 before ids "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        for(j = 0; j < npart_file; j++)
        {
            if(!(my_fread(&Temp_P[j].id, sizeof(int), 1, fd)))
            {
                sprintf(error_message_local, "Error reading id "\
                    "(file %d, part %d)!\n", i, j);
                exit_code_local = 0;
                fclose(fd);
                free(Temp_P);
                return;
            }
        }

        if(!(my_fread(&blksize2, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize2 after ids "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        blocknr = IDS;
        if(!(block_check(blocknr, blksize1, blksize2, header)))
        {
            sprintf(error_message_local, "Error block checking ids (file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        // Read Masses
        n_with_mass = get_with_masses(header);
        if(n_with_mass > 0)
        {
            if(!(my_fread(&blksize1, sizeof(int), 1, fd)))
            {
                sprintf(error_message_local, "Error reading blksize1 before mass "\
                    "(file %d)!\n", i);
                exit_code_local = 0;
                fclose(fd);
                free(Temp_P);
                return;
            }
        }

        for(k = 0, l = 0; k < 6; k++)
        {
            for(j = 0; j < header.npart[k]; j++, l++)
            {
                // Assign type
                Temp_P[l].type = k;

                if(header.mass[k] == 0)
                {
                    if(!(my_fread(&Temp_P[l].Mass, sizeof(float), 1, fd)))
                    {
                        sprintf(error_message_local, "Error reading mass "\
                            "(file %d, part %d)!\n", i, j);
                        exit_code_local = 0;
                        fclose(fd);
                        free(Temp_P);
                        return;
                    }
                }
                else
                {
                    Temp_P[l].Mass = header.mass[k];
                }
            }
        }

        if(n_with_mass > 0)
        {
            if(!(my_fread(&blksize2, sizeof(int), 1, fd)))
            {
                sprintf(error_message_local, "Error reading blksize2 after mass "\
                    "(file %d)!\n", i);
                exit_code_local = 0;
                fclose(fd);
                free(Temp_P);
                return;
            }

            blocknr = MASS;
            if(!(block_check(blocknr, blksize1, blksize2, header)))
            {
                sprintf(error_message_local, "Error block checking mass "\
                    "(file %d)!\n", i);
                exit_code_local = 0;
                fclose(fd);
                free(Temp_P);
                return;
            }
        }

        // Read SPH properties. Need to be able to handle gas and dm
        // only snapshots.
        if(header.npartTot[0] > 0)
        {
            n = header.npart[0];
        }
        else if(header.npartTot[0] == 0)
        {
            n = header.npart[1];
        }

        // This should never happen, but should is the operative word there
        if(n <= 0)
        {
            sprintf(error_message_local, "Error, number of gas or dm particles is 0 "\
                "(file %d!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        // Read Temps
        if(!(my_fread(&blksize1, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize1 before temps "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        for(j = 0; j < n; j++)
        {
            if(!(my_fread(&Temp_P[j].temperature, sizeof(float), 1, fd)))
            {
                sprintf(error_message_local, "Error reading temp "\
                    "(file %d, part %d)!\n", i, j);
                exit_code_local = 0;
                fclose(fd);
                free(Temp_P);
                return;
            }

            if(isnan(Temp_P[j].temperature) == 1)
            {
                sprintf(error_message_local, "Error, temperature is a nan "\
                    "(file %d, part %d)!\n", i, j);
                exit_code_local = 0;
                free(Temp_P);
                fclose(fd);
                return;
            }

            if(isfinite(Temp_P[j].temperature) == 0)
            {
                sprintf(error_message_local, "Error, temperature is infinite "\
                    "(file %d, part %d)!\n", i, j);
                exit_code_local = 0;
                free(Temp_P);
                fclose(fd);
                return;
            }
        }

        if(!(my_fread(&blksize2, sizeof(float), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize2 after temps "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        blocknr = TEMP;
        if(!(block_check(blocknr, blksize1, blksize2, header)))
        {
            sprintf(error_message_local, "Error block checking temps (file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        // Read Densities
        if(!(my_fread(&blksize1, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize1 before rho "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        for(j = 0; j < n; j++)
        {
            if(!(my_fread(&Temp_P[j].Rho, sizeof(float), 1, fd)))
            {
                sprintf(error_message_local, "Error reading rho "\
                    "(file %d, part %d)!\n", i, j);
                exit_code_local = 0;
                fclose(fd);
                free(Temp_P);
                return;
            }
        }

        if(!(my_fread(&blksize2, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize2 after rho "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        blocknr = RHO;
        if(!(block_check(blocknr, blksize1, blksize2, header)))
        {
            sprintf(error_message_local, "Error block checking rho (file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        // Read smoothing lengths
        if(!(my_fread(&blksize1, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize1 before hsml "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        for(j = 0; j < n; j++)
        {
            if(!(my_fread(&Temp_P[j].Hsml, sizeof(float), 1, fd)))
            {
                sprintf(error_message_local, "Error reading hsml "\
                    "(file %d, part %d)!\n", i, j);
                exit_code_local = 0;
                fclose(fd);
                free(Temp_P);
                return;
            }
        }

        if(!(my_fread(&blksize2, sizeof(int), 1, fd)))
        {
            sprintf(error_message_local, "Error reading blksize2 after hsml "\
                "(file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        blocknr = HSML;
        if(!(block_check(blocknr, blksize1, blksize2, header)))
        {
            sprintf(error_message_local, "Error block checking hsml (file %d)!\n", i);
            exit_code_local = 0;
            fclose(fd);
            free(Temp_P);
            return;
        }

        // Close the file
        fclose(fd);

        // Copy data to master particle data struct
        if(!(copy_to_master(Temp_P, header, n)))
        {   
            sprintf(error_message_local, "Error when copying to master (file %d)!\n", i);
            free(Temp_P);
            exit_code_local = 0;
            return;
        }

        // Free memory for current file's particle data container
        free(Temp_P);
    }
}



/********************************************
              copy_to_master
********************************************/
int copy_to_master(PARTICLE_DATA *T, SNAP_HEADER h, int n)
{
    // This function copies the data from the current file's
    // particle data container into the master particle data
    // container. n holds the number of either gas or dm
    // particles in current file.

    int i;
    static int master_ind = 0;

    for(i = 0; i < n; i++, master_ind++)
    {
        P[master_ind].Pos[0]        = T[i].Pos[0];
        P[master_ind].Pos[1]        = T[i].Pos[1];
        P[master_ind].Pos[2]        = T[i].Pos[2];
        P[master_ind].Vel[0]        = T[i].Vel[0];
        P[master_ind].Vel[1]        = T[i].Vel[1];
        P[master_ind].Vel[2]        = T[i].Vel[2];
        P[master_ind].id            = T[i].id;
        P[master_ind].Mass          = T[i].Mass;
        P[master_ind].temperature   = T[i].temperature;
        P[master_ind].Rho           = T[i].Rho;
        P[master_ind].Hsml          = T[i].Hsml;
        P[master_ind].type          = T[i].type;

        // Set the neutral fraction to be negative
        P[master_ind].XHI = -1.0;

        // Set tot_vel_flag to be 0
        P[master_ind].tot_vel_flag = 0;
    }

    // Check our master index
    if(master_ind > numpart)
    {
        return 0;
    }

    return 1;
}



/********************************************
                  my_fread
********************************************/
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    // This function makes sure the correct number of bytes are read

    size_t nread;

    if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
        return 0;
    } 

    return nread;
}



/********************************************
                block_check
********************************************/
int block_check(enum fields blocknr, int blksize1, int blksize2, SNAP_HEADER h)
{
    // This function gets the size of the given block and ensures that it's size
    // matches the padding values surrounding the block

    int size = 0;

    // Check padding values
    if(blksize1 != blksize2)
    {
        sprintf(error_message_local, "Error, block paddings don't match\n");
        return 0;
    }

    // Make sure the padding matches the actual size of the block
    size = get_block_size(blocknr, h);

    if(blksize1 != size)
    {
        sprintf(error_message_local, "Paddings don't match actual block size!\n");
        return 0;
    }

    return 1;
}



/********************************************
               get_block_size
********************************************/
int get_block_size(enum fields blocknr, SNAP_HEADER h)
{
    // This function gets the actual block size to make sure that the correct
    // values are written to the padding.

    int bsize = 0;
    int i;
    int nmass = 0;

    switch(blocknr)
    {
        case HEADER:
            bsize = sizeof(SNAP_HEADER);
            break;

        case POS:
        case VEL:
            bsize = sizeof(float) * 3 * (h.npart[0] + h.npart[1] + h.npart[2] + 
                    h.npart[3] + h.npart[4] + h.npart[5]); 
            break;

        case IDS:
            bsize = sizeof(int) * (h.npart[0] + h.npart[1] + h.npart[2] + h.npart[3] +
                    h.npart[4] + h.npart[5]);
            break;

        case MASS:
            for(i = 0; i < 6; i++)
            {
                if(h.mass[i] == 0)
                {
                    nmass += h.npart[i];
                }
            }

            bsize = sizeof(float) * nmass;
            break;

        case TEMP:
        case RHO:
        case HSML:
            if(h.npartTot[0] == 0)
            {
                bsize = sizeof(float) * h.npart[1];
            }

            if(h.npartTot[0] > 0)
            {
                bsize = sizeof(float) * h.npart[0];
            }
            break;

        default:
            bsize = -1;
    }

    return bsize;
}



/********************************************
               get_with_mass
********************************************/
int get_with_masses(SNAP_HEADER h)
{
    // Returns the number of particles with variable masses
 
    int i;
    int npart = 0;

    for(i = 0; i < 6; i++)
    {
        if(h.mass[i] == 0)
        {
            npart += h.npart[i];
        }
    }

    return npart;
}
