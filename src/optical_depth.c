/*****************************************************************************************
Title:   optical_depth.c
Purpose: Evaluates Eqs. B1,B2 and B3 in Bertone05
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"



/********************************************
               optical_depth
********************************************/
void optical_depth(void)
{
    // This is the driver routine for evaluating eqs.
    // B1,B2, and B3 and B6 in Bertone05

    int i             = 0;
    int j             = 0;
    int ngroups       = 0;
    int n_to_send     = 0;
    int nleft         = 0;
    int current_group = 0;
    int offset        = 0;
    float dist2       = 0.0;
    float t           = 0.0; // This is the parameter used in getting the perpendicular 
                             // distance I pass it to the function so I can "return" it 
                             // (I need it in the get total velocity function, and this 
                             // way I don't need to recalc it)
    PARTICLE_DATA *Temp_P = NULL;

    // If we have multiple processors, we need to distribute the particles from root to
    // other processors in chunks. The reason for this is because for large numpart,
    // MPI_Bcast cannot actually pass that much information at once. We also can't 
    // allocate memory to hold all of P at once on the other processors because we do 
    // not have enough RAM to do that.
    // So, I've set a buffer_size variable. Gadget uses 100MB, so that's what I'm going 
    // to use. Essentially, I'll send 100MB of particles to each processor at a time.    
    if(ntasks > 1)
    {
        // Check that the buffer size is large enough to hold at least one particle
        if(buffer_size < sizeof(PARTICLE_DATA))
        {
            if(thistask == 0)
            {
                sprintf(error_message_local, "Error, buffer size is not large enough to "\
                    "hold even one particle!\n");
                exit_code_local = 0;
            }

            // The return is outside the if statement so that all processors return
            return;
        }

        // Get the number of particles to send in each pass. This is the number of 
        // particles that will fit within the buffer, so n * sizeof(Pdata) = buffer. 
        // We have to floor it to ensure that it fits within the buffer. If we ceil, 
        // we would go over the buffer limit.
        n_to_send = floor((float)buffer_size / (float)sizeof(PARTICLE_DATA));

        // Now we have to get the number of passes that we need to do (the number of 
        // groups of n_to_send particles that we need in order to include all numpart 
        // particles)
        ngroups = ceil(((float)numpart * (float)sizeof(PARTICLE_DATA)) / \
            (float)buffer_size);

        #ifdef MPI_DEBUGGING
            printf("n_to_send = %d\n", n_to_send);
            printf("ngroups = %d\n", ngroups);
        #endif

        // Make sure we're sending enough particles. That is, 
        // n_to_send * ngroups >= numpart. If it's not, we need to add groups until the 
        // above condition is true
        while((n_to_send * ngroups) < numpart)
        {
            ngroups++;
        }

        // Initialize the number of particles left to be sent
        nleft = numpart;

        // Initialize the current group. We need to keep track of this because on the 
        // last group n_to_send will be different than what it is above if 
        // (n_to_send * ngroups) % numpart != 0 and it will need to be changed. It will
        //  be less than the above n_to_send, I believe. 
        current_group = 0;
    }

    else
    {  
        // We can just use P on root in the serial case
        ngroups = 1;
        n_to_send = numpart;
        nleft = numpart;
    }

    // Loop over every particle group
    for(current_group = 0; current_group < ngroups; current_group++)
    {
        // Set up P for other processors and send data, if necessary
        if(ntasks > 1)
        {
            // Get the point in P to start at for copying particles over. This is
            // is done before changing n_to_send because the offset depends on the
            // number of particles sent in the previous groups, not the number to be
            // sent in the current group. If there are 5 particles total and we send 2
            // per group with 3 groups, the first time we send particles 0 and 1, then
            //  2,3, then, on the third group, we want to start at 4. If we change 
            // n_to_send before getting offset, it would be 1 and 1 * (current_group = 2)
            //  = 2, not 4, so you start in the wrong place.
            offset = current_group * n_to_send;

            // If we're on the last group, we need to change n_to_send so we don't 
            // overshoot numpart
            if(current_group == (ngroups - 1))
            {
                n_to_send = nleft; 
            }

            // Allocate memory
            if(!(Temp_P = calloc(n_to_send, sizeof(PARTICLE_DATA))))
            {
                sprintf(error_message_local, "Error, could not allocate memory for "\
                    "Temp_P (task %d(\n", thistask);
                exit_code_local = 0;
                return;
            }

            // Now we need to copy over the particles from P to Temp_P on root
            if(thistask == 0)
            {
                for(i = offset, j = 0; j < n_to_send; i++, j++)
                {
                    // I think I can get away with a shallow copy. The members of P
                    // are never actually changed in any way in this code, so this 
                    // should be safe. 
                    Temp_P[j] = P[i];
                }
            }

            // Check for errors. This is because if someting goes wrong allocating memory
            // on root for Temp_P, we need all of the other processors to stop before the
            // bcast
            check_exit_code();

            // Broadcast the particles in the current group
            MPI_Bcast(&Temp_P[0], n_to_send, mpi_particle_type, 0, MPI_COMM_WORLD);
        }

        // Serial case
        else if(ntasks == 1)
        {
            Temp_P = &P[0];
        }
      
        // Update the number of particles that are left to transfer
        nleft = nleft - n_to_send; 

        // Loop over every particle in the group
        for(i = 0; i < n_to_send; i++)
        {
            // Get perpendicular distance between particle and LOS (returns d^2 as sqrt 
            // is slow) dist2 has units of comoving GUL^2
            dist2 = get_perpendicular_dist(Temp_P[i], &t);

            // Compare the distance to the particle's smoothing length (h has units of 
            // comoving GUL)
            if(dist2 <= (Temp_P[i].Hsml * Temp_P[i].Hsml))
            {
                // We're here if the particle overlaps the LOS, so we need to calculate
                // the neutral hydrogen fraction and the total velocity for particle
         
                // Get neutral fraction (only need to calculate if we haven't already 
                // done so)
                if(Temp_P[i].XHI < 0)
                {
                    Temp_P[i].XHI = get_neutral_frac(Temp_P[i]);
                }

                // Get total particle velocity in km/s (only need to calculate if we 
                // haven't already). This is only true if using LOS that are all oriented
                // the same way (all || to x-axis, for example). In the GENERAL_LOS case,
                // this value, which is actually the component of the total vel along
                // the LOS, does need to be recalculated.
                #ifndef GENERAL_LOS
                    if(Temp_P[i].tot_vel_flag == 0)
                    {
                        Temp_P[i].tot_vel = get_total_vel(Temp_P[i], t);
                    }
                #else
                    Temp_P[i].tot_vel = get_total_vel(Temp_P[i], t);
                #endif

                // Now get contribution of the particle to the pixels it overlaps
                get_particle_contributions(Temp_P[i], dist2);
            }

            // Check for errors
            if(exit_code_local == 0)
            {
                if(ntasks > 1)
                {
                    free(Temp_P);
                }
                return;
            }
        }
   
        // Free resources
        if(ntasks > 1)
        {
            free(Temp_P);
        }
    }

    // Check to make sure that all of the particles were transferred
    if(nleft != 0)
    {
        sprintf(error_message_local, "Error, nleft does not equal zero "\
            "(task %d)!\n", thistask);
        exit_code_local = 0;
        return;
    }

    // Use the results of B1, B2, and B3 to get tau for each pixel (eval B6.)
    pixel_optical_depth();
}



/********************************************
           get_perpendicular_dist
********************************************/
float get_perpendicular_dist(PARTICLE_DATA pd, float *param)
{
    // This function calculates the perpedicular distance between the particle and
    // the LOS. 
    // See: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    float t_numerator   = 0.0;
    float t_denominator = 0.0;
    float t             = 0.0;
    float dist2         = 0.0;
    int j               = 0;

    // Get t
    for(j = 0; j < 3; j++)
    {
        t_numerator +=  (pixels[0].Pos[j] - pd.Pos[j]) *
                        (pixels[npixels-1].Pos[j] - pixels[0].Pos[j]);

        t_denominator +=  (pixels[npixels-1].Pos[j] - pixels[0].Pos[j]) *
                          (pixels[npixels-1].Pos[j] - pixels[0].Pos[j]);
    }

    t = -1.0 * (t_numerator / t_denominator);

    // Use t to get perpendicular distance between particle and LOS
    // Both the particle and pixel positions are comoving GUL.
    for(j = 0; j < 3; j++)
    {
        dist2 +=    ((pixels[0].Pos[j] - pd.Pos[j]) + ((pixels[npixels-1].Pos[j] - \
                    pixels[0].Pos[j]) * t)) * ((pixels[0].Pos[j] - pd.Pos[j]) + \
                    ((pixels[npixels-1].Pos[j] - pixels[0].Pos[j]) * t));
    }

    // "Return" t so I can reuse it in velocity calc without having to refind it there
    *param = t;

    // Check for nans and infs
    if(isnan(t) == 1)
    {
        sprintf(error_message_local, "Error, t is a nan (task %d, part %d)!\n", \
            thistask, pd.id);
        exit_code_local = 0;
        return t;
    }

    if(isfinite(t) == 0)
    {
        printf(error_message_local, "Error, t is infinite (task %d, part %d)!\n", \
            thistask, pd.id);
        exit_code_local = 0;
        return t;
    }

    if(isnan(dist2) == 1)
    {
        printf(error_message_local, "Error, dist2 is a nan (task %d, part %d)!\n", \
            thistask, pd.id);
        exit_code_local = 0;
        return dist2;
    }

    if(isfinite(dist2) == 0)
    {
        printf(error_message_local, "Error, dist2 is infinite (task %d, part %d)!\n", \
            thistask, pd.id);
        exit_code_local = 0;
        return dist2;
    }

    return dist2;
}



/********************************************
               get_total_vel
*********************************************/
float get_total_vel(PARTICLE_DATA pd, float t)
{
    // This function calculates the total velocity of the particle, v = v_hub + v_pec. 
    // Returns v in km/s. See Amendola eq. 2.31: The total speed of an object along the
    // direction from the observer (located at the origin) to the object is given by:
    // v = \vec{\dot{r}} \cdot \frac{\vec{r}}{r} = Hr + \vec{v}_p \cdot \frac{\vec{r}}{r}
    // where \frac{\vec{r}}{r} is just \hat{l}, which is the unit vector that points along
    // the LOS. Here r in the Hubble velocity term, I think, should be the distance
    // between los_end and the point where the perpendicular distance between the
    // particle and the LOS intersects the LOS. That way, the Hubble velocity is wrt to
    // the observer. I know my LOS are parallel to the x-axis, but I'm going to keep 
    // this general.

    float v_pec_los = 0.0;
    float v_hub     = 0.0;
    float x_inter   = 0.0;
    float y_inter   = 0.0;
    float z_inter   = 0.0;
    float r         = 0.0;

    // Get the component of the particle's peculiar velocity along the LOS. Because the
    // magnitude of l hat is 1, all that matters for this dot product is the direction,
    // not the actual origin of r. 

    // NOTE: The gadget snapshot contains u. v_p = sqrt(a) * u. Also, they are written in
    // GUVs, which, if the default system is chosen, makes them km/s. I always use the 
    // default. So P.Vel has units of km/s. dx, dy, dz, and los_len have units of 
    // comoving GUL. The comoving GUL cancels to leave km/s. 
    v_pec_los = ((pd.Vel[0] * dx) + (pd.Vel[1] * dy) + (pd.Vel[2] * dz)) \
            * sqrt(header.time) / los_len;

    // Get point where perp dist intersects LOS. See the page for the perpendicular 
    // distance calc, but basically what's going on is in the perpendicular dist func, 
    // I found the t that corresponds to the intersection point on the LOS, and I've 
    // passed it here so I don't need to refind it. Can just plug it into the 3d equation 
    // for the LOS

    // The units here are comoving GUL. t is dimensionless
    x_inter = pixels[0].Pos[0] + ((pixels[npixels - 1].Pos[0] - pixels[0].Pos[0]) * t);   
    y_inter = pixels[0].Pos[1] + ((pixels[npixels - 1].Pos[1] - pixels[0].Pos[1]) * t);   
    z_inter = pixels[0].Pos[2] + ((pixels[npixels - 1].Pos[2] - pixels[0].Pos[2]) * t);

    // Get dist from los_end to intersection point. This will have units of comoving 
    // GUL^2
    r = ((los_end[0] - x_inter) * (los_end[0] - x_inter)) + \
        ((los_end[1] - y_inter) * (los_end[1] - y_inter)) + \
        ((los_end[2] - z_inter) * (los_end[2] - z_inter));

    // Units of comoving GUL
    r = sqrt(r);

    // Get Hubble velocity along LOS. Units of comoving GUL/s
    v_hub = H_sim * r;

    // Convert from comoving GUL/s to physical GUL/s
    v_hub = v_hub * header.time;

    // Convert from physical GUL/s to physical km/s
    v_hub = v_hub * GUL_in_cm / (header.HubbleParam * 100000.0);

    // Flag the particle as being done
    pd.tot_vel_flag = 1;

    // Check for nans and infs
    if(isnan(v_hub) == 1)
    {
        sprintf(error_message_local, "Error, v_hub is a nan (task %d, part %d)!\n", \
            thistask, pd.id);
        exit_code_local = 0;
        return v_hub;
    }

    if(isfinite(v_hub) == 0)
    {
        sprintf(error_message_local, "Error, v_hub is infinite "\
            "(task %d, part %d)!\n", thistask, pd.id);
        exit_code_local = 0;
        return v_hub;
    }

    if(isnan(v_pec_los) == 1)
    {
        sprintf(error_message_local, "Error, v_pec_los is a nan (task %d, part %d)!\n", \
            thistask, pd.id);
        exit_code_local = 0;
        return v_pec_los;
    }

    if(isfinite(v_pec_los) == 0)
    {
        printf(error_message_local, "Error, v_pec_los is infinite "\
            "(task %d, part %d)!\n", thistask, pd.id);
        exit_code_local = 0;
        return v_pec_los;
    }

    return (v_hub + v_pec_los);
}



/********************************************
        get_particle_contributions
********************************************/
void get_particle_contributions(PARTICLE_DATA pd, float perp_dist2)
{
    // This function evaluates eqs. B1, B2, and B3 from Bertone05. We don't need to 
    // loop over every pixel in the LOS, just those that could possibly be overlapped
    //  by the particle's smoothing sphere. This max range of pixels is achieved by 
    // pretending the particle is directly on the LOS and finding how many pixels the
    // smoothing length spans. 

    int i    = 0;
    int n    = 0;   // Index of pixel that particle projects into
    int hn   = 0;   // Number of pixels spanned by smoothing length
    float s  = 0.0; // The line segment connecting particle and los_start
    float b  = 0.0; // The line segment between pixels[0] and where d
                    // intersects the LOS
    float w  = 0.0; // Smoothing kernel as defined in Gadget2 paper (eq. 4)
   
    // Figue out which pixel the particle projects down into. To do this, I set up a
    // right triangle formed by the points A (los_start), P (particle location), and
    // inter (the point where the perpendicular distance line segment between P and
    // the LOS intersects the LOS). This gives three lengths: s = A->P, d = P->inter,
    // and b = A-> inter.

    // d is the perpendicular distance between the particle and the LOS and is known.
    // Since the locations of both A and P are known, this allows s to be found via
    // the distance equation. Then, since s, d, and b form a right triangle, Pythagoras
    // gives me b.

    // b = n * delta_pix, where n is the number of pixels between the start of the LOS 
    // and the pixel that the particle projects down into. With b and delta_pix known, 
    // this gives me n, which is the index of the pixel the particle projects down into

    // Units of comoving GUL^2
    s = ((pd.Pos[0] - los_start[0]) * (pd.Pos[0] - los_start[0])) + \
        ((pd.Pos[1] - los_start[1]) * (pd.Pos[1] - los_start[1])) + \
        ((pd.Pos[2] - los_start[2]) * (pd.Pos[2] - los_start[2]));

    // Units of comoving GUL
    b = sqrt(s - perp_dist2);

    // Both b and delta_pix have units of comoving GUL
    n = (int)(b / delta_pix);

    // Now get the size of the smoothing length in pixel units (that is, h = x # of 
    // pixels) This is the same thing that was used to get n, but now we round up 
    // for safety's sake

    // Both h and delta_pix have units of comoving GUL 
    hn = (int)ceil(pd.Hsml / delta_pix);

    // Loop over possibly overlapped pixels
    for(i = (n - hn); i <= (n + hn); i++)
    {
        // Make sure we're in a valid pixel range first
        if((i >= 0) && (i < npixels))
        {
            // Get smoothing kernel in comoving GUL^-3
            w = get_smoothing_kernel(pd, i);

            // Get Bertone05 eq. B4, but without the 1/h^3 that she has as that factor 
            // is already included in the definition of w.

            // Units of comoving GUM / GUL^3. The units on w don't actually matter, 
            // because they get normalized out when actually calculating tau in
            // pixel_optical_depth() 
            w = pd.Mass * hydrogen_mass_frac * pd.XHI * w;

            // Multiply by baryon fraction if we're using dm only
            if(header.npartTot[0] == 0)
            {
                w = baryon_fraction * w;
            }

            // Eq. B1
            // Units of comoving GUM / GUL^3
            pixels[i].Rho_HI  += w;
            pixels[i].RhoTot  += w / (hydrogen_mass_frac * pd.XHI);
         
            // Eq. B2
            // Units of comoving GUM * K / GUL^3
            pixels[i].Temp_HI += w * pd.temperature;

            // Eq. B3
            // Units of comoving GUM / GUL^3 * (km/s)
            pixels[i].Vel_HI += w * pd.tot_vel; 

            // In order to get the Doppler redshift for each pixel according to
            // eq. 11 in the trident paper, I need v_i, which is the magnitude of
            // the velocity of the gas in the pixel. I know I've already done this calc
            // in get_total_vel, but with 1024^3 particles, saving even one extra float
            // in the P structure equates to 4 extra GB of RAM, which is a lot, so I
            // figured it was better to just redo the calculation of v_pec_los and v_pec
            // here rather than use that much extra memory.
            #ifdef LOS_REDSHIFT
                // This is the magnitude of the peculiar velocity for the HI gas (since
                // I'm using w, which has XHI built in) weighted by the HI density.
                // Units of comoving GUM / GUL^3 * (km/s). The density is normalized out
                // later, but this kind of bothers me. The standard sph way of getting
                // a quantity is A = \sum m * A * W / rho, and I'm not dividing by rho
                // in any of the quantities in this function. But those are the equations
                // given in Bertone05 and Theuns98. In Bertone's code she does the same
                // thing I am and divides rho out later using the full pixel quantity
                pixels[i].v_pec += w * sqrt(header.time * ((pd.Vel[0]*pd.Vel[0]) + \
                    (pd.Vel[1] * pd.Vel[1]) + (pd.Vel[2] * pd.Vel[2])));

                pixels[i].v_pec_los += w * ((pd.Vel[0] * dx) + (pd.Vel[1] * dy) + \
                    (pd.Vel[2] * dz)) * sqrt(header.time) / los_len;

                // Check for nans and infs
                if(isnan(pixels[i].v_pec) == 1)
                {
                    sprintf(error_message_local, "Error, v_pec is a nan "\
                        "(task %d, part %d)!\n", thistask, pd.id);
                    exit_code_local = 0;
                    return;
                }

                if(isfinite(pixels[i].v_pec) == 0)
                {
                    sprintf(error_message_local, "Error, v_pec is infinite "\
                        "(task %d, part %d)!\n", thistask, pd.id);
                    exit_code_local = 0;
                    return;
                }

                if(isnan(pixels[i].v_pec_los) == 1)
                {
                    printf(error_message_local, "Error, pixels.v_pec_los is a nan "\
                        "(task %d, part %d)!\n", thistask, pd.id);
                    exit_code_local = 0;
                    return;
                }

                if(isfinite(pixels[i].v_pec_los) == 0)
                {
                    sprintf(error_message_local, "Error, pixels.v_pec_los is infinite "\
                        "(task %d, part %d)!\n", thistask, pd.id);
                    exit_code_local = 0;
                    return;
                }
            #endif

            // Check for nans and infs
            if(isnan(w) == 1)
            {
                sprintf(error_message_local, "Error, w is a nan "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }

            if(isfinite(w) == 0)
            {
                printf(error_message_local, "Error, w is infinite "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }

            if(isnan(pixels[i].Rho_HI) == 1)
            {
                sprintf(error_message_local, "Error, rho_HI is a nan "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }

            if(isfinite(pixels[i].Rho_HI) == 0)
            {
                sprintf(error_message_local, "Error, rho_HI is infinite "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }

            if(isnan(pixels[i].RhoTot) == 1)
            {
                sprintf(error_message_local, "Error, rhoTot is a nan "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }

            if(isfinite(pixels[i].RhoTot) == 0)
            {
                sprintf(error_message_local, "Error, rhoTot is infinite "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }

            if(isnan(pixels[i].Temp_HI) == 1)
            {
                sprintf(error_message_local, "Error, Temp_HI is a nan "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }

            if(isfinite(pixels[i].Temp_HI) == 0)
            {
                sprintf(error_message_local, "Error, Temp_HI is infinite "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }

            if(isnan(pixels[i].Vel_HI) == 1)
            {
                sprintf(error_message_local, "Error, Vel_HI is a nan "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }

            if(isfinite(pixels[i].Vel_HI) == 0)
            {
                sprintf(error_message_local, "Error, Vel_HI is infinite "\
                    "(task %d, part %d)!\n", thistask, pd.id);
                exit_code_local = 0;
                return;
            }
        }
    }
}



/********************************************
           get_smoothing_length
********************************************/
float get_smoothing_kernel(PARTICLE_DATA pd, int i)
{
    // Evaluates eq. 4 in the Gadget2 paper

    float d = 0.0;
    float q = 0.0;
    float w = 0.0;

    // Get distance between pixel and particle
    // Units of comoving GUL^2
    d = ((pd.Pos[0] - pixels[i].Pos[0]) * (pd.Pos[0] - pixels[i].Pos[0])) + \
        ((pd.Pos[1] - pixels[i].Pos[1]) * (pd.Pos[1] - pixels[i].Pos[1])) + \
        ((pd.Pos[2] - pixels[i].Pos[2]) * (pd.Pos[2] - pixels[i].Pos[2]));

    // Units of comoving GUL
    d = sqrt(d);

    // Get ratio between particle distance and smoothing length
    // Both d and h have units of comoving GUL, so q is dimensionless
    q = d / pd.Hsml;

    // Ensure this is > 0 (it should be)
    if(q < 0.0)
    {
        sprintf(error_message_local, "Error, q < 0 "\
            "(task %d, part %d)\n", thistask, pd.id);
        exit_code_local = 0;
        return w;
    }

    // Get w
    if(q <= 0.5)
    {
        w = 1.0 - (6.0 * q * q) + (6.0 * q * q * q);
    }

    else if(q <= 1.0)
    {
        w = 2.0 * pow(1.0 - q, 3.0);
    }

    else if(q > 1.0)
    {
        w = 0.0;
    }

    // Units of comoving GUL^-3
    w = w * 8.0 / (M_PI * pd.Hsml * pd.Hsml * pd.Hsml);

    return w;
}




/********************************************
            pixel_optical_depth
********************************************/
void pixel_optical_depth(void)
{
    // Evaluates Bertone05 eq. B6

    int i          = 0;
    int j          = 0;
    double d_pix   = 0.0;
    double doppler = 0.0;
    double tauc    = 0.0;
    double tau     = 0.0;

    #ifdef LOS_REDSHIFT
        double z_eff_plus_one = 0.0;
        double dist           = 0.0;
    #endif

    // Normalize the HI density weighted temperature and velocity for each pixel
    // by dividing out the HI density. This is why the units of pixel.Rho_HI don't
    // matter in get_particle_contributions().
    for(i = 0; i < npixels; i++)
    {
        // Make sure the pixel has gas in it
        if(pixels[i].Rho_HI > 0.0)
        {
            // Temperature in Kelvin
            // Temp_HI has units of comoving GUM * K / GUL^3, and
            // Rho_HI has units of comoving GUM / GUL^3.
            pixels[i].Temp_HI = pixels[i].Temp_HI / pixels[i].Rho_HI;

            // Velocity in km/s proper
            // Vel_HI has units of comoving GUM / GUL^3 * (km/s)
            // Rho_HI has units of comoving GUM / GUL^3.
            pixels[i].Vel_HI  = pixels[i].Vel_HI / pixels[i].Rho_HI;

            // Get the Doppler redshift for this pixel, if desired. This is eq. 11
            // from the trident paper. Vel_HI is already along the LOS
            // (see get total_vel()) so I don't need to worry about the cos(\theta).
            // Also get the effective redshift and, while I'm at it, get the lambda and
            // z_vel quantities for the pixel, as well, since z_doppler was the last
            // missing ingredient. z_vel is used when getting tau, so this is actually
            // needed here.
            #ifdef LOS_REDSHIFT
                // Divide out rho_hi to give the velocities units of km/s
                pixels[i].v_pec /= pixels[i].Rho_HI;
                pixels[i].v_pec_los /= pixels[i].Rho_HI;

                // Evaluate eq. 11 from the trident paper. This is actually 1 + z_dopp
                pixels[i].z_doppler = (1.0 + (pixels[i].v_pec_los / speed_of_light)) / \
                    sqrt(1.0 - pow(pixels[i].v_pec / speed_of_light, 2.0));

                // Get the effective redshift (trident paper eq. 12)
                z_eff_plus_one = pixels[i].z_doppler * (1.0 + pixels[i].z_cos);

                // Now we shift the lambdas by this amount
                pixels[i].lambda = LYMANA * z_eff_plus_one;
            #endif

            // Total density in proper cgs (proper g/cm^3)
            // Right now it's in comoving GUM / GUL^3
            // Convert to physical
            pixels[i].RhoTot = pixels[i].RhoTot / (header.time * header.time * \
                header.time);

            // Convert to proper cgs from proper GUM / GUL^3
            pixels[i].RhoTot = pixels[i].RhoTot * header.HubbleParam * \
                header.HubbleParam;
            pixels[i].RhoTot = pixels[i].RhoTot * GUD_in_cgs;

            // Get the HI number density in proper cgs (proper cm^-3)
            // Convert Rho_HI from comoving GUM / GUL^3 to proper GUM / GUL^3
            pixels[i].nHI = pixels[i].Rho_HI / (header.time * header.time * header.time);

            // Convert from proper GUM / GUL^3 to proper g / cm^3
            pixels[i].nHI = pixels[i].nHI * header.HubbleParam * header.HubbleParam;
            pixels[i].nHI = pixels[i].nHI * GUD_in_cgs;

            // Convert from proper g / cm^3 to cm^-3 by dividing by mass of hydrogen in g
            pixels[i].nHI = pixels[i].nHI / hydrogen_mass;
        }

        #ifdef LOS_REDSHIFT
            // Get the Hubble velocity for each pixel (z_vel). This is just H(z_i) * 
            // pixel dist from observer. This is in comoving GUL. This needs to happen for
            // the pixel regardless of whether or not there is gas in it, which is why
            // this is outside the above if statement on Rho_HI.
            dist = sqrt(pow(pixels[i].Pos[0] - los_end[0], 2.0) + \
                pow(pixels[i].Pos[1] - los_end[1], 2.0) + \
                pow(pixels[i].Pos[2] - los_end[2], 2.0));

            // Convert to physical GUL
            dist *= 1.0 / (1.0 + pixels[i].z_cos);

            // Convert to km
            dist *= GUL_in_cm / (header.HubbleParam * 100000.0);

            // Multiply by H(z_i) to get z_vel (really v_H for the pixel) in km/s
            pixels[i].z_vel = get_hubble(pixels[i].z_cos) * dist;
        #endif
    }

    // Evaluate tau
    for(i = 0; i < npixels; i++)
    {
        // Get current pixel's doppler parameter (Bertone05 eq. B7) in proper cm/s
        // Need Boltzmann in erg / K, and protonmass in g. 
        // Gives doppler in cm/s
        doppler = sqrt(2.0 * Boltzmann * pixels[i].Temp_HI / protonmass);

        // Convert to km/s
        doppler = doppler / 100000.0;

        // If there's gas in pixel, get the non-exponential part of B7
        if(doppler > 0.0)
        {
            // Doppler is in proper km/s, so we need speed of light in proper km/s
            // sigma_alpha is in proper cm^2, nHI is in proper cm^-3, so I need delta_pix
            // in proper cm to make the whole thing dimensionless.
         
            // Convert delta_pix from comoving GUL to proper GUL
            d_pix = delta_pix * header.time;

            // Convert proper GUL to proper cm
            d_pix = d_pix * GUL_in_cm / header.HubbleParam;
         
            tauc = sigma_alpha * speed_of_light * pixels[i].nHI * d_pix;
            tauc = tauc / (sqrt(M_PI) * doppler);
        }

        else
        {
            tauc = 0.0;
        }

        // Get this pixel's contribution to all other pixel's optical depth
        for(j = 0; j < npixels; j++)
        {
            // Only contributes if there's gas in the current pixel
            if(doppler > 0.0)
            {
                // z_vel is in km/s, so tau is dimensionless
                tau = -1.0 * ((pixels[j].z_vel - pixels[i].Vel_HI) * \
                      (pixels[j].z_vel - pixels[i].Vel_HI)) / (doppler * doppler);

                tau = exp(tau);

                tau = tauc * tau; 
            }

            else
            {
                tau = 0.0;
            }

            // Add to the total optical depth of the pixel
            pixels[j].tau += tau;

            if(isnan(pixels[j].tau) == 1)
            {
                sprintf(error_message_local, "Error, pixel %d is a nan (task %d)\n",\
                     j, thistask);
                exit_code_local = 0;
                return;
            }
        }
    }

    // Check for nans and infs in tau
    for(j = 0; j < npixels; j++)
    {
        if(isnan(pixels[j].tau) == 1)
        {
            sprintf(error_message_local, "Error, pixel %d tau is a nan (task %d)\n", \
                j, thistask);
            exit_code_local = 0;
            return;
        }

        if(isfinite(pixels[j].tau) == 0)
        {
            sprintf(error_message_local, "Error, pixel %d tau is inf (task %d)!\n", \
                j, thistask);
            exit_code_local = 0;
            return;
        }
    }
}
