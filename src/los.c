/*****************************************************************************************
Title: los.c
Purpose: Generates an los and sets up the pixels along it.
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"



/********************************************
                    los
********************************************/
void los(void)
{
    // This function picks an los and sets up the
    // pixels along it. All los are parallel to the
    // x axis of the simualation volume. 

    int   i = 0;

    #ifdef GENERAL_LOS
        int cube_face_start = 0;
    #endif

    #ifdef LOS_REDSHIFT
        float *dredshift = NULL;
        double z_end     = 0.0;
        int    j         = 0;
    #else
        float dist = 0.0;
    #endif

    // Get los start point in COMOVING GUL
    #ifndef GENERAL_LOS
        los_start[0] = 0.0; 
        los_start[1] = gsl_rng_uniform(rng) * header.BoxSize;
        los_start[2] = gsl_rng_uniform(rng) * header.BoxSize;

        // Get los end point in COMOVING GUL
        los_end[0] = header.BoxSize;
        los_end[1] = los_start[1];
        los_end[2] = los_start[2];

        // Choose center of box if we're debugging so we have the
        // same LOS every time
        #ifdef DEBUGGING
            los_start[1] = header.BoxSize / 2.0;
            los_start[2] = header.BoxSize / 2.0;
            los_end[1]   = los_start[1];
            los_end[2]   = los_start[2];
        #endif

        #ifdef MPI_DEBUGGING
            los_start[1] = header.BoxSize / 2.0;
            los_start[2] = header.BoxSize / 2.0;
            los_end[1]   = los_start[1];
            los_end[2]   = los_start[2];
        #endif
    #else
        // With this option (GENERAL_LOS), the LOS does not have to be parallel to the
        // x-axis. Here we randomly choose a starting face of the cube. Then we randomly
        // choose a starting point on that face. Then we choose a random position on the
        // opposite face (for now, and for simplicity).

        // Recall that the origin is the back left corner of the cube ---> x, up is y, and
        // out of the screen is z. I really only need to choose either the left, bottom,
        // or back faces, since, because I'm always going to the opposite face, going
        // right to left is functionally the same as left to right.
        // Faces: 0: x = 0 y-z face (the left face)
        // 1: y = 0, x-z face (the bottom)
        // 2: z = 0, x-y face (the back face)

        // Get the starting cube face. gsl_rng_uniform_int(0,n) produces integers in the 
        // range [0, n-1], so I need to have n = 3 to get values from 0 to 2.
        cube_face_start = gsl_rng_uniform_int(rng, 3);

        // Now get the opposite face. There's probably a better way of doing this, but,
        // whatever. Also get the starting and ending point. Positions in COMOVING GUL
        switch(cube_face_start)
        {
            // Start = left, end = right
            case 0:
                // Starting point on left face
                los_start[0] = 0.0;
                los_start[1] = gsl_rng_uniform(rng) * header.BoxSize;
                los_start[2] = gsl_rng_uniform(rng) * header.BoxSize;

                // End point on opposite face (right)
                los_end[0] = header.BoxSize;
                los_end[1] = gsl_rng_uniform(rng) * header.BoxSize;
                los_end[2] = gsl_rng_uniform(rng) * header.BoxSize;
                break; 

            // Start = bottom, end = top
            case 1:
                // Starting point on bottom face
                los_start[0] = gsl_rng_uniform(rng) * header.BoxSize;
                los_start[1] = 0.0;
                los_start[2] = gsl_rng_uniform(rng) * header.BoxSize;

                // End point on opposite face (top)
                los_end[0] = gsl_rng_uniform(rng) * header.BoxSize;
                los_end[1] = header.BoxSize;
                los_end[2] = gsl_rng_uniform(rng) * header.BoxSize;
                break; 

            // Start = back, end = front
            case 2:
                // Starting point on back face
                los_start[0] = gsl_rng_uniform(rng) * header.BoxSize;
                los_start[1] = gsl_rng_uniform(rng) * header.BoxSize;
                los_start[2] = 0.0;

                // End point on opposite face (front)
                los_end[0] = gsl_rng_uniform(rng) * header.BoxSize;
                los_end[1] = gsl_rng_uniform(rng) * header.BoxSize;
                los_end[2] = header.BoxSize;
                break; 

            // Default
            default:
                sprintf(error_message_local, "Error, incorrect cube face chosen "\
                    "(task %d)!\n", thistask);
                exit_code_local = 0;
                return;
        }
    #endif

    // Get the length and projections of the LOS in comoving GUL
    // Use actual start and endpoints of the LOS instead of first and last pixels
    dx = los_end[0] - los_start[0];
    dy = los_end[1] - los_start[1];
    dz = los_end[2] - los_start[2];

    los_len = sqrt((dx * dx) + (dy * dy) + (dz * dz));

    // If using GENERAL_LOS, then the number of pixels is a variable from LOS to LOS,
    // so we figure it out here (the idea is to keep the resolution (pixel width) the
    // same from LOS to LOS, but since each LOS will have a different length, this
    // necessitates changing npixels). Also, this is where I have to allocate memory
    // for the pixels structure when using general los
    #ifdef GENERAL_LOS
        npixels = (int)floor(los_len / delta_pix);

        // This should never happen because the los crosses from one cube face to the 
        // opposite face, but I'm going to check for it anyway
        if(npixels == 0)
        {
            sprintf(error_message_local, "Error, npixels is 0 (task %d)!\n", thistask);
            exit_code_local = 0;
            return;
        }

        // Allocate memory for pixels
        if(!(pixels = calloc(npixels, sizeof(PIXEL))))
        {
            sprintf(error_message_local, "Error, could not allocate memory for pixels "\
                "(task %d)!\n", thistask);
            exit_code_local = 0;
            return;
        }
    #endif

    // Allocate memory for dz, if needed
    #ifdef LOS_REDSHIFT
        if(!(dredshift = calloc(npixels, sizeof(float))))
        {
            sprintf(error_message_local, "Error, could not allocate memory for dz "\
                "(task %d)!\n", thistask);
            exit_code_local = 0;
            return;
        }
    #endif

    // If using LOS_REDSHIFT, each pixel along the LOS is assigned it's own redshift z_i.
    // In order to do this, we need to know the redshift range spanned by los_start and
    // los_end. That's what this does. los_start corresponds to the redshift of the
    // snapshot.
    #ifdef LOS_REDSHIFT
        z_end = redshift_range(los_len);

        // Check for errors
        if(z_end < 0.0)
        {
            return;
        }
    #endif

    // Initialize the pixels along los
    for(i = 0; i < npixels; i++)
    {
        // Pixel centers in comoving GUL. 
        #ifndef GENERAL_LOS
            // Because los is || to x-axis, y and z don't change. The x coordinate 
            // is the number of half delta multiples it takes to get to that pixel 
            // number (but my pixels are indexed from 0, hence the + 1)
            pixels[i].Pos[0] = los_start[0] + (0.5 * delta_pix * ((2 * i) + 1));
            pixels[i].Pos[1] = los_start[1];
            pixels[i].Pos[2] = los_start[2];
        #else
            // Pixel centers for a general los. For a general LOS, we have \vec{l} that
            // goes from A to B, so \vec{l} = (Bx - Ax)\hat{i} + (By - Ay)\hat{j} + 
            // (Bz - Az)\hat{k} and l = |\vec{l}| = sqrt((Bx-Ax)^2 + (By-Ay)^2 + 
            // (Bz-Az)^2). Then \hat{l} = \vec{l} / l. 

            // The vector \vec{P} goes from the origin to the point P, which is the 
            // center of the i'th pixel. The length of the segment |AP| is the number 
            // of half deltas to P, so |AP| = (2i + 1) * delta / 2. 

            // In general, Px = Ax + |AP|_x, where |AP|_x = \vec{AP} \cdot \hat{i}.
            // \vec{AP} = |AP|\hat{l}. So, |AP|_x = (2i + 1) * delta * dx / (2l). We then
            // have something similar for the y and z components.
            pixels[i].Pos[0] = los_start[0] + (((2.0 * i) + 1.0) * delta_pix * dx / \
                (2.0 * los_len)); 

            pixels[i].Pos[1] = los_start[1] + (((2.0 * i) + 1.0) * delta_pix * dy / \
                (2.0 * los_len));

            pixels[i].Pos[2] = los_start[2] + (((2.0 * i) + 1.0) * delta_pix * dz / \
                (2.0 * los_len));
        #endif

        // Tau components
        pixels[i].Rho_HI     = 0.0;
        pixels[i].RhoTot     = 0.0;
        pixels[i].Temp_HI    = 0.0;
        pixels[i].Vel_HI     = 0.0;
        pixels[i].nHI        = 0.0;
        pixels[i].tau        = 0.0;
        pixels[i].z_vel      = 0.0;
        pixels[i].lambda     = 0.0;

        #ifdef LOS_REDSHIFT
            pixels[i].z_cos     = 0.0;
            pixels[i].z_doppler = 0.0;
            pixels[i].v_pec     = 0.0;
            pixels[i].v_pec_los = 0.0;

            // Evaluate equation 8 from the trident paper. In this case, since delta_pix
            // is the same for every pixel, this will be the same value for each pixel.
            // However, I'm going to leave it in array form in case I decide to change
            // that in the future for some reason. Additionally, l in the trident paper
            // is negative because of their limits of integration on eq. 6 
            // (see los_redshift.c). Since z_A > z_B, I need to multiply by -1 in order
            // to make this whole thing positive. Redshift is dimensionless, los_len has
            // units of comoving GUL and delta_pix has units of comoving GUL.
            dredshift[i] = -1.0 * (delta_pix / los_len) * (z_end - header.redshift);
        #endif
    }


    // The reason this is separated from the above for loop is because I need z_vel to
    // start at 0 at pixels[npixels - 1] (i.e. the observer), so it's easiest to start
    // at the end of the pixels array and loop dowards to 0. If this were included in
    // the above for loop then when you get past half of the pixels, you end up
    // overwriting the value of z_vel saved in that pixel with 0 due to the above
    // initialization. Rather than rework the formula, I just added a second loop
    // because it's conceptually easier.
    #ifndef LOS_REDSHIFT
        for(i = 0; i < npixels; i++)
        {
            // Get z_vel (velocity in redshift space) and lambda (the corresponding
            // wavelength, used for plotting). From Theuns (private communication):
            // "Suppose you only had gas in a single bin (pixel i of N, say), and that 
            // the peculiar velocity of that bin is zero. The gas in this bin will 
            // then produce absorption, with optical depth a Gaussian centred at the 
            // location of that bin (ignoring natural broadening). So if v=Hr is the 
            // velocity extent of your box (H is the Hubble constant, r the physical 
            // extent of the box at that redshift, z), the absorption line is centred 
            // at velocity (i/N) * v, or in wavelength space, the absorption is centred 
            // at wavelength 1215.67 * (1+z) * (1+i/N * v/c) AA for HI."

            // The reason that z_vel is zeroed out at pix[npix-1] is because pix[npix-1]
            // corresponds to the end (read: close to point B) of the LOS, which is where
            // the observer is. Since the spectrum must be made from the POV of the
            // observer, I need to have the Hubble velocity increasing away from them,
            // towards point A. It just worked out that my pixel indexing started at
            // point A, which is the quasar.

            // Get the Hubble velocity for each pixel (z_vel). This is just H(z) * 
            // pixel dist from observer. This is in comoving GUL
            dist = sqrt(pow(pixels[i].Pos[0] - los_end[0], 2.0) + \
                pow(pixels[i].Pos[1] - los_end[1], 2.0) + \
                pow(pixels[i].Pos[2] - los_end[2], 2.0));

            // Convert dist from comoving GUL to physical GUL
            dist = dist / (1.0 + header.redshift);

            // Get z_vel in physical GUL / s. I commented out the line below because
            // it gives pixel[n - 1] a zvel = 0, which is not right, since the observer
            // is at B, not pixel[n - 1]. Those two points differ by delta_pix / 2.
            //pixels[(npixels - 1) - i].z_vel = i * header.time * H_sim * delta_pix;
            pixels[(npixels - 1) - i].z_vel = H_sim * dist;

            // Convert from physical GUL/s to physical km/s
            pixels[(npixels - 1) - i].z_vel = pixels[(npixels - 1) - i].z_vel * \
                GUL_in_cm / (header.HubbleParam * 100000.0);

            // Get lambda in cm
            pixels[(npixels - 1) - i].lambda = LYMANA * (1.0 + header.redshift) * \
                              (1.0 + (pixels[(npixels - 1) - i].z_vel / speed_of_light));
        }   
    #endif

    // If using LOS_REDSHIFT, we need to get pixels[i].lambda in a slightly different
    // way than we did above. The general idea is the same, though. The 1 + z_vel/c term
    // is what trident calls 1 + z_doppler. They use the relativisitic doppler equation.
    // Additionally, instead of using 1 + header.z, we need to use 1 + z_i. However, z_i
    // depends on all of the dz that come before it, which is why we couldn't calculate
    // this quantity in the above loop (see the sum in trident paper eq. 9).
    #ifdef LOS_REDSHIFT
        for(i = 0; i < npixels; i++)
        {
            // Initialize z_cos
            pixels[i].z_cos = header.redshift;

            // This is a little annoying. los_start and los_end (points A and B), refer
            // to the points on the cube faces. However, my pixel positions are given by
            // the centers of these bins. As such, pixel[0].pos != A, but offset by 1/2
            // delta_pix along the \hat{l} direction. However, the redshift interval I
            // found is z_A - z_B (i.e. between the points on the faces of the cube, not
            // the first and last pixel centers). As such, for the sum in trident eq. 9,
            // I think I need to go up to i-1 and then add 1/2 dz in order to "land" at
            // the location of the pixel center in redshift space.
            for(j = 0; j < i; j++)
            {
                pixels[i].z_cos -= dredshift[j];
            }

            // Out here, j == i, so we can tack on the extra half dz. I did it this way
            // so as to avoid having to check an if statement at every iteration
            pixels[i].z_cos -= 0.5 * dredshift[j];
        }

        // Clean up
        free(dredshift);
    #endif
}
