/*****************************************************************************************
Title:   allvars.h
Purpose: Structure and global variables
Notes:
*****************************************************************************************/
#ifndef ALLVARS_H
    #define ALLVARS_H

    #include <gsl/gsl_rng.h>
    #include <gsl/gsl_integration.h>
    #include <mpi.h>

    /********************
        mpi variables
    ********************/
    extern int ntasks;
    extern int thistask;
    extern int nlos_local;
    extern int buffer_size;
    extern MPI_Datatype mpi_header_type;
    extern MPI_Datatype mpi_particle_type;



    /********************
          Constants
    ********************/
    extern const double LYMANA;
    extern const double PLANCK;
    extern const double Boltzmann;
    extern const double hydrogen_mass;
    extern const double protonmass;
    extern const double sigma_alpha;
    extern const double speed_of_light;



    /********************
           Globals
    ********************/
    extern char     snap_file[256];
    extern char     hm_file[256];
    extern char     parent_dir[256];
    extern char     spec_dir[256];
    extern char     *error_message_local;
    extern int      err_msg_size;
    extern int      npixels;
    extern int      nlos;
    extern int      numpart;
    extern int      dark_param;
    extern int      max_iters;
    extern int      start_index;
    extern int      end_index;
    extern int      exit_code_local;
    extern int      exit_code_global;
    extern int      write_error_to_log;
    extern double   baryon_fraction;
    extern double   delta_pix;
    extern double   GUM_in_g;
    extern double   GUL_in_cm;
    extern double   GUD_in_cgs;
    extern double   w0;
    extern double   wa;
    extern double   H_sim;
    extern double   GammaXHI;
    extern double   dx;
    extern double   dy;
    extern double   dz;
    extern double   los_len;
    extern double   los_end[3];
    extern double   los_start[3];
    extern double   hydrogen_mass_frac;
    extern double   helium_mass_frac;
    extern double   tolerance;
    extern double   molweight;
    extern double   min_pixels_frac;
    extern gsl_rng  *rng;
    extern FILE     *log_file;



    /********************
         Structures
    ********************/
    // Header
    typedef struct SNAP_HEADER
    {
        int    npart[6];
        double mass[6];
        double time;
        double redshift;
        int    flag_sfr;
        int    flag_feedback;
        int    npartTot[6];
        int    flag_cooling;
        int    num_files;
        double BoxSize;
        double Omega0;
        double OmegaLambda;
        double HubbleParam;
        char   fill[96];
    } SNAP_HEADER;



    // Particle
    typedef struct PARTICLE_DATA
    {
        float Pos[3];
        float Vel[3];
        float Hsml;
        float Rho;
        float temperature;
        float Mass;
        float XHI;
        float tot_vel;
        int   id;
        int   type;
        int   tot_vel_flag;
    } PARTICLE_DATA;



    // Pixel
    typedef struct PIXEL
    {
        double Pos[3];
        double Rho_HI;
        double RhoTot;
        double Temp_HI;
        double Vel_HI;
        double nHI;
        double z_vel;
        double lambda;
        double tau;

        #ifdef LOS_REDSHIFT
            double z_cos;
            double z_doppler;
            double v_pec;
            double v_pec_los;
        #endif
    } PIXEL;

    // Timing Container
    typedef struct TimingContainer
    {
        double init_time;
        double read_snapshot_time;
        double gamma_time;
        double spectra_time;
        double run_time;
        int    is_instantiated;
    } TimingContainer;



    // Block fields
    enum fields
    {
        HEADER,
        POS,
        VEL,
        IDS,
        MASS,
        TEMP,
        RHO,
        HSML
    };



    // Integration parameters
    typedef struct Integration_Params
    {
        double upper_lim;
        double answer;
        double hubble_dist;
        gsl_integration_workspace *w;
    } Integration_Params;



    extern SNAP_HEADER header;
    extern PARTICLE_DATA *P;
    extern PIXEL *pixels;
    extern TimingContainer timer;
#endif
