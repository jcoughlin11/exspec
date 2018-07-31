/*****************************************************************************************
Title:   allvars.c
Purpose: Globals!
Notes:
*****************************************************************************************/
#include <gsl/gsl_rng.h>
#include <mpi.h>
#include "allvars.h"



// MPI variables
int ntasks;
int thistask;
int nlos_local;
int buffer_size = 100000000; // Size of buffer in bytes. Gadget uses 100 MB, 
                             // so I will too (ignores the 1024 stuff)
MPI_Datatype mpi_header_type;
MPI_Datatype mpi_particle_type;



// Constants
const double LYMANA         = 1.2156e-5;  // In cm
const double PLANCK         = 6.626e-27;  // In erg s
const double Boltzmann      = 1.38e-16;   // In erg / K
const double hydrogen_mass  = 1.6737e-24; // In g
const double protonmass     = 1.6726e-24; // In g
const double sigma_alpha    = 1.12e-18;   // In cm^2
const double speed_of_light = 2.99e5;     // In km/s



char   snap_file[256];
char   hm_file[256];
char   parent_dir[256];
char   spec_dir[256];
char   *error_message_local = NULL;
int    err_msg_size         = 256;
int    npixels              = 0;
int    nlos                 = 0;
int    numpart              = 0;
int    dark_param           = 0;
int    max_iters            = 0;
int    start_index          = 0;
int    end_index            = 0;
int    exit_code_local      = 1;
int    exit_code_global     = 1;
int    write_error_to_log   = 1;
double baryon_fraction      = 0.0;
double delta_pix            = 0.0;
double GUM_in_g             = 0.0;
double GUL_in_cm            = 0.0;
double GUD_in_cgs           = 0.0;
double w0                   = 0.0;
double wa                   = 0.0;
double H_sim                = 0.0;
double GammaXHI             = 0.0;
double dx                   = 0.0;
double dy                   = 0.0;
double dz                   = 0.0;
double los_len              = 0.0;
double los_end[3]           = {0.0, 0.0, 0.0};
double los_start[3]         = {0.0, 0.0, 0.0};
double hydrogen_mass_frac   = 0.0;
double helium_mass_frac     = 0.0;
double tolerance            = 0.0;
double molweight            = 0.0;
double min_pixels_frac      = 1.0 / 3.0;
FILE *log_file              = NULL;
gsl_rng *rng                = NULL;
SNAP_HEADER header;
PARTICLE_DATA *P            = NULL;
PIXEL *pixels               = NULL;
TimingContainer timer;
