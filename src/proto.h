/*****************************************************************************************
Title:   proto.h
Purpose: Contains all function prototypes
Notes:
*****************************************************************************************/
#ifndef ALLVARS_H
   #include "allvars.h"
#endif



/********************************************
                   de.c
********************************************/
void   de_check(void);
double expansion_factor(double);
double get_dark_factor(double);
double get_hubble(double);



/********************************************
                  endrun.c
********************************************/
void check_exit_code(void);
void free_globals(void);



/********************************************
                  gamma.c
********************************************/
double **get_gamma_integrand(double **, double *, int, int);
double **get_HM_J(long *, int, int);
double *get_HM_lambdas(long *, int *, int *);
double *get_HM_redshifts(long *, int *);
double *integrate_gamma(double **, double *, int, int);
double interp_gamma(double *, double *, int);
void   photoionization_rate(void);



/********************************************
                   init.c
********************************************/
void init(int, char **);
int  checkforerror(char *, int, char *);
void make_custom_mpi_type(void);



/********************************************
                log_write.c
********************************************/
void log_write_compile_options(void);
void log_write_errors(void);
void log_write_parameters(void);
void log_write_timings(void);



/********************************************
                   los.c
********************************************/
void los(void);



/********************************************
              los_redshift.c
********************************************/
double integrand(double, void *);
float  redshift_range(float);
double redshift_roots(double, void *);



/********************************************
                  neutral.c
********************************************/
float get_neutral_frac(PARTICLE_DATA);
float get_electron_number_density(PARTICLE_DATA);
float get_recombination_rate(PARTICLE_DATA);
float get_collisional_rate(PARTICLE_DATA);



/********************************************
                  normal.c
********************************************/
float get_tausim(void);
void  normalize_spectrum(void);



/********************************************
              optical_depth.c
********************************************/
void  optical_depth(void);
float get_perpendicular_dist(PARTICLE_DATA, float *);
float get_total_vel(PARTICLE_DATA, float);
void  get_particle_contributions(PARTICLE_DATA, float);
float get_smoothing_kernel(PARTICLE_DATA, int);
void  pixel_optical_depth(void);



/********************************************
               read_snapshot
********************************************/
void   load_snapshot(void);
int    copy_to_master(PARTICLE_DATA *, SNAP_HEADER, int);
size_t my_fread(void *, size_t, size_t, FILE *);
int    block_check(enum fields, int, int, SNAP_HEADER);
int    get_block_size(enum fields, SNAP_HEADER);
int    get_with_masses(SNAP_HEADER);



/********************************************
                  write.c
********************************************/
void write_spectrum(int);
