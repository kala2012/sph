#ifndef _SYSTEM_H_
#define _SYSTEM_H_
#include <malloc.h>



typedef struct __system {
  int dim; // Dimension of the system
  int *itype;
  int nx;
  int ny;
  int NumberOfInteractingParticles;
  int MaxNumberOfParticles;
  int ntotal; /* total number of fluid particles */
  int nvirt;  /* number of virtual particles */
  int nvll;   /* number of boundary particles */
  double **Position; /* initialize the position vector */
  double **Velocity; /* initialize the velocity vector */
  double **Velocity_xsph; /* move particles with this average velocity */
  double *mass; /* initialize particle mass */
  double *Pressure; /* initialize particle pressure */
  double *Energy; /* initialize internal energy */
  double *hsml; /* initialize smoothing length h */
  double *rho; /* initialize particle density */
  double **Velocity_min;
  double *drhodt;
  double **dvdt;
  double **av;
  double *rho_min;
  double rho0;
  // temporary variables used for the calculations
  double alpha_d;
  double alpha;
  double SpeedOfSound;
  int *i_pair, *j_pair;
  int *i_inflow;
  double *rij2, *w;
  double **dwdx;
  int inflowParticles;
  double CompressionFactor;
  double AdiabaticConstant;
} System;

System *CreateSystem(const int dim, 
		     const int nx, 
		     const int ny, 
		     const int MaxNumberOfParticles, 
		     const double alpha, 
		     const double SpeedOfSound, 
		     const double rho0,
		     const double AdiabaticConstant);
void DestroySystem(System *a);
#endif
