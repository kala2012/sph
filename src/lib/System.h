#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <malloc.h>
#include <stdlib.h>
#include "defs.h"
#include "dists.h"
#include "rbc.h"
#include "utils.h"

typedef struct __system {
  int dim; // Dimension of the system
  int *itype;
  int nx;
  int ny;
  double Lx, Ly;
  double dx, dy;
  int NumberOfInteractingParticles;
  int MaxNumberOfParticles;
  int ntotal; /* total number of fluid particles */
  int NBoundaries;  /* number of virtual particles */
  int NumberOfVirtualParticles;   /* number of boundary particles */
  double LeftBoundary, RightBoundary, WallVelocity, WallAcceleration;
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
  double *gradVxx; /* components of the velocity gradient*/
  double *gradVxy;
  double *gradVyx;
  double *gradVyy;
  double *strainRatexx; /* components of the strain rate tensor by 2*/
  double *strainRatexy;
  double *strainRateyx;
  double *strainRateyy;
  double *divV;
  double *deviatoricxx; /*deviatoric stress tensor */ 
  double *deviatoricxy;
  double *deviatoricyx;
  double *deviatoricyy;
  double *KinVisc; /* kinematic viscosity*/
  double *numberDensity;
  double **dummyVelocity;
  double *dummyPressure;
  double **dummyAcceleration;
  double **fluctuation;
  double *vorticity;
  double **temps; /* used for checking correctness of dvdt[i][d]*/
  // temporary variables used for the calculations
  double alpha_cubic, alpha_quintic, alpha_wendland;
  double alpha;
  double zeta; /* h = zeta* dx */
  int NearestNeighbor; /* Average number of neighbors per particle */
  double SpeedOfSound;
  int *i_pair, *j_pair;
  int *i_inflow;
  double *rij2, *w;
  double **dwdx, **dwdx_r;
  double *w_r; /* scalar part of grad(w)*/
  double **dwdx;
  double **dwdx_r; /* grad(dwdr/r) for computing the hessian in vectorial form */
  int **Neighbors;
  double **DistanceNeighbors;
  matrix r;
  int numReps;
  rep *ri;
  int NumberOfRepresentatives;
  int inflowParticles;
  double CompressionFactor;
  double AdiabaticConstant;
  int niac;
  unint HowManyNeighbors;
  matrix x,q;
} System;

double **CreateMatrix(int x, int y);
void DestroyMatrix(double **a);
System *CreateSystem(const int dim, 
		     const int nx, 
		     const int ny, 
		     const double Lx, 
		     const double Ly,
		     const double alpha, 
		     const double SpeedOfSound, 
		     const double rho0,
		     const double AdiabaticConstant,
                     const int NearestNeighbor,
		     const int HowManyNeighbors,
		     const int numReps);
void DestroySystem(System *a);
void SearchNeighbors(System *sys);
#endif
