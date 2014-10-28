#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <string.h>
#include "param.h"
#include "System.h"


double **CreateMatrix(int x, int y)
{
  double **a = (double **)malloc(sizeof(double *)*x);
  posix_memalign(a, 16, sizeof(double)*x*y);
  for(int i=1;i<x;i++)
    a[i] = a[i-1] + y;
  return a;
}

void DestroyMatrix(double **a)
{
  free(a[0]);
  free(a);
}


System *CreateSystem(const int dim, 
		     const int nx, 
		     const int ny, 
		     const int MaxNumberOfParticles, 
		     const double alpha, 
		     const double SpeedOfSound, 
		     const double rho0,
		     const double AdiabaticConstant)
{
  System *a = (System *)malloc(sizeof(System));
  memset(a, 0, sizeof(System));
  a->nx = nx;
  a->ny = ny;
  a->dim = dim;
  a->ntotal = nx*ny;
  //  a->nvirt = nvirt;
  a->MaxNumberOfParticles = MaxNumberOfParticles;
  a->alpha = alpha;
  a->Position = CreateMatrix(MaxNumberOfParticles, dim);
  a->Velocity = CreateMatrix(MaxNumberOfParticles, dim);
  a->Velocity_xsph = CreateMatrix(MaxNumberOfParticles, dim);
  a->Velocity_min = CreateMatrix(MaxNumberOfParticles, dim);
  a->dvdt = CreateMatrix(MaxNumberOfParticles, dim);
  a->av = CreateMatrix(MaxNumberOfParticles, dim);
  a->dwdx = CreateMatrix(MaxNumberOfParticles, dim);
  
  posix_memalign(&a->Energy, 16, sizeof(double)*MaxNumberOfParticles);
  posix_memalign(&a->hsml, 16, sizeof(double)*MaxNumberOfParticles);
  posix_memalign(&a->Pressure, 16, sizeof(double)*MaxNumberOfParticles);
  posix_memalign(&a->itype, 16, sizeof(int)*MaxNumberOfParticles);
  posix_memalign(&a->rho, 16, sizeof(double)*MaxNumberOfParticles);
  posix_memalign(&a->rho_min, 16, sizeof(double)*MaxNumberOfParticles);
  posix_memalign(&a->drhodt, 16, sizeof(double)*MaxNumberOfParticles);
  posix_memalign(&a->mass, 16, sizeof(double)*MaxNumberOfParticles);
  posix_memalign(&a->i_pair, 16, sizeof(int)*MaxNumberOfParticles);
  posix_memalign(&a->j_pair, 16, sizeof(int)*MaxNumberOfParticles);
  posix_memalign(&a->rij2, 16, sizeof(double)*MaxNumberOfParticles);
  posix_memalign(&a->w, 16, sizeof(double)*MaxNumberOfParticles);
  posix_memalign(&a->i_inflow, 16, sizeof(int)*MaxNumberOfParticles);

  a->alpha_d = 0;
  a->SpeedOfSound = SpeedOfSound;
  a->NumberOfInteractingParticles = 0;
  a->CompressionFactor = (1./AdiabaticConstant) * rho0 * SpeedOfSound * SpeedOfSound;
  a->AdiabaticConstant = AdiabaticConstant;
  a->rho0 = rho0;
  return a;
}

void DestroySystem(System *a)
{
  DestroyMatrix(a->Position);
  DestroyMatrix(a->Velocity);
  DestroyMatrix(a->Velocity_xsph);
  DestroyMatrix(a->Velocity_min);
  DestroyMatrix(a->dvdt);
  DestroyMatrix(a->av);
  DestroyMatrix(a->dwdx);
  
  free(a->Energy);
  free(a->hsml);
  free(a->Pressure);
  free(a->itype);
  free(a->rho);
  free(a->rho_min);
  free(a->drhodt);  
  free(a->mass);
  free(a->i_pair);
  free(a->j_pair);
  free(a->rij2);
  free(a->w);
  free(a->i_inflow);
  free(a);
}

