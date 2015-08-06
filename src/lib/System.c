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
  posix_memalign((void **)a, 16, sizeof(double)*x*y);
    
  for(int i=1;i<x;i++)
    a[i] = a[i-1] + y;
  return a;
}



int **CreateMatrixInt(int x, int y)
{
  int **a = (int **)malloc(sizeof(int *)*x);
  posix_memalign((void **)a, 16, sizeof(int)*x*y);
    
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
		     const double Lx, 
		     const double Ly,
		     const double alpha, 
		     const double SpeedOfSound, 
		     const double rho0,
		     const double AdiabaticConstant,
                     const int NearestNeighbor,
		     const int HowManyNeighbors,
		     const int numReps)
{
  System *a = (System *)malloc(sizeof(System));
  memset(a, 0, sizeof(System));
  a->nx = nx;
  a->ny = ny;
  a->Lx = Lx;
  a->Ly = Ly;
  a->dx = Lx / nx;
  a->dy = Ly / ny;
  a->dim = dim;
//   a->ntotal = nx*ny;
  //  a->nvirt = nvirt;
  a->MaxNumberOfParticles = 0.1995 * 0.1995 * 3.14 / (a->dx*a->dy);
  a->MaxNumberOfParticles += 10;
  a->MaxNumberOfParticles *= 2*a->nx*ny;
  a->alpha = alpha;
  a->NearestNeighbor = NearestNeighbor;
  a->zeta = sqrt(NearestNeighbor/(4*M_PI));
  
  a->Position = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->Velocity = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->Velocity_xsph = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->Velocity_min = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->dvdt = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->av = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->dwdx = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->dwdx_r = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->dummyVelocity = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->dummyAcceleration = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->fluctuation = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->temps = CreateMatrix(a->MaxNumberOfParticles, dim);
  
  posix_memalign(&a->Energy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->hsml, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->Pressure, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->itype, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign(&a->rho, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->rho_min, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->drhodt, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->mass, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->i_pair, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign(&a->j_pair, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign(&a->rij2, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->w, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->i_inflow, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign(&a->gradVxx, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->gradVxy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->gradVyx, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->gradVyy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->strainRatexx, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->strainRatexy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->strainRateyx, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->strainRateyy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->deviatoricxx, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->deviatoricxy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->deviatoricyx, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->deviatoricyy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->divV, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->KinVisc, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->numberDensity, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->dummyPressure, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign(&a->vorticity, 16, sizeof(double)*a->MaxNumberOfParticles);
  
  a->alpha_cubic = 0;
  a->alpha_quintic = 0;
  a->alpha_wendland = 0;
  a->SpeedOfSound = SpeedOfSound;
  a->NumberOfInteractingParticles = 0;
  a->CompressionFactor = (1./AdiabaticConstant) * rho0 * SpeedOfSound * SpeedOfSound;
  a->AdiabaticConstant = AdiabaticConstant;
  a->rho0 = rho0;

  
  // We only need to calculate the neighbors of the particles in the box
  a->HowManyNeighbors = HowManyNeighbors;
  a->Neighbors = CreateMatrixInt(a->ntotal, a->HowManyNeighbors);
  a->DistanceNeighbors = CreateMatrixInt(a->ntotal, HowManyNeighbors);
  a->numReps = numReps;
  a->ri = (rep*)calloc( CPAD(a->numReps), sizeof(*a->ri) ); //data struct for RBC
  
  initMat( &a->q, a->ntotal, 2);
  a->q.mat=a->Position[0];
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
  DestroyMatrix(a->dwdx_r);
  DestroyMatrix(a->dummyVelocity);
  DestroyMatrix(a->dummyAcceleration);
  DestroyMatrix(a->fluctuation);
  DestroyMatrix(a->temps);
  DestroyMatrix(a->DistanceNeighbors);
  DestroyMatrix((double **)a->Neighbors);
  
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
  free(a->gradVxx);
  free(a->gradVxy);
  free(a->gradVyx);
  free(a->gradVyy);
  free(a->strainRatexx);
  free(a->strainRatexy);
  free(a->strainRateyx);
  free(a->strainRateyy);
  free(a->deviatoricxx);
  free(a->deviatoricxy);
  free(a->deviatoricyx);
  free(a->deviatoricyy);
  free(a->divV);
  free(a->KinVisc);
  free(a->numberDensity);
  free(a->dummyPressure);
  free(a->vorticity);
  free(a);
  cleanup();
  
}

void SearchNeighbors(System *sys)
{
  initMat( &sys->x, sys->ntotal+sys->NBoundaries+sys->NumberOfVirtualParticles, 2);
  sys->x.mat = sys->Position[0];
  buildOneShot(sys->x, &sys->r, sys->ri, sys->numReps);
  searchOneShotK(sys->q, sys->x, sys->r, sys->ri, sys->Neighbors, sys->DistanceNeighbors, sys->HowManyNeighbors);
}
