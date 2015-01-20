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
  a->ntotal = nx*ny;
  //  a->nvirt = nvirt;
  a->MaxNumberOfParticles = 0.14 * 0.14 * 3.14 / (a->dx*a->dy);
  a->MaxNumberOfParticles += 10;
  a->MaxNumberOfParticles *= 2*a->ntotal;
  a->alpha = alpha;

  a->Position = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->Velocity = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->Velocity_xsph = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->Velocity_min = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->dvdt = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->av = CreateMatrix(a->MaxNumberOfParticles, dim);
  a->dwdx = CreateMatrix(a->MaxNumberOfParticles, dim);
  
  posix_memalign((void **)&a->Energy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->hsml, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->Pressure, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->itype, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->rho, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->rho_min, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->drhodt, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->mass, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->i_pair, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->j_pair, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->rij2, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->w, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->i_inflow, 16, sizeof(int)*a->MaxNumberOfParticles);

  a->alpha_d = 0;
  a->SpeedOfSound = SpeedOfSound;
  a->NumberOfInteractingParticles = 0;
  a->CompressionFactor = (1./AdiabaticConstant) * rho0 * SpeedOfSound * SpeedOfSound;
  a->AdiabaticConstant = AdiabaticConstant;
  a->rho0 = rho0;

  
  // We only need to calculate the neighbors of the particles in the box
  a->HowManyNeighbors = HowManyNeighbors;
  a->Neighbors = CreateMatrixInt(a->ntotal, a->HowManyNeighbors);
  a->DistanceNeighbors = CreateMatrix(a->ntotal, a->HowManyNeighbors);
  a->numReps = numReps;
  
  
  initMat( &a->q, a->ntotal, a->dim);
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
  free(a);
  cleanup();
  
}

void SearchNeighbors(System *sys)
{
  sys->ri = (rep*)calloc( CPAD(sys->numReps), sizeof(*sys->ri) ); //data struct for RBC
  initMat(&sys->x, sys->ntotal+sys->NBoundaries+sys->NumberOfVirtualParticles, sys->dim);
  sys->x.mat = sys->Position[0];
  /* buildExact(sys->x, &sys->r, sys->ri, sys->numReps); */
  /* searchExactK(sys->q, sys->x, sys->r, sys->ri, (unint **)sys->Neighbors, sys->DistanceNeighbors, sys->HowManyNeighbors); */
  /* for(int c=0;c<32;c++) */
  /*   printf("%.5lf ", sys->DistanceNeighbors[0][c]); */
  /* freeRBC(sys->r, sys->ri); */
  //sys->ri = (rep*)calloc( CPAD(sys->numReps), sizeof(*sys->ri) ); //data struct for RBC
  buildOneShot(sys->x, &sys->r, sys->ri, sys->numReps);
  searchOneShotK(sys->q, sys->x, sys->r, sys->ri, (unint **)sys->Neighbors, sys->DistanceNeighbors, sys->HowManyNeighbors);
  freeRBC(sys->r, sys->ri);
  /* printf("\n"); */
  /* for(int c=0;c<32;c++) */
  /*   printf("%.5lf ", sys->DistanceNeighbors[0][c]); */
}
