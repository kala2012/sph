/* ******************************************************
 * This is the main program.
 * It calls " input.c" once, starts time loop with 
 * "derivatives.c" and leap frog time integration every
 * time step. Calls "output.c" if requested and then 
 * generates output on the screen.
 * ****************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <inttypes.h>
#include "param.h"
#include "input.h"
#include "derivatives.h"
#include "output.h"
#include "System.h"

int main(int argc, char **argv)
{
  const int dim=2;
  double dt = 1e-3; /* set timestep */
  unsigned int maxtimestep = 139; /* set maximum timestep */
  int printstep = 1;
  int screenstep = 1;
  int itime; /* iteration time */
  int nx = 100; /* number of particles in x-dirextion */
  int ny = 50;  /* number of particles in y-direction */

  if(argc == 1)
    {
      printf("Give a maximum number of time step\n");
      exit(1);
    }
  maxtimestep= atoi(argv[1]);
  // parameters for water
  System *sys = CreateSystem(dim, nx, ny, 10, 5, 0.01, 50.0, 1000, 7, 64, 200);
  double xl = 1.0;
  double yl = 1.0;

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  input(sys, xl, yl); 

  //SearchNeighbors(sys);
 
  
  
  memset(sys->dvdt[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  memset(sys->drhodt, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memcpy(sys->rho_min, sys->rho, sizeof(double)*sys->ntotal);
  
  for(itime = 1; itime <= maxtimestep; itime++)
    {
      double t = itime * dt;
      if(itime!=1) {
      	memcpy(sys->rho_min, sys->rho, sizeof(double)*sys->ntotal);
	
  	/* cblas_daxpy(sys->ntotal, 0.5*dt, sys->drhodt, 1, sys->rho, 1); */
  	/* memcpy(sys->Velocity_min[0], sys->Velocity[0], sizeof(double)*2*sys->ntotal); */
  	/* cblas_daxpy(sys->ntotal, 0.5*dt, sys->dvdt[0], 1, sys->Velocity_xsph[0], 1); */
  	/* cblas_daxpy(sys->ntotal, 0.5*dt, sys->dvdt[0], 1, sys->Velocity[0], 1); */
	
      	for(int i=0;i<sys->ntotal ;i++) {
      	  sys->rho[i] += 0.5*dt*sys->drhodt[i];
      	  for(int d=0;d<sys->dim;d++) {
      	     sys->Velocity_min[i][d] = sys->Velocity[i][d];
      	    sys->Velocity_xsph[i][d] += 0.5 * dt * sys->dvdt[i][d];
      	    sys->Velocity[i][d] +=  0.5 * dt * sys->dvdt[i][d];
      	  }
      	}	
      }
      
      if(itime==1) {
  	memcpy(sys->Velocity_xsph[0], sys->Velocity[0], sizeof(double)*sys->dim*sys->MaxNumberOfParticles);
  	memcpy(sys->Velocity_min[0], sys->Velocity[0], sizeof(double)*sys->dim*sys->MaxNumberOfParticles);
      }

      derivatives(sys, t);

      
      
      if(itime==1) {
      	for(int i=0;i<sys->ntotal ;i++) {
      	  sys->rho[i] += 0.5*dt*sys->drhodt[i];
      	  for(int d=0;d<sys->dim;d++) {
      	    sys->Velocity_xsph[i][d] = sys->Velocity[i][d] + 0.5 * dt * sys->dvdt[i][d] + sys->av[i][d];
      	    sys->Velocity[i][d] +=  0.5 * dt * sys->dvdt[i][d];
      	    sys->Position[i][d] +=  dt * sys->Velocity_xsph[i][d];
      	  }
      	}
      } else {
      	for(int i=0;i<sys->ntotal ;i++) {
      	  sys->rho[i] = sys->rho_min[i] + dt*sys->drhodt[i];
      	  for(int d=0;d<sys->dim;d++) {
      	    sys->Velocity_xsph[i][d] = sys->Velocity_min[i][d] + dt * sys->dvdt[i][d] + sys->av[i][d];
      	    sys->Velocity[i][d] =  sys->Velocity_min[i][d] + dt * sys->dvdt[i][d];
      	    sys->Position[i][d] +=  dt * sys->Velocity_xsph[i][d];
      	  }
      	}
      }
      if( (itime % printstep) == 0)
  	{
  	  printf("%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
  		 sys->Velocity[0][0],
  		 sys->Velocity[0][1],
  		 sys->dvdt[0][0],
  		 sys->dvdt[0][1],
  		 sys->Pressure[0],
  		 sys->rho[0],
  		 sys->drhodt[0]);
  	  output(sys, itime);
  	}
      if( (itime % screenstep) == 0)
  	{
  	  /* write useful info to the console */
  	  printf("******************************************\n");
	  
  	  /* verify %lld using with mingw's gcc */
  	  printf("\nTimestep %d of  %d\n", itime, maxtimestep);
  	  printf("Interaction pairs: %d\n", sys->NumberOfInteractingParticles);
  	  printf("******************************************\n");
  	}
    }

  end = clock();
  cpu_time_used = (double)(end - start)/CLOCKS_PER_SEC;
  printf("CPU time %lf s\n", cpu_time_used);
  
  DestroySystem(sys);
  return 0;
}
