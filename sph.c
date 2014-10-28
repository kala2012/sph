/* ******************************************************
 * This is the main program.
 * It calls " input.c" once, starts time loop with 
 * "derivatives.c" and leap frog time integration every
 * time step. Calls "output.c" if requested and then 
 * generates output on the screen.
 * ****************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <inttypes.h>
#include "param.h"
#include "input.h"
#include "derivatives.h"
#include "output.h"
#include "System.h"

int main(void)
{
  const int dim=2;
  double dt = 1e-3; /* set timestep */
  unsigned long long int maxtimestep = 10; /* set maximum timestep */
  int printstep = 1;
  int screenstep = 1;
  int itime; /* iteration time */
  int niac = 0; /* number of interacting particles */
  int nx = 100; /* number of particles in x-dirextion */
  int ny = 50;  /* number of particles in y-direction */
  
  // parameters for water
  System *sys = CreateSystem(dim, nx, ny, 10*maxn, 0.01, 50.0, 1000, 7);
  double xl = 1.0;
  double yl = 1.0;

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  input(sys, xl, yl); 
  derivatives(sys, 1, dt);
  /* for(itime = 1; itime <= maxtimestep; itime++) */
  /*   { */
  /*     //      do{ */
  /* 	for(int i=0; i < sys->ntotal; i++) */
  /* 	  { */
  /* 	    sys->rho_min[i] = sys->rho[i]; */
  /* 	    sys->rho[i] += 0.5*dt*sys->drhodt[i]; */
  /* 	    for(int d = 0; d < sys->dim; d++) */
  /* 	      { */
  /* 		sys->Velocity_min[i][d]  = sys->Velocity[i][d]; */
  /* 		sys->Velocity_min[i][d] += 0.5*dt*sys->dvdt[i][d]; */
  /* 		sys->Velocity_xsph[i][d] += 0.5*dt*sys->dvdt[i][d]; */
  /* 	      } */
  /* 	  } */
  /* 	// }while(itime != 1); */
  /*     /\* update the fluid variables *\/ */
  /*     if(itime >= 1) */
  /* 	{ */
  /* 	  derivatives(sys, itime, dt);	 */
  /* 	} */
  /*     if(itime == 1) */
  /* 	{ */
  /* 	  for(int i = 0; i < sys->ntotal; i++) */
  /* 	    { */
  /* 	      sys->rho[i] += 0.5*dt*sys->drhodt[i]; */
  /* 	      for(int d = 0; d < dim; d++) */
  /* 		{ */
  /* 		  sys->Velocity_xsph[i][d] += 0.5*dt*sys->dvdt[i][d] + sys->av[i][d]; */
  /* 		  sys->Velocity[i][d] += 0.5*dt*sys->dvdt[i][d]; */
  /* 		  sys->Position[i][d] += dt*sys->Velocity_xsph[i][d]; */
  /* 		} */
  /* 	    } */
  /* 	} */
  /*     else */
  /* 	{ */
  /* 	  for(int i = 0; i < sys->ntotal; i++) */
  /* 	    { */
  /* 	      sys->rho[i] = sys->rho_min[i] + dt*sys->drhodt[i]; */
  /* 	      for(int d = 0; d < sys->dim; d++) */
  /* 		{ */
  /* 		  sys->Velocity_xsph[i][d] = sys->Velocity_min[i][d] + dt*sys->dvdt[i][d] + sys->av[i][d]; */
  /* 		  sys->Velocity[i][d] = sys->Velocity_min[i][d] + dt*sys->dvdt[i][d]; */
  /* 		  sys->Position[i][d] += dt*sys->Velocity_xsph[i][d]; */
  /* 		}		 */
  /* 	    }	    */
  /* 	} */
  /*     if( (itime % printstep) == 0) */
  /* 	{ */
  /* 	  output(sys, itime); */
  /* 	} */
  /*     if( (itime % screenstep) == 0) */
  /* 	{ */
  /* 	  /\* write useful info to the console *\/ */
  /* 	  printf("******************************************\n"); */
		    
  /* 	  /\* verify %lld using with mingw's gcc *\/ */
  /* 	  printf("\nTimestep %d of  %"PRIu64"\n", itime, maxtimestep); */
  /* 	  printf("Interaction pairs: %d\n", niac); */
  /* 	  printf("******************************************"); */
  /* 	} */
  /*   } */
  /* end = clock(); */
  /* cpu_time_used = (double)(end - start)/CLOCKS_PER_SEC; */
  /* printf("CPU time %lf s\n", cpu_time_used); */
	
  DestroySystem(sys);
  return 0;
}
