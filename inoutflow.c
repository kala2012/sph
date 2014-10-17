/*
****************************************************************
* This subprogram handles inflow and outflow of particles
* inflow particles are generated at inflow boundary and 
* added to the existing arrays. At the outflow boundary, 
* particles are removed. Only horizontal velocity component
* is possible at inflow.
****************************************************************
*/
#include <stdio.h>
#include "param.h"
#include "System.h"

void inoutflow(System *sys)
{

  int ii;
  double inputmass, inputvol;

  /*Left inflow*/
  ii = sys->ntotal-1;
  inputvol = 0.0;
  inputmass = 0.0;
  for(int k = 0; k < sys->inflowParticles; k++)
    {
      int i = sys->i_inflow[k];
      if(sys->Position[i][1] > 0.1)
	{
	  sys->Position[ii][0] = sys->Position[i][0] - 0.1;
	  sys->Position[ii][1] = sys->Position[i][1];
	  sys->Velocity[ii][0] = sys->Velocity[i][0];
	  sys->Velocity[ii][1] = 0.0;
	  sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
	  sys->Velocity_xsph[1][ii] = 0.0;
	  sys->rho[ii] = sys->rho[i];
	  sys->mass[ii] = sys->mass[i];
	  sys->Pressure[ii] = sys->Pressure[i];
	  sys->Energy[ii] = sys->Energy[i];
	  sys->hsml[ii] = sys->hsml[i];
	  sys->rho_min[ii] = sys->rho_min[i];
	  sys->Velocity_min[ii][0] = sys->Velocity_min[i][0];
	  sys->Velocity_min[ii][1] = 0.0;
	  /*Track input discharge*/
	  inputmass += sys->mass[ii];
	  inputvol += sys->mass[ii]/sys->rho[ii];
	  ii +=1;
	}
    }

  if(inputmass != 0)
    {
      printf("input mass: %lf kg/m\t and input vol: %lf m^2\n", inputmass, inputvol);
    }

  /*Right boundary*/
  /* review this because it fails badly */ 

  int ntotal = ii;
  for(int i = 0; i < ntotal; i++)
    {
      if(sys->Position[i][1] > 50.0)
	{
	  sys->Position[i][0] = sys->Position[ntotal][0];
	  sys->Position[i][1] = sys->Position[ntotal][1];
	  sys->Velocity[i][0] = sys->Velocity[ntotal][0];
	  sys->Velocity[i][1] = sys->Velocity[ntotal][1];
	  sys->Velocity_xsph[i][0] = sys->Velocity_xsph[ntotal][0];
	  sys->Velocity_xsph[i][1] = sys->Velocity_xsph[ntotal][1];
	  sys->mass[i] = sys->mass[ntotal];
	  sys->rho[i] = sys->rho[ntotal];
	  sys->Pressure[i] = sys->Pressure[ntotal];
	  sys->Energy[i] = sys->Energy[ntotal];
	  sys->hsml[i] = sys->hsml[ntotal];
	  sys->rho_min[i] = sys->rho_min[ntotal];
	  sys->Velocity_min[i][0] = sys->Velocity_min[0][ntotal];
	  sys->Velocity_min[i][1] = sys->Velocity_min[1][ntotal];
	  sys->ntotal--;
	}
    }
}	
