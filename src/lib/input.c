/*
*********************************************************************
******* Input file *********
* This is the subprogram for the initial input
* Generates or reads from files all the physical 
* properties of fluid particles only					 
*********************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <malloc.h>
#include "param.h"
#include "System.h"

void input(System *sys, double xl, double yl)
{
  FILE *fip, *fop, *fin;
  double space;
  int tmp = 0;
  bool input_file = false;
  /* load initial fluid particles from file */
  if(input_file)
    {
      fip = fopen("ini_txt.dat","a");
      fop = fopen("ini_state","a");
      fin = fopen("ini_other","a");
		
      printf("***************************************************\n");
      printf("Loading initial particle configuration\n");
      printf("total number of particles: %d\n", sys->ntotal);
      printf("***************************************************\n");
      fscanf(fip, "%d", &tmp);
      if(tmp != sys->ntotal)
	printf("Error the files do not correspond to the system you want to simulate\n");

      for(int i = 0; i < sys->ntotal; i++)
	{
	  for(int d = 0; d < sys->dim; d++)
	    {
	      fread(sys->Position[i], sizeof(double), 2, fip);
	      fread(sys->Velocity[i], sizeof(double), 2, fip);
	    }
	}

      fread(sys->mass, sizeof(double), sys->ntotal, fop);
      fread(sys->rho, sizeof(double), sys->ntotal, fop);
      fread(sys->Pressure, sizeof(double), sys->ntotal, fop);
      fread(sys->Energy, sizeof(double), sys->ntotal, fop);
      fread(sys->itype, sizeof(int), sys->ntotal, fin);
      fread(sys->hsml, sizeof(double), sys->ntotal, fin);
    }	
  else
    {
      /* Generate fluid particles in this subprogram */
		
      fip = fopen("ini_xv.dat","w+");
      fop = fopen("ini_state.dat","w+");
      fin = fopen("ini_other.dat","w+");
      
      /* outer grid */
      double count = 0;
      for(int iy = 0; iy < sys->ny; iy++)
	{
	  for(int ix = 0; ix < sys->nx; ix++)
	    {
	      sys->Position[iy*sys->nx+ix][0] = (ix+1)*sys->dx;
	      sys->Position[iy*sys->nx+ix][1] = (iy+1)*sys->dy;
              count++;
	    }
	}
      sys->ntotal = count;
      memset(sys->Velocity[0], 0, sizeof(double)*2*sys->MaxNumberOfParticles);
//       memset(sys->Pressure, 0, sizeof(double)*sys->MaxNumberOfParticles);
      memset(sys->Energy, 0, sizeof(double)*sys->MaxNumberOfParticles);
      /* initial condtions */
      double g = 9.81;
      for(int i = 0; i < sys->ntotal; i++)
	{
	  sys->rho[i] = sys->rho0*pow(1. + (g*sys->rho0/sys->CompressionFactor)*(sys->Ly-sys->Position[i][1]), 1./sys->AdiabaticConstant);
	  sys->mass[i] = sys->rho0*sys->dx*sys->dy; /* density is mass over area in 2D*/
	  sys->hsml[i] = sys->zeta*sys->dx; //0.129;
	  sys->itype[i] = 1;
	  sys->Pressure[i] = sys->rho0*g*(sys->Ly-sys->Position[i][1]);
          sys->KinVisc[i] = 1.0e-6;
	}

      /* write just the derived particle info into file */
      /* fprintf(fip, "%d\n", sys->ntotal); */
      for(int i = 0; i < sys->ntotal; i++)
	{
	  fwrite(sys->Position[i], sizeof(double), 2, fip);
	  fwrite(sys->Velocity[i], sizeof(double), 2, fip);
	  /* fwrite(&sys->mass[i], sizeof(double), 1, fop); */
	  /* fwrite(&sys->rho[i], sizeof(double), 1, fop); */
	  /* fwrite(&sys->Pressure[i], sizeof(double), 1, fop); */
	  /* fwrite(&sys->Energy[i], sizeof(double), 1, fop); */
	  
	}

      fwrite(sys->mass, sizeof(double), sys->ntotal, fop);
      fwrite(sys->rho, sizeof(double), sys->ntotal, fop);
      fwrite(sys->Pressure, sizeof(double), sys->ntotal, fop);
      fwrite(sys->Energy, sizeof(double), sys->ntotal, fop);
      fwrite(sys->itype, sizeof(int), sys->ntotal, fin);
      fwrite(sys->hsml, sizeof(double), sys->ntotal, fin);

      /* fprintf(fip, "%d ", i); */
      /* for(int d = 0; d < sys->dim; d++) */
      /*   { */
      /*     fprintf(fip, "\t%e ", sys->Position[i][d]); */
      /*   } */
      /* for(int d = 0; d < sys->dim; d++) */
      /*   { */
      /*     fprintf(fip, "\t%e ", sys->Velocity[i][d]); */
      /*   } */
	  
      /* fprintf(fip,"\n"); */
      /* fprintf(fop, "%d\t%e\t%e\t%e\t%e\n ", i, sys->mass[i], sys->rho[i], sys->Pressure[i], sys->Energy[i]); */
      /* fprintf(fin, "%d\t%d\t%e\n ", i, sys->itype[i], sys->hsml[i]); */
    /* } */
  printf("****************************************************\n");
  printf("Initial particle info generated\n");
  printf("Total number of particles:%d\n", sys->ntotal);
  printf("****************************************************\n");
}
fclose(fip);
fclose(fop);
fclose(fin);
}	
