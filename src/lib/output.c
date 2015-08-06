/* 
*****************************************************************************
* subprogram to write data to output file
* writes output in text file with 8 columns
* containing particle number i, 
* position vector components x, y
* velocity field components u, v, 
* particle density rho,
* pressure p, and 
* mass m.  
*****************************************************************************
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "param.h"
#include "output.h"

void output(System *sys, int itime)
{
  FILE *fip;
  int nsum;
  char filename_xv[128];
  // verify this fprintf statement
  memset(filename_xv, 0, sizeof(char)*128);
  sprintf(filename_xv, "f_xv.%d.dat", itime);
  
  nsum = sys->ntotal + sys->NBoundaries + sys->NumberOfVirtualParticles; /* fluid + boundary + virtual particles */
  fip = fopen(filename_xv,"w+");
  for(int i = 0; i < nsum; i++)
    {
      
      fprintf(fip, "%d ",i);
      for(int d = 0; d < sys->dim; d++)
	{
	  fprintf(fip, "\t%.10lf", sys->Position[i][d]);
	  
	}
      for(int d = 0; d < sys->dim; d++)
	{
	  fprintf(fip, "\t%e", sys->Velocity[i][d]);
	}
      fprintf(fip, "\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", sys->rho[i], sys->Pressure[i], sys->mass[i], sys->hsml[i], sys->temps[i][0], sys->temps[i][1]);
    }
  fclose(fip);	
}	
	
	
