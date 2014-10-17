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
#include "param.h"
#include "System.h"

void output(System *sys, int itime)
{
  FILE *fip;
  int d, i, nsum;
  char filename_xv[30];
  // verify this fprintf statement
  for(i = 1; i <itime; i++)
    {
      sprintf(filename_xv, "f_xv.%d.dat", itime);
    }
	
  nsum = sys->ntotal + sys->nvirt + sys->nvll; /* fluid + boundary + virtual particles */
  fip = fopen("filename_xv","w+");
  nsum = sys->ntotal + sys->nvirt + sys->nvll;
  for(int i = 0; i < nsum; i++)
    {
      for( int d = 0; d < sys->dim; d++)
	{
	  fprintf(fip, "%d\t%e\t%e\t%e\t%e\t%e\n", i, sys->Position[i][d],
		  sys->Velocity[i][d], sys->rho[i], sys->Pressure[i], sys->Energy[i]);
	}
    }
  fclose(fip);
	
}	
	
	
