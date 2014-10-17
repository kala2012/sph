/*
***********************************************************
* Input Boundary Particles
* This subprogram generates boundary particles.
***********************************************************
*/
#include "inputbp.h"

void inputbp(System *sys, double xl, double vl, int itime, double dt)
{
  FILE *fip, *fop, *fin;
  int i, j, d;
  double dx; /* dx = space*/
  bool vp_file = false;
  /* load boundary particles from file */
  if(vp_file)
    {
      fip = fopen("xv_vp.dat","r");
      fop = fopen("ini_state.dat","r");
      fin = fopen("ini_other.dat","r");
      fscanf(fip, "%d", sys->nvirt);
      
      for(int j = 0; j < sys->nvirt; j++)
		{
		  i = j + sys->ntotal;
		  for(int d = 0; d < sys->dim; d++)
		    {
		      fprintf(fip, "%d\t%e\t%e\n ", i, sys->Position[i][d], sys->Velocity[i][d]);
		    }
		  fprintf(fop, "%d\t%e\t%e\t%e\t%e\n ", i, sys->mass[i], 
			  sys->rho[i], sys->Pressure[i], sys->Energy[i]);
		  fprintf(fin, "%d\t%d\t%e\n ", i, sys->itype[i], sys->hsml[i]);
		}
      fclose(fip);
      fclose(fop);
      fclose(fin);
      /* if input of boundary particles is not from file,
       * then generate boundary particles here */
    } 
  else
    {
      sys->nvirt = 0;
      dx = 0.1;
      
      /* boundary particles in bottom boundary */
      for(int i = 0; i < 120; i++)
	{
	  sys->Position[sys->ntotal+sys->nvirt+i][0] = (i-1)*dx;
	  sys->Position[sys->ntotal+sys->nvirt+i][1] = 0.0;
	  sys->Velocity[sys->ntotal+sys->nvirt+i][0] = 0.0;
	  sys->Velocity[sys->ntotal+sys->nvirt+i][1] = 0.0;
	}
      /* boundary particles on left side */
      for(int i = 0; i < 100; i++)
	{
	  sys->Position[sys->ntotal+sys->nvirt+i][0] = 0.2575 - 0.2575*cos(2.445*itime*dt);
	  sys->Position[sys->ntotal+sys->nvirt+i][1] = i*dx;
	  sys->Velocity[sys->ntotal+sys->nvirt+i][0] = 0.2575*2.445*sin(2.445*itime*dt);
	  sys->Velocity[sys->ntotal+sys->nvirt+i][1] = 0.0;
	  sys->nvirt ++;
	  
	}
      xl = 0.2575 - 0.2575*cos(2.445*itime*dt);
      vl = 0.2575*2.445*sin(2.445*itime*dt);
      /* boundary particles on the right side */
      for(int i = 0; i < 100; i++)
	{
	  sys->Position[sys->ntotal+sys->nvirt][0] = 0.2575 - 0.2575*cos(2.445*itime*dt) + 10.1;
	  sys->Position[sys->ntotal+sys->nvirt][1] = i*dx;
	  sys->Velocity[sys->ntotal+sys->nvirt][0] = 0.2575*2.445*sin(2.445*itime*dt);
	  sys->Velocity[sys->ntotal+sys->nvirt][1] = 0.0;
	  sys->nvirt++;
	}
      xl = 0.2575 - 0.2575*cos(2.445*itime*dt);
      vl = 0.2575*2.445*sin(2.445*itime*dt);
      
      /* initialize fluid variables */
      for(int i = 0; i < sys->nvirt; i++)
	{
	  sys->rho[sys->ntotal+i] = 1000.0;
	  sys->mass[sys->ntotal+i] = sys->rho[sys->ntotal+i]*dx*dx;
	  sys->Pressure[sys->ntotal+i] = 0.0;
	  sys->Energy[sys->ntotal+i] = 200.0;
	  sys->itype[sys->ntotal+i] = -2;
	  sys->hsml[sys->ntotal+i] = 0.14;
	}
    }
  /* write data virt part to file */
  if(itime == 1)
    {
      fip = fopen("xv_vp.dat","w+");
      fop = fopen("state_vp.dat","w+");
      fin = fopen("other_vp.dat","w+");
      fscanf(fip, "%d", sys->nvirt);
      for(int i = sys->ntotal; i <sys->ntotal+sys->nvirt; i++)
	{
	  for(int d = 0; d < sys->dim; d++)
	    {
	      fprintf(fip, "%d\t%e\t%e\n ", i, sys->Position[i][d], sys->Velocity[i][d]);
	    }
	  fprintf(fop, "%d\t %e\t%e\t%e\t%e\n ", i, sys->mass[i], 
		  sys->rho[i], sys->Pressure[i], sys->Energy[i]);
	  fprintf(fin, "%d\t %d\t%e\n ", i, sys->itype[i], sys->hsml[i]);	
	}
      fclose(fip);
      fclose(fop);
      fclose(fin);
    }
}	
