/*
***********************************************************
* Input Boundary Particles
* This subprogram generates boundary particles.
***********************************************************
*/
#include "inputbp.h"

void inputbp(System *sys, double xl, double vl, double t)
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
      fscanf(fip, "%d", &sys->NBoundaries);
      for(int i = sys->ntotal; i <sys->ntotal+sys->NBoundaries; i++)
  	{
	  fread(sys->Position[i], sizeof(double), sys->dim, fip);
	  fread(sys->Velocity[i], sizeof(double), sys->dim, fip);
	}
      fread(sys->mass+sys->ntotal, sizeof(double), sys->NBoundaries, fop);
      fread(sys->rho+sys->ntotal, sizeof(double), sys->NBoundaries, fop);
      fread(sys->Pressure+sys->ntotal, sizeof(double), sys->NBoundaries, fop);
      fread(sys->Energy+sys->ntotal, sizeof(double), sys->NBoundaries, fop);
      
      fread(sys->itype+sys->ntotal, sizeof(int), sys->NBoundaries, fin);
      fread(sys->hsml+sys->ntotal, sizeof(double), sys->NBoundaries, fin);

      /* for(int j = 0; j < sys->NBoundaries; j++) */
      /* 	{ */
      /* 	  i = j + sys->ntotal; */
      /* 	  for(int d = 0; d < sys->dim; d++) */
      /* 	    { */
      /* 	      fprintf(fip, "%d\t%e\t%e\n ", i, sys->Position[i][d], sys->Velocity[i][d]); */
      /* 	    } */
      /* 	  fprintf(fop, "%d\t%e\t%e\t%e\t%e\n ", i, sys->mass[i],  */
      /* 		  sys->rho[i], sys->Pressure[i], sys->Energy[i]); */
      /* 	  fprintf(fin, "%d\t%d\t%e\n ", i, sys->itype[i], sys->hsml[i]); */
      /* 	} */
      fclose(fip);
      fclose(fop);
      fclose(fin);
      /* if input of boundary particles is not from file,
       * then generate boundary particles here */
    } 
  else
    {
      sys->NBoundaries = 0;
      dx = 0.1;
      /* t = 0; */  
      sys->LeftBoundary = 0.2575*(1.- cos(2.445*t));
      /* boundary particles in bottom boundary */
      for(int i = -10; i < 112; i++)
	{
	  sys->Position[sys->ntotal+sys->NBoundaries][0] = i*dx + sys->LeftBoundary;
	  sys->Position[sys->ntotal+sys->NBoundaries][1] = 0.0;
	  sys->Velocity[sys->ntotal+sys->NBoundaries][0] = 0.0;
	  sys->Velocity[sys->ntotal+sys->NBoundaries][1] = 0.0;
	  sys->NBoundaries++;
	}
      /* boundary particles on left side */
      /* sys->LeftBoundary = 0.2575 - 0.2575*cos(2.445*t); */
      for(int i = -10; i < 100; i++)
	{
	  sys->Position[sys->ntotal+sys->NBoundaries][0] = sys->LeftBoundary;
	  sys->Position[sys->ntotal+sys->NBoundaries][1] = (i+1)*dx;
	  sys->Velocity[sys->ntotal+sys->NBoundaries][0] = 0.2575*2.445*sin(2.445*t);
	  sys->Velocity[sys->ntotal+sys->NBoundaries][1] = 0.0;
	  sys->NBoundaries ++;	  
	}
      xl = 0.2575 - 0.2575*sin(2.445*t);
      vl = 0.2575*2.445*cos(2.445*t);
      /* boundary particles on the right side */
      sys->RightBoundary = sys->LeftBoundary + 10.1;
      for(int i = -10; i < 100; i++)
	{
	  sys->Position[sys->ntotal+sys->NBoundaries][0] = sys->RightBoundary;
	  sys->Position[sys->ntotal+sys->NBoundaries][1] = (i+1)*dx;
	  sys->Velocity[sys->ntotal+sys->NBoundaries][0] = 0.2575*2.445*sin(2.445*t);
	  sys->Velocity[sys->ntotal+sys->NBoundaries][1] = 0.0;
	  sys->NBoundaries++;
	}
      xl = 0.2575 - 0.2575*cos(2.445*t) + 10.1;
      vl = 0.2575*2.445*sin(2.445*t);
      
      /* initialize fluid variables */
      for(int i = 0; i < sys->NBoundaries; i++)
	{
	  sys->rho[sys->ntotal+i] = 1000.0;
	  sys->mass[sys->ntotal+i] = sys->rho[sys->ntotal+i]*dx*dx;
	  sys->Pressure[sys->ntotal+i] = 0.0;
	  sys->Energy[sys->ntotal+i] = 357.1;
	  sys->itype[sys->ntotal+i] = -2;
	  sys->hsml[sys->ntotal+i] = 0.14;
	}
    }
  /* write data virt part to file */
  // FIX THIS !!!!!!!!!!

  if(t ==0 )
    {
      fip = fopen("xv_vp.dat","w+");
      fop = fopen("state_vp.dat","w+");
      fin = fopen("other_vp.dat","w+");
      fprintf(fip, "%d\n", sys->NBoundaries); 
      for(int i = sys->ntotal; i <sys->ntotal+sys->NBoundaries; i++)
  	{
	  fwrite(sys->Position[i], sizeof(double), sys->dim, fip);
	  fwrite(sys->Velocity[i], sizeof(double), sys->dim, fip);
	}
      fwrite(sys->mass+sys->ntotal, sizeof(double), sys->NBoundaries, fop);
      fwrite(sys->rho+sys->ntotal, sizeof(double), sys->NBoundaries, fop);
      fwrite(sys->Pressure+sys->ntotal, sizeof(double), sys->NBoundaries, fop);
      fwrite(sys->Energy+sys->ntotal, sizeof(double), sys->NBoundaries, fop);
      
      fwrite(sys->itype+sys->ntotal, sizeof(int), sys->NBoundaries, fin);
      fwrite(sys->hsml+sys->ntotal, sizeof(double), sys->NBoundaries, fin);
      fclose(fip);
      fclose(fop);
      fclose(fin);
    }
}	
