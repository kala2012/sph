/* ****************************************************
 * derivative solver
 * the hydrodynamic equation are solved every time step
 * virtual particles are also generated within 2h away
 * from the boundary.
 * makes function calls to inoutflow.h and inputbp.h
 * interaction pairs, kernel W & grad(W) and drho/dt
 * are found.
 * use equation of state to compute pressure p(rho).
 * compute pressure term, gravity, artificial viscous
 * force and boundary forces.
 * compute dv/dt.
 * influence of average velocity calculated with xSPH.
 * at first timestep a function call to parameterfile.h
 * is made.
 *******************************************************
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "param.h"
#include "inoutflow.h"
#include "inputbp.h"
#include "parameterfile.h"
#include "derivatives.h"
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

void CalculateKernel(System *sys)
{
  double dx[3];
  int niac = 0;
  for(int i = 0; i < sys->ntotal; i++) // loop over all particles
    {
      for(int j = i+1; j < (sys->ntotal+sys->NBoundaries+sys->NumberOfVirtualParticles); j++)
	{
	  int check = (sys->Position[i][0]>= 0.1)&&(j>(sys->ntotal+sys->NBoundaries-1));
	  if(!check) {
	    for(int d = 0; d < sys->dim; d++)
	      {
		dx[d] = sys->Position[i][d] - sys->Position[j][d]; /*compute distance between 2 particles*/
	      }
	    double r2 = dx[0]*dx[0] + dx[1]*dx[1]; /*square the initial separation dist*/
	    double h = 0.5*(sys->hsml[i]+sys->hsml[j]); /* average smoothing length*/
	    if(r2 <= 4.*h*h)
	      {
		sys->i_pair[niac] = i;
		sys->j_pair[niac] = j;
		sys->rij2[niac] = r2;
		/*kernel is piecewise cubic spline*/
		double r = sqrt(r2);
		double q = r/h;
		sys->alpha_d = 10.0/(7.*M_PI*h*h); /*2D normalization const for spline*/
		
		/*define cubic spline*/
		if(q >= 0.0 && q < 1.0)
		  {
		    double q2 = q*q;
		    double q3 = q*q*q;
		    sys->w[niac] = sys->alpha_d*(1- 1.5*q2 + 0.75*q3);
		    for(int d = 0; d < sys->dim; d++)
		      {
			sys->dwdx[niac][d] = sys->alpha_d*(-3 + 2.25*q)*dx[d]/(h*h);
		      }
		  }
		else if(q >= 1.0 && q < 2.0)
		  {
		    double qp2 = (2-q)*(2-q);
		    double qp3 = qp2 * (2-q);
		    sys->w[niac] = sys->alpha_d*0.25*qp3;
		    for(int d = 0; d < sys->dim; d++)
		      {
			sys->dwdx[niac][d] = -sys->alpha_d*0.75*qp2*dx[d]/(h*r);
		      }
		  }
		else
		  {
		    sys->w[niac] = 0;
		    for(int d = 0; d < sys->dim; d++)
		      {
			sys->dwdx[niac][d] = 0.0;
		      }
		  }
		niac++;
	      }
	  }
	}
    }
  sys->NumberOfInteractingParticles = niac;
}

void MomentumEquations(System *sys)
{
  const double g = 9.81;
  for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
      int i = sys->i_pair[k];
      int j = sys->j_pair[k];
      double RR = 0.0;
      double f = 0.0;

      double ppart = sys->Pressure[i]/(sys->rho[i]*sys->rho[i]) + sys->Pressure[j]/(sys->rho[j]*sys->rho[j]);
      if(virt_pres) /*virtual pressure for removing tensile instability*/
	{             /*happens if p[i]<0 and p[j]<0*/
	  f = sys->w[k]/(0.5*sys->alpha_d);
	  if(sys->Pressure[i] < 0)
	    RR = -0.2*sys->Pressure[i]/(sys->rho[i]*sys->rho[i]);

	  if(sys->Pressure[j] < 0)
	    RR = RR - 0.2*sys->Pressure[j]/(sys->rho[j]*sys->rho[j]);

	  if(sys->Pressure[i] > 0 && sys->Pressure[j] > 0)
	    RR = 0.01*ppart;
	}
      // f^2
      f *= f;
      // f^4
      f *= f;
      for(int d = 0; (d < sys->dim)&&(i<sys->ntotal); d++)
	sys->dvdt[i][d] -= sys->mass[j]*(ppart + RR*f)*sys->dwdx[k][d];
      for(int d = 0; (d < sys->dim)&&(j<sys->ntotal); d++)
	sys->dvdt[j][d] += sys->mass[i]*(ppart + RR*f)*sys->dwdx[k][d];
    }
  for(int i=0;i<sys->ntotal;i++)
    sys->dvdt[i][1] -= g;
}

/*Leonnard- Jones boundary force*/

 void LeonnardJonesBoundaryForces(System *sys)
{
  const double r0 = 0.1; /*initial particle spacing*/
  const double DD = 10.0;
  double dx[2];
  for(int j = sys->ntotal; j < sys->ntotal + sys->NBoundaries; j++)
    {
      for(int i = 0; i < sys->ntotal; i++)
	{
	  for(int d = 0; d < sys->dim; d++)
	    {
	      dx[d] = sys->Position[i][d] - sys->Position[j][d];
	    }
	  double r2 = dx[0]*dx[0] + dx[1]*dx[1];

	  if(r2 < r0*r0)
	    {
	      double r = r0/sqrt(r2);
	      r *= r; // r^2
	      r *= r; // r^4
	      double pb = DD*r*(r*r-1.)/r2;
	      for(int d = 0; d < sys->dim; d++)
		{
		  sys->dvdt[i][d] += pb*dx[d];
		}
	    }
	}
    }
}

void DensityEquation(System *sys)
{
  int niac = sys->NumberOfInteractingParticles;
  /* solve the density equation */
  for(int k = 0; k < sys->NumberOfInteractingParticles; k++) // niac = Total number of interacting particles
    {
      int i = sys->i_pair[k]; /* interaction pair : particle i and particle j*/
      int j = sys->j_pair[k];
      for(int d = 0; d < sys->dim; d++)
	{
	  double vijdwdx = (sys->Velocity_xsph[i][d] - sys->Velocity_xsph[j][d])*sys->dwdx[k][d];
	  sys->drhodt[i] += sys->mass[j]*vijdwdx;
	  if(j<sys->ntotal)
	    sys->drhodt[j] += sys->mass[i]*vijdwdx; /* double change in sign = OK */
	}
    }
}

 void EquationOfState(System *sys)
{
  for(int i = 0; i < (sys->ntotal + sys->NBoundaries); i++)
    {
      const double rho = sys->rho[i]/sys->rho0;
      sys->Pressure[i] = sys->CompressionFactor*(pow(rho, sys->AdiabaticConstant) - 1.);
    }
}

 void ArtificialViscosity (System *sys)
{
  for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
      int i = sys->i_pair[k];
      int j = sys->j_pair[k];
      if(j < sys->ntotal)
	{
	  double vijrij = 0.0;
	  double mrho = 0.5*(sys->rho[i] + sys->rho[j]);
	  for(int d = 0; d < sys->dim; d++)
	    {
	      vijrij += (sys->Velocity[i][d] - sys->Velocity[j][d])*(sys->Position[i][d] - sys->Position[j][d]);
	    }
	  if(vijrij < 0.0) /*only switch on artificial viscosity if 2*/
	    {                /* interacting particles are approaching each other*/
	      double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
	      double eta2 = 0.01*h*h;
	      double nu = h*vijrij/(sys->rij2[k] + eta2);
	      double artvisc = -sys->alpha*sys->SpeedOfSound*nu/mrho;
	      for(int d = 0; d < sys->dim; d++)
		sys->dvdt[i][d] -= sys->mass[j]*artvisc*sys->dwdx[k][d];
	      
	      for(int d = 0; (d < sys->dim)&&(j<sys->ntotal); d++)
		sys->dvdt[j][d] += sys->mass[i]*artvisc*sys->dwdx[k][d];
	    }

	}
    }
}

 void AverageVelocity(System *sys)
{
  const int niac = sys->NumberOfInteractingParticles;
  const double epsilon = 0.5;
  for(int k = 0; k < niac; k++)
    {
      int i = sys->i_pair[k];
      int j = sys->j_pair[k];
      if(j < sys->ntotal)
	{
	  double mrho = 0.5*(sys->rho[i]+sys->rho[j]);
	  for(int d = 0; d < sys->dim; d++)
	    {
	      // be careful it might be m_i v_i - m_j v_j
	      // luckely m_i = m_j
	      double avgpart = epsilon*(sys->Velocity_xsph[j][d]-sys->Velocity_xsph[i][d])*sys->w[k]/mrho;
	      sys->av[i][d] += sys->mass[j]*avgpart;
	      sys->av[j][d] -= sys->mass[i]*avgpart;
	    }
	}
    }
}

void AverageInflowVelocity(System *sys)
{

  /*get average inflow velocity*/

  double sumdvdt = 0.0;
  double sumav	= 0.0;
  /* for(int k = 0; k < sys->inflowParticles; k++) */
  /*   { */
  /*     int i = sys->i_inflow[k]; */
  /*     sumdvdt += sys->dvdt[i][1]; */
  /*     sumav += sys->av[i][1]; */
  /*   }	 */

    /* particles within h from the inflow boundary all have same velocity */
  /* for(int k = 0; k <sys->inflowParticles; k++) */
  /*   { */
  /*     int i = sys->i_inflow[k]; */
  /*     sys->dvdt[i][1] = 0.0; */
  /*     sys->av[i][1] = 0.0; */
  /*     sys->dvdt[i][0] = sumdvdt/sys->inflowParticles; */
  /*     sys->av[i][0]+= sumav/sys->inflowParticles; */
  /*     /\* a particle within h is not allowed to flow back*\/ */
  /*     if((sys->Velocity[i][0] <= 0.0) && (sys->dvdt[i][0] <= 0.0)) */
  /* 	{ */
  /* 	  sys->dvdt[i][1] = 0.0; */
  /* 	} */
  /*     if((sys->Velocity[i][0] <= 0.0) && (sys->av[i][0] <= 0.0)) */
  /* 	{ */
  /* 	  sys->av[i][0] = 0.0; */
  /* 	} */
  /*     sys->drhodt[i] = 0.0; /\* combination of depth with p ( & rho) is fixated*\/ */
  /*   } */

}
void derivatives(System *sys, double t)
{
  int ii;
  double xl = 1.0, vl = 1.0;
  double epsilon = 0.5;

  /* after moving particles inflow and outflow */

  /* figure out why this if-statement gives an error when " if{}"*/
  /* if(t > 0.00000001) */
  memset(sys->drhodt, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dvdt[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  memset(sys->av[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);

  if(t>0) {
    inoutflow(sys);
  }


  /* with inflow/outflow boundaries, call inputbp.h every timestep */
  inputbp(sys, xl, vl, t);

  /* track particles within space from inflow boundary */
  sys->inflowParticles = 0;
  for(int i = 0; i < sys->ntotal; i++)
    {
      if(sys->Position[0][i] <= 0.1)
	{
	  sys->i_inflow[sys->inflowParticles] = i;
	  sys->inflowParticles++;
	}
    }


  /* boundary particles also need vxsph */

  for(int i=0;i<sys->NBoundaries;i++) {
    for(int d=0;d<sys->dim;d++)
      sys->Velocity_xsph[sys->ntotal + i][d] = sys->Velocity[sys->ntotal+i][d];
  }

  /* create virtual particles within 2h*/
  ii = sys->ntotal + sys->NBoundaries;

  //  left boundary
  for(int i = 0; i < sys->ntotal; i++) {
    double dx = fabs(sys->Position[i][0]-sys->LeftBoundary);
    if(dx < 2 * sys->hsml[i]) {
      sys->Position[ii][0] = sys->LeftBoundary-sys->Position[i][0];
      sys->Position[ii][1] = sys->Position[i][1];
      sys->Velocity[ii][0] = -sys->Velocity[i][0];
      sys->Velocity[ii][1] = sys->Velocity[i][1];
      sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
      sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
      sys->mass[ii] = sys->mass[i];
      sys->rho[ii] = sys->rho[i];
      sys->Pressure[ii] = sys->Pressure[i];
      sys->Energy[ii] = sys->Energy[i];
      sys->hsml[ii] = sys->hsml[i];
      ii ++;
    }
  }
  //  Right boundary 
  for(int i = 0; i < sys->ntotal; i++) {
    if(fabs(sys->Position[i][0]-sys->RightBoundary) < 2 * sys->hsml[i]) {
      sys->Position[ii][0] = sys->Position[i][0] + 2.*fabs(sys->Position[i][0]-sys->RightBoundary);
      sys->Position[ii][1] = sys->Position[i][1];
      sys->Velocity[ii][0] = -sys->Velocity[i][0];
      sys->Velocity[ii][1] = sys->Velocity[i][1];
      sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
      sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
      sys->mass[ii] = sys->mass[i];
      sys->rho[ii] = sys->rho[i];
      sys->Pressure[ii] = sys->Pressure[i];
      sys->Energy[ii] = sys->Energy[i];
      sys->hsml[ii] = sys->hsml[i];
      ii ++;
    }
  }

  /* for(int k = 0; k < sys->inflowParticles; k++) */
  /*   { */
  /*     int i = sys->i_inflow[k]; */
  /*     for(int j = 0; j <2; j++) */
  /*	{ */
  /*	  sys->Position[ii][0] = -j*0.1 + sys->Position[i][0]; */
  /*	  sys->Position[ii][1] = sys->Position[i][1]; */
  /*	  sys->Position[ii][0] = sys->Position[i][0]; */
  /*	  sys->Velocity[ii][1] = sys->Velocity[i][1]; */
  /*	  sys->Velocity[ii][0] = sys->Velocity[i][0]; */
  /*	  sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0]; */
  /*	  sys->mass[ii] = sys->mass[i]; */
  /*	  sys->Pressure[ii] = sys->Pressure[i]; */
  /*	  sys->Energy[ii] = sys->Energy[i]; */
  /*	  sys->hsml[ii] = sys->hsml[i]; */
  /*	  ii ++; */
  /*	} */
  /*   } */

  //Bottom boundary
  for(int i = 0; i < sys->ntotal; i++)
    {
      if((sys->Position[i][1] <= 2*sys->hsml[i]) && (sys->Position[i][0] <= 12.0))
  	{
  	  sys->Position[ii][0] = sys->Position[i][0];
  	  sys->Position[ii][1] = -sys->Position[i][1];
  	  sys->Velocity[ii][0] = sys->Velocity[i][0];
  	  sys->Velocity[ii][1] = -sys->Velocity[i][1];
  	  sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
  	  sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
  	  sys->mass[ii] = sys->mass[i];
  	  sys->rho[ii] = sys->rho[i];
  	  sys->Pressure[ii] = sys->Pressure[i];
  	  sys->Energy[ii] = sys->Energy[i];
  	  sys->hsml[ii] = sys->hsml[i];
  	  ii ++;
  	}
    }

  // edges 
  for(int i = 0; i < sys->ntotal; i++) {
    double dx = sys->Position[i][0] - sys->LeftBoundary;
    double dy = sys->Position[i][1];

    if((dx < 2*sys->hsml[i])&&(dy < 2*sys->hsml[i]))
      {
  	sys->Position[ii][0] = sys->LeftBoundary-sys->Position[i][0];
  	sys->Position[ii][1] = -sys->Position[i][1];
  	sys->Velocity[ii][0] = -sys->Velocity[i][0];
  	sys->Velocity[ii][1] = -sys->Velocity[i][1];
  	sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
  	sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
  	sys->mass[ii] = sys->mass[i];
  	sys->rho[ii] = sys->rho[i];
  	sys->Pressure[ii] = sys->Pressure[i];
  	sys->Energy[ii] = sys->Energy[i];
  	sys->hsml[ii] = sys->hsml[i];
  	ii ++;
      }

    //  right corner

    dx = fabs(sys->Position[i][0] - sys->RightBoundary);
    if((dx < 2*sys->hsml[i])&&(dy < 2*sys->hsml[i])) {
      sys->Position[ii][0] = sys->RightBoundary+dx;
      sys->Position[ii][1] = -sys->Position[i][1];
      sys->Velocity[ii][0] = -sys->Velocity[i][0];
      sys->Velocity[ii][1] = -sys->Velocity[i][1];
      sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
      sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
      sys->mass[ii] = sys->mass[i];
      sys->rho[ii] = sys->rho[i];
      sys->Pressure[ii] = sys->Pressure[i];
      sys->Energy[ii] = sys->Energy[i];
      sys->hsml[ii] = sys->hsml[i];
      ii ++;
    }
  }

  /* Initialize derivatives: must be after inoutflow.h */
  sys->NumberOfVirtualParticles = ii - sys->ntotal - sys->NBoundaries;

  /* Interaction plus kernel definition*/
  CalculateKernel(sys);

  // up to that it is correct.

  /* solve the density equation */

  DensityEquation(sys);

  /* equation of state: pressure p(rho) */

  EquationOfState(sys);

  /* Momentum equation */
  MomentumEquations(sys);

  /* artificial viscosity for modelling viscosity */
  if(art_visc)
    ArtificialViscosity(sys);

  /* gravity*/
  /* for(int i = 0; i < sys->ntotal; i++) */
  /*   sys->dvdt[i][sys->dim-1] -= 9.81; */




  LeonnardJonesBoundaryForces(sys);

  /*Average velocity : xSPH*/

  AverageVelocity(sys);


  /* AverageInflowVelocity(sys);  */

  /* generate data file with useful info about the simulation*/
  /* if(itime ==1) */
  /*   parameterfile(sys, epsilon, dt); */

}
