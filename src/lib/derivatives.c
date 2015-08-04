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

__inline__ void Kernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  sys->alpha_cubic = 10.0/(7.*M_PI*h*h); /*2D normalization const for cubic spline*/
  sys->alpha_quintic = 7.0/(478.0*M_PI*h*h); /*2D normalization const for quintic spline*/
  sys->alpha_wendland = 7.0/(64.0*M_PI*h*h);
  const double q = r/h;
 
  /*define cubic spline*/
  if(q >= 0.0 && q < 1.0)
    {
      double qp2 = (2-q)*(2-q);
      double q2 = (2-q)*qp2;
      double qp1 = (1-q)*(1-q);
      double q1 = (1-q)*qp1;
      sys->w[niac] = sys->alpha_cubic*(0.25*q2 - q1);
#pragma ivdep
      for(int d = 0; d < sys->dim; d++)
	sys->dwdx[niac][d] = sys->alpha_cubic*(-3 + 2.25*q)*dx[d]/(h*h);
    }
  else if(q >= 1.0 && q < 2.0)
    {
      double qp2 = (2-q)*(2-q);
      double q2 = (2-q)*qp2;
      sys->w[niac] = sys->alpha_cubic*0.25*q2;
      for(int d = 0; d < sys->dim; d++)
	sys->dwdx[niac][d] = -sys->alpha_cubic*0.75*qp2*dx[d]/(h*h*q);
    }

    
     
//   /*define quintic spline*/
//   if(q >= 0.0 && q <= 1.0)
//     { 
//       double qp1 = (1-q)*(1-q)*(1-q)*(1-q);
//       double q1  = (1-q)*qp1;
//       double qp2 = (2-q)*(2-q)*(2-q)*(2-q);
//       double q2  = (2-q)*qp2;
//       double qp3 = (3-q)*(3-q)*(3-q)*(3-q);
//       double q3  = (3-q)*qp3;
//       double qq2 = q*q;
//       double qq3 = q*q*q;
//       sys->w[niac] = sys->alpha_quintic*(q3 - 6*q2 + 15*q1);
// #pragma ivdep
//       for(int d = 0; d < sys->dim; d++)
//       {
// 	sys->dwdx[niac][d] = -sys->alpha_quintic*(50*qq3 - 120*qq2 +120)*dx[d]/(h*h);
//         sys->dwdx_r[niac][d] = sys->alpha_quintic*(240 - 150*q)*dx[d]/(h*h*h*h);
//       }
//     }
//   else if(q >= 1.0 && q <= 2.0)
//     {
//       double qpp2 =(2-q)*(2-q)*(2-q); 
//       double qp2 = (2-q)*qpp2;
//       double q2  = (2-q)*qp2;
//       double qpp3 = (3-q)*(3-q)*(3-q);
//       double qp3  = (3-q)*qpp3;
//       double q3  = (3-q)*qp3;
//   
//       sys->w[niac] = sys->alpha_quintic*(q3 -6*q2);
//       for(int d = 0; d < sys->dim; d++)
//       {
// 	sys->dwdx[niac][d] = -sys->alpha_quintic*(5*qp3 - 30*qp2)*dx[d]/(q*h*h);
//         sys->dwdx_r[niac][d] = ((5*qp3 - 30*qp2)/(q*q*q)+(20*qpp3 - 120*qpp2)/(q*q))*dx[d]/(h*h*h*h);
//       }
//     }
//    else if(q >= 2.0 && q <= 3.0)
//     {
//       double qpp2 =(2-q)*(2-q)*(2-q); 
//       double qp2 = (2-q)*qpp2;
//       double q2  = (2-q)*qp2;
//       double qpp3 = (3-q)*(3-q)*(3-q);
//       double qp3  = (3-q)*qpp3;
//       double q3  = (3-q)*qp3;
//   
//       sys->w[niac] = sys->alpha_quintic*(q3);
//       for(int d = 0; d < sys->dim; d++)
//       {
// 	sys->dwdx[niac][d] = -sys->alpha_quintic*(5*qp3)*dx[d]/(q*h*h);
//         sys->dwdx_r[niac][d] = ((5*qp3)/(q*q*q)+(20*qpp3)/(q*q))*dx[d]/(h*h*h*h);
//       }
//     }  
//   // else sys->dwdx = 0 and sys->w = 0 . They are
//     // initialized to zero by default at the beginning
//   
//     /*define wendland kernel */
//   if(q >= 0.0 && q <= 2.0)
//     {
//       double q3 = (2-q)*(2-q)*(2-q);
//       double q4 = (2-q)*q3;
//       sys->w[niac] = sys->alpha_wendland*(1+2*q)*q4;
// #pragma ivdep
//       for(int d = 0; d < sys->dim; d++)
// 	sys->dwdx[niac][d] = -sys->alpha_wendland*10*q3*dx[d]/(h*h);
//     }
}

void CalculateKernel(System *sys)
{
  double dx[3];
  const int particles = sys->ntotal+sys->NBoundaries+sys->NumberOfVirtualParticles;
  memset(sys->w, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dwdx[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  memset(sys->dwdx_r[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  int niac = 0;
  
  
  for(int i = 0; i < sys->ntotal; i++) // loop over all particles
    {
      int test = (sys->Position[i][0]>= sys->dx);
      for(int j = i+1; j < particles; j++)
	{
	  int check = test&&(j>(sys->ntotal+sys->NBoundaries-1));
	  if(!check) {
	    
	    /*compute distance between 2 particles*/
#pragma ivdep
	    for(int d = 0; d < sys->dim; d++)
	      dx[d] = sys->Position[i][d] - sys->Position[j][d]; 
	    
	    double r2 = dx[0]*dx[0] + dx[1]*dx[1]; /*square the initial separation dist*/
	    double h = 0.5*(sys->hsml[i]+sys->hsml[j]); /* average smoothing length*/
	    if(r2 <= 4.*h*h)
	      {
		sys->i_pair[niac] = i;
		sys->j_pair[niac] = j;
		sys->rij2[niac] = r2;
		/*kernel is piecewise cubic spline*/
		double r = sqrt(r2);
		Kernel(sys, dx, r, h, niac);
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
	  f = sys->w[k]/(0.5*sys->alpha_cubic);
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
      for(int d = 0; d < sys->dim; d++)
      {
	sys->dvdt[i][d] -= sys->mass[j]*(ppart + 0.0*RR*f)*sys->dwdx[k][d];
//       for(int d = 0; (sys->dim-d)&&(j<sys->ntotal); d++)
	sys->dvdt[j][d] += sys->mass[i]*(ppart + 0.0*RR*f)*sys->dwdx[k][d];
      }
    }

#pragma omp parallel for
  for(int i=0;i<sys->ntotal;i++)
    sys->dvdt[i][1] -= g;
}

/*Leonnard- Jones boundary force*/

 void LeonnardJonesBoundaryForces(System *sys)
{
  const double r0 = 0.1; /*initial particle spacing*/
  const double DD = 0.0;
  double dx[2];

  for(int j = sys->ntotal; j < sys->ntotal + sys->NBoundaries; j++)
    {
      for(int i = 0; i < sys->ntotal; i++)
	{
#pragma ivdep	  
	  for(int d = 0; d < sys->dim; d++)
	    dx[d] = sys->Position[i][d] - sys->Position[j][d];

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

void FluctuationComponent(System *sys)
{
    double dv[2];
    double dx[2];
     for(int k = 0; k < sys->NumberOfInteractingParticles; k++) 
     {
      int i = sys->i_pair[k]; 
      int j = sys->j_pair[k];
      for(int d = 0; d  < sys->dim; d++)
      {
          dv[d] = sys->Velocity_xsph[i][d] - sys->Velocity_xsph[j][d];
          dx[d] = sys->Position[i][d] - sys->Position[j][d];
      }
      for(int d = 0; d < sys->dim; d++)
      {
          double vijdwdx = dv[d]*sys->dwdx[k][d];
          sys->fluctuation[i][d] += sys->mass[j]*vijdwdx*dx[d];
          sys->fluctuation[j][d] -= sys->mass[i]*vijdwdx*dx[d]; 
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
            double fluctParts = (sys->fluctuation[i][d]- sys->fluctuation[i][d])*sys->dwdx[k][d];
	    sys->drhodt[i] += sys->mass[j]*(vijdwdx - fluctParts/sys->rho[i]);
            sys->drhodt[j] += sys->mass[i]*(vijdwdx - fluctParts/sys->rho[j]);
	}
    }
}

 void EquationOfState(System *sys)
{
#pragma omp parallel for
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

void DivergenceOfVelocity(System *sys)
{
    double dv[2];
    double dp[2];
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++) // niac = Total number of interacting particles
    {
        int i = sys->i_pair[k]; /* interaction pair : particle i and particle j*/
        int j = sys->j_pair[k];
        for(int d = 0; d < sys->dim; d++)
        {
            dv[d] = sys->Velocity_xsph[i][d] - sys->Velocity_xsph[j][d];
            dp[d] = sys->fluctuation[i][d] - sys->fluctuation[j][d];
        }
        for(int d = 0; d < sys->dim; d++)
	{
            double vijdwdx = dv[d]*sys->dwdx[k][d];
            double    fluc = dp[d]*sys->dwdx[k][d];
            sys->divV[i] -= sys->mass[j]*(vijdwdx-fluc/sys->rho[i])/sys->rho[i];
            sys->divV[j] -= sys->mass[i]*(vijdwdx-fluc/sys->rho[j])/sys->rho[j]; /* double change in sign = OK */
	}
    }
}
void ExtrapolationFromFluidToWallParticles(System *sys)
{
    double dx[2];
    double g = 9.81;
    for(int i = sys->ntotal ; i < sys->ntotal+sys->NBoundaries; i++)
    {
        for(int j = 0; j < sys->ntotal; j++)
        {
            for(int d = 0; d < sys->dim; d++)
                dx[d] = sys->Position[j][d] - sys->Position[i][d];
            double r2 = dx[0]*dx[0] + dx[1]*dx[1];
            double h = 0.5*(sys->hsml[i]+sys->hsml[j]);
            h *= h;
            if(r2 < 4*h)
            {
                sys->numberDensity[i] += sys->w[j];
                /* smoothed velocity of fluid particles within 2h from boundary*/
                for(int d = 0; d < sys->dim; d++)
                {
                    sys->dummyVelocity[i][d] += sys->Velocity[j][d]*sys->w[j]/sys->numberDensity[i];
                }
                /*extrapolate velocity of these fluid particles to adjacent wall particles*/
                sys->Velocity[i][0] = 2*sys->WallVelocity - sys->dummyVelocity[i][0];
                sys->Velocity[i][1] = -sys->dummyVelocity[i][1];

                /*extrapolate the pressure of these fluid particles to adjacent wall particles*/
                double accl = sys->WallAcceleration*dx[0];
                sys->Pressure[i] += (sys->Pressure[j]-sys->rho[j]*g*dx[1]-accl)*sys->w[j]/sys->numberDensity[i];
                
                /*extrapolate the density in using eos*/
                sys->rho[i] = sys->rho0*pow((1+sys->Pressure[i]/sys->CompressionFactor),1.0/sys->AdiabaticConstant);
            }
        }
    }
}
void DeviatoricStressTensor(System *sys)
{
    double dv[2];
    double dp[2];
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        for(int d = 0; d < sys->dim; d++)
        {
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dp[d] = sys->fluctuation[i][d] - sys->fluctuation[j][d];
        }
        sys->gradVxx[i] -= (dv[0]-dp[0]/sys->rho[i])*sys->dwdx[k][0]/sys->rho[i];
        sys->gradVxx[j] -= (dv[0]-dp[0]/sys->rho[j])*sys->dwdx[k][0]/sys->rho[j];
        
        sys->gradVxy[i] -= (dv[0]-dp[0]/sys->rho[i])*sys->dwdx[k][1]/sys->rho[i];
        sys->gradVxy[j] -= (dv[0]-dp[0]/sys->rho[j])*sys->dwdx[k][1]/sys->rho[j];
        
        sys->gradVyx[i] -= (dv[1]-dp[1]/sys->rho[i])*sys->dwdx[k][0]/sys->rho[i];
        sys->gradVyx[j] -= (dv[1]-dp[1]/sys->rho[j])*sys->dwdx[k][0]/sys->rho[j];
        
        sys->gradVyy[i] -= (dv[1]-dp[1]/sys->rho[i])*sys->dwdx[k][1]/sys->rho[i];
        sys->gradVyy[j] -= (dv[1]-dp[1]/sys->rho[j])*sys->dwdx[k][1]/sys->rho[j];
        
        /* define the strain rate tensor by 2 */
        sys->strainRatexx[i] = sys->gradVxx[i] + sys->gradVxx[i];
        sys->strainRatexx[j] = sys->gradVxx[j] + sys->gradVxx[j];
        
        sys->strainRatexy[i] = sys->gradVxy[i] + sys->gradVyx[i];
        sys->strainRatexy[j] = sys->gradVxy[j] + sys->gradVyx[j];
        
        sys->strainRateyx[i] = sys->gradVyx[i] + sys->gradVxy[i];
        sys->strainRateyx[j] = sys->gradVyx[j] + sys->gradVxy[j];
        
        sys->strainRateyy[i] = sys->gradVyy[i] + sys->gradVyy[i];
        sys->strainRateyy[j] = sys->gradVyy[j] + sys->gradVyy[j];
        
        /*deviatoric stress tensor*/
        sys->deviatoricxx[i] = sys->KinVisc[i]*sys->rho[i]*(sys->strainRatexx[i] - 2*sys->divV[i]/3);
        sys->deviatoricxx[j] = sys->KinVisc[j]*sys->rho[j]*(sys->strainRatexx[j] - 2*sys->divV[j]/3);
        
        sys->deviatoricxy[i] = sys->KinVisc[i]*sys->rho[i]*sys->strainRatexy[i];
        sys->deviatoricxy[j] = sys->KinVisc[j]*sys->rho[j]*sys->strainRatexy[j];
        
        sys->deviatoricyx[i] = sys->KinVisc[i]*sys->rho[i]*sys->strainRateyx[i];
        sys->deviatoricyx[j] = sys->KinVisc[j]*sys->rho[j]*sys->strainRateyx[j];
        
        sys->deviatoricyy[i] = sys->KinVisc[i]*sys->rho[i]*(sys->strainRateyy[i] - 2*sys->divV[i]/3);
        sys->deviatoricyy[j] = sys->KinVisc[j]*sys->rho[j]*(sys->strainRateyy[j] - 2*sys->divV[j]/3);
    }
}
void ViscousForce(System *sys)
{
  for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double sums_x  = (sys->deviatoricxx[i]/(sys->rho[i]*sys->rho[i]) + sys->deviatoricxx[j]/(sys->rho[j]*sys->rho[j]))*sys->dwdx[k][0];
               sums_x += (sys->deviatoricxy[i]/(sys->rho[i]*sys->rho[i]) + sys->deviatoricxy[j]/(sys->rho[j]*sys->rho[j]))*sys->dwdx[k][1];
        double sums_y  = (sys->deviatoricyx[i]/(sys->rho[i]*sys->rho[i]) + sys->deviatoricyx[j]/(sys->rho[j]*sys->rho[j]))*sys->dwdx[k][0];
               sums_y += (sys->deviatoricyy[i]/(sys->rho[i]*sys->rho[i]) + sys->deviatoricyy[j]/(sys->rho[j]*sys->rho[j]))*sys->dwdx[k][1];     
        /*x-component of viscous force*/
        sys->dvdt[i][0] += sys->mass[j]*sums_x;
        sys->dvdt[j][0] -= sys->mass[i]*sums_x;
        /*y-component of viscous force*/
        sys->dvdt[i][1] += sys->mass[j]*sums_y;
        sys->dvdt[j][1] -= sys->mass[i]*sums_y;
    }
}

void VorticityFunction(System *sys)
{
    double dv[2];
    double dp[2];
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        for(int d = 0; d < sys->dim; d++)
        {
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dp[d] = sys->fluctuation[i][d] - sys->fluctuation[j][d];
        }
        double classic_sum = dv[0]*sys->dwdx[k][1];
               classic_sum = dv[1]*sys->dwdx[k][0];
        double correction_sum = dp[0]*sys->dwdx[k][1];
               correction_sum = dp[1]*sys->dwdx[k][0];
        sys->vorticity[i] += sys->mass[j]*(classic_sum - correction_sum/sys->rho[i])/sys->rho[i];  
        sys->vorticity[j] += sys->mass[i]*(classic_sum - correction_sum/sys->rho[j])/sys->rho[j];
    }
}
void AparentForces(System *sys)
{
    double dv[2];
    double dx[2];
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        for(int d = 0; d <sys->dim; d++)
        {
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dv[d] = sys->Position[i][d] - sys->Position[j][d];
        }
        for(int d = 0; d < sys->dim; d++)
        {
            double vijdwdx = dv[d]*sys->dwdx[k][d];
            sys->dvdt[i][d] -= 2*sys->mass[j]*vijdwdx*dv[d]/sys->rho[i];
            sys->dvdt[j][d] += 2*sys->mass[i]*vijdwdx*dv[d]/sys->rho[j];
        }
        /* second apparent force m_jvec(r_ij)*( (u_dotgradV)dot(gradw)*/
        double sum_x_i = dv[0]*sys->gradVxx[i] + dv[1]*sys->gradVyx[i];
        double sum_y_i = dv[0]*sys->gradVxy[i] + dv[1]*sys->gradVyy[i];
        double sum_x_j = dv[0]*sys->gradVxx[j] + dv[1]*sys->gradVyx[j];
        double sum_y_j = dv[0]*sys->gradVxy[j] + dv[1]*sys->gradVyy[j];
        double hijdwdx  = (sum_x_i-sum_x_j)*sys->dwdx[k][0];
               hijdwdx += (sum_y_i-sum_y_j)*sys->dwdx[k][1];
        for(int d = 0 ;d < sys->dim; d++)
        {
            sys->dvdt[i][d] -= sys->mass[j]*hijdwdx*dx[d]/sys->rho[i];   
            sys->dvdt[j][d] += sys->mass[i]*hijdwdx*dx[d]/sys->rho[j];
        }
    }
}
        
void AverageInflowVelocity(System *sys)
{

  /*get average inflow velocity*/

  double sumdvdt = 0.0;
  double sumav	= 0.0;
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
  memset(sys->gradVxx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->gradVxy, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->gradVyx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->gradVyy, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->strainRatexx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->strainRatexy, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->strainRateyx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->strainRateyy, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->deviatoricxx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->deviatoricxy, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->deviatoricyx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->deviatoricyy, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->divV, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->numberDensity, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dummyVelocity[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  memset(sys->dummyAcceleration[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  memset(sys->fluctuation[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  memset(sys->vorticity, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->temps[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);

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

//   /* create virtual particles within 2h*/
//   ii = sys->ntotal + sys->NBoundaries;
// 
//   //  left boundary
// #pragma omp parallel for
//   for(int i = 0; i < sys->ntotal; i++) {
//     double dx = fabs(sys->Position[i][0]-sys->LeftBoundary);
//     if(dx < 2 * sys->hsml[i]) {
//       sys->Position[ii][0] = sys->LeftBoundary-dx;
//       sys->Position[ii][1] = sys->Position[i][1];
//       sys->Velocity[ii][0] = -sys->Velocity[i][0];
//       sys->Velocity[ii][1] = sys->Velocity[i][1];
//       sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
//       sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
//       sys->mass[ii] = sys->mass[i];
//       sys->rho[ii] = sys->rho[i];
//       sys->Pressure[ii] = sys->Pressure[i];
//       sys->Energy[ii] = sys->Energy[i];
//       sys->hsml[ii] = sys->hsml[i];
//       ii ++;
//     }
//   }
//   //  Right boundary 
// #pragma omp parallel for
//   for(int i = 0; i < sys->ntotal; i++) {
//     if(fabs(sys->Position[i][0]-sys->RightBoundary) < 2 * sys->hsml[i]) {
//       sys->Position[ii][0] = sys->Position[i][0] + 2.*fabs(sys->Position[i][0]-sys->RightBoundary);
//       sys->Position[ii][1] = sys->Position[i][1];
//       sys->Velocity[ii][0] = -sys->Velocity[i][0];
//       sys->Velocity[ii][1] = sys->Velocity[i][1];
//       sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
//       sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
//       sys->mass[ii] = sys->mass[i];
//       sys->rho[ii] = sys->rho[i];
//       sys->Pressure[ii] = sys->Pressure[i];
//       sys->Energy[ii] = sys->Energy[i];
//       sys->hsml[ii] = sys->hsml[i];
//       ii ++;
//     }
//   }
// 
//   /* for(int k = 0; k < sys->inflowParticles; k++) */
//   /*   { */
//   /*     int i = sys->i_inflow[k]; */
//   /*     for(int j = 0; j <2; j++) */
//   /*	{ */
//   /*	  sys->Position[ii][0] = -j*0.1 + sys->Position[i][0]; */
//   /*	  sys->Position[ii][1] = sys->Position[i][1]; */
//   /*	  sys->Position[ii][0] = sys->Position[i][0]; */
//   /*	  sys->Velocity[ii][1] = sys->Velocity[i][1]; */
//   /*	  sys->Velocity[ii][0] = sys->Velocity[i][0]; */
//   /*	  sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0]; */
//   /*	  sys->mass[ii] = sys->mass[i]; */
//   /*	  sys->Pressure[ii] = sys->Pressure[i]; */
//   /*	  sys->Energy[ii] = sys->Energy[i]; */
//   /*	  sys->hsml[ii] = sys->hsml[i]; */
//   /*	  ii ++; */
//   /*	} */
//   /*   } */
// 
//   //Bottom boundary
// 
// #pragma omp parallel for
//   for(int i = 0; i < sys->ntotal; i++)
//     {
//       if((sys->Position[i][1] <= 2*sys->hsml[i]) && (sys->Position[i][0] <= 12.0))
//   	{
//   	  sys->Position[ii][0] = sys->Position[i][0];
//   	  sys->Position[ii][1] = -sys->Position[i][1];
//   	  sys->Velocity[ii][0] = sys->Velocity[i][0];
//   	  sys->Velocity[ii][1] = -sys->Velocity[i][1];
//   	  sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
//   	  sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
//   	  sys->mass[ii] = sys->mass[i];
//   	  sys->rho[ii] = sys->rho[i];
//   	  sys->Pressure[ii] = sys->Pressure[i];
//   	  sys->Energy[ii] = sys->Energy[i];
//   	  sys->hsml[ii] = sys->hsml[i];
//   	  ii ++;
//   	}
//     }
// 
//   // edges 
// #pragma omp parallel for
//   for(int i = 0; i < sys->ntotal; i++) {
//     double dx = fabs(sys->Position[i][0] - sys->LeftBoundary);
//     double dy = sys->Position[i][1];
// 
//     if((dx < 2*sys->hsml[i])&&(dy < 2*sys->hsml[i]))
//       {
//   	sys->Position[ii][0] = sys->LeftBoundary-dx;
//   	sys->Position[ii][1] = -sys->Position[i][1];
//   	sys->Velocity[ii][0] = -sys->Velocity[i][0];
//   	sys->Velocity[ii][1] = -sys->Velocity[i][1];
//   	sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
//   	sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
//   	sys->mass[ii] = sys->mass[i];
//   	sys->rho[ii] = sys->rho[i];
//   	sys->Pressure[ii] = sys->Pressure[i];
//   	sys->Energy[ii] = sys->Energy[i];
//   	sys->hsml[ii] = sys->hsml[i];
//   	ii ++;
//       }
// 
//     //  right corner
// 
//     dx = fabs(sys->Position[i][0] - sys->RightBoundary);
//     if((dx < 2*sys->hsml[i])&&(dy < 2*sys->hsml[i])) {
//       sys->Position[ii][0] = sys->RightBoundary+dx;
//       sys->Position[ii][1] = -sys->Position[i][1];
//       sys->Velocity[ii][0] = -sys->Velocity[i][0];
//       sys->Velocity[ii][1] = -sys->Velocity[i][1];
//       sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
//       sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
//       sys->mass[ii] = sys->mass[i];
//       sys->rho[ii] = sys->rho[i];
//       sys->Pressure[ii] = sys->Pressure[i];
//       sys->Energy[ii] = sys->Energy[i];
//       sys->hsml[ii] = sys->hsml[i];
//       ii ++;
//     }
//   }
// 
//   /* Initialize derivatives: must be after inoutflow.h */
//   sys->NumberOfVirtualParticles = ii - sys->ntotal - sys->NBoundaries;

  /* Interaction plus kernel definition*/
  CalculateKernel(sys);

  // up to that it is correct.

  FluctuationComponent(sys);
  /* solve the density equation */

  DensityEquation(sys);

  /* equation of state: pressure p(rho) */

  EquationOfState(sys);
  
  /* Momentum equation */
  MomentumEquations(sys);

  /* artificial viscosity for modelling viscosity */
  if(art_visc)
    ArtificialViscosity(sys);
  
 /* apparent forces*/
 AparentForces(sys);
 
  
 
 VorticityFunction(sys);

  
 /* implement the free-slip and no-slip boundary
  conditions here. For free slip, omit
  the extrapolated quanties from the viscous 
  forces. Whereas for no-slip include the extrapolated
  values of wall particles*/ 
 /*___________________________________________________________*/
 ExtrapolationFromFluidToWallParticles(sys); 

 DivergenceOfVelocity(sys);
 
 DeviatoricStressTensor(sys);
 
 ViscousForce(sys);
/*_____________________________________________________________*/

 /* this is a repulsive force for implenting the  
    no pernetration boundary conditions*/
/*_____________________________________________________________ */
 LeonnardJonesBoundaryForces(sys);
/*_____________________________________________________________*/
  /*Average velocity : xSPH*/

  AverageVelocity(sys);
/*_____________________________________________________________*/

  /* AverageInflowVelocity(sys);  */

  /* generate data file with useful info about the simulation*/
  /* if(itime ==1) */
  /*   parameterfile(sys, epsilon, dt); */

}
