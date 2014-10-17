/*
***********************************************************************
* subprogram to write useful info about the simulation
***********************************************************************
*/
#include <stdio.h>
#include "param.h"
#include "System.h"

void parameterfile(System *sys, double c,
double B, double g, double alpha, double epsilon, double dt)
{
	FILE *fip;
	char filename[30];
	int file_no = 1;
	sprintf(filename, "simulation%d.txt", file_no);
	fip = fopen("filename", "w+");
	fputs(" **************************************** \n", fip);
	fputs("Smoothed Particle Hydrodynamics 2D Model\n", fip);
	fputs("Boundary with leonnard-Jones force & virtual particles\n", fip);
	fputs("info about the simulation:\n", fip);
	fputs(" **************************************** \n", fip);
	fprintf(fip, "number of fluid particles: %d\n", sys->ntotal);
	fprintf(fip, "number of boundary particles: %d\n", sys->nvirt);
	fprintf(fip, "timestep: %lf [s]\n", dt);
	fputs("Processes included:\n", fip);
	if(art_visc)
	{
		fprintf(fip, "artificial viscosity: %d\n", art_visc);
		fprintf(fip, "alpha: %lf\n", alpha);
	}
	else
	{
		fprintf(fip, "artificial viscosity: %d\n", art_visc);
	}	
	if(avg_vel)
	{
		fprintf(fip, "Average velocity: %d\n", avg_vel);
		fprintf(fip, "epsilon: %lf\n", epsilon);
	}
	else
	{
		fprintf(fip, "Average velocity: %d\n", avg_vel);
	}
	fprintf(fip, "Prevent tensile instability: %d\n", virt_pres);
	fputs("Particle info:\n", fip);
	fprintf(fip, "fluid article Smoothing length: %lf [m]\n", sys->hsml[1]);
	fprintf(fip, "Boundary particle smoothing length: %lf [m]\n", sys->hsml[1+sys->ntotal]);
	fprintf(fip, "mass of fluid particle: %lf [kg/m]\n", sys->mass[1]);
	fprintf(fip, "mass of boundary particle: %lf [kg/m]\n", sys->mass[1+sys->ntotal]);
	fputs("parameters:\n", fip);
	fprintf(fip, "gravitational acceleration: %lf [m/s^2]\n", g);
	fprintf(fip, "sound speed: %lf [m/2]\n", c);
	fprintf(fip, "Compression factor B = rho*c^2/7 = %lf [kg/m^2s^2]\n", B);
	fputs("Written by Kalale Chola\n", fip);
	fputs(" **************************************** \n", fip);
}	
	
	
	
	
	
	
