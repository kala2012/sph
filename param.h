/*
**************************************************
* including file for parameters and constants used
* in all "sph.c" files.
**************************************************
*/
#ifndef param_h
#define _param_h_

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define pi 3.1415926535897932385
#define maxn 12000
#define max_int  12000

int itype[maxn];
//int ntotal; /* total number of fluid particles */
//int nvirt;  /* number of virtual particles */
//int nvll;   /* number of boundary particles */
/* double **x[dim][maxn]; /\* initialize the position vector *\/ */
/* double **v[dim][maxn]; /\* initialize the velocity vector *\/ */
/* double **vxsph[dim][maxn]; /\* move particles with this average velocity *\/ */
/* double *mass[maxn]; /\* initialize particle mass *\/ */
/* double *p[maxn]; /\* initialize particle pressure *\/ */
/* double *u[maxn]; /\* initialize internal energy *\/ */
/* double *hsml[maxn]; /\* initialize smoothing length h *\/ */
/* double *rho[maxn]; /\* initialize particle density *\/ */
//bool virt_pres, art_visc, avg_vel, input_file, vp_file;

bool virt_pres; /* virtual pressure to prevent tensile instability */
bool art_visc; /* artificial viscosity for viscosity simul */
bool avg_vel; /* average velocity for xSPH variant */
bool input_file; /* input fluid particles from file */
bool vp_file ; /* input boundary particles from file */

#endif
