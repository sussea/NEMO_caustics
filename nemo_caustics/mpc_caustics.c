/*
 * mpc.c: procedures for intializing and calculating the acceleration and
 *             	potential of a combination of 3 potentials:
 *			Miyamoto Nagai Disk
 *			Plummer Sphere
 *			Caustic Ring Halo
 *
 * This file was adapted from mpl4.c created by bwillett. Below is the file history.
 *
 *      Refs: BT pp. 43-44; Miyamoto and Nagai PASJ 27, 533 (1975)
 *	Miyamoto Potential Phi_m (R, z) = -GMdisk / sqrt (R^2+(a+sqrt(z^2+b^2))^2)
 * 	        Parameters: a, b (shape parameters), M (mass); G=1
 *              Names used: miya_ascal, miya_bscal, miya_mass
 *	Plummer Potential Phi_p (r) = -GMsph / (r + rc) = -GMsph/ (sqrt(R^2 + z^2) + rc)
 *		Parameters: rc, Msph; G=1
 *		Names used: plu_rc, plu_mass
 * 	Logarithmic Halo Phi_l (R, z) = vhalo^2 ln(R^2 + (z^2/q^2) + d^2)
 *		Parameters: vhalo, q, d
 *		Names used: vhalo, q, d
 *
 *  March 90 Stefano Casertano, University of Pittsburgh
 * 10-Nov-90 inserted omega as first parameter for standard Nemo  PJT
 *  6-oct-91 fixed bug - and made code accept both XYZ and XZY versions (pjt)
 *  7-mar-92 merged sun and 3b1 versions once more			 pjt
 *    oct-93 get_pattern
 *  12-mar-07 bwillett changed from miyamoto.c to mpl.c
 *  27-apr-07 bwillett modifying mpl3.c to change the acceleration fields
 *  1-may-07 bwillett created mpl4.c - took out gravitational constant
 *	all masses scaled as 1 M.U. = 222288.47 Ms 
 *	length unit: 1 kpc time unit: 1 Gyr
 *
 *
 * The Caustic Ring Halo acceleration is calculated using functions f1-f5, T1-T4, gfield_far, and gfield_close

 * Refs: Tam, H. 2012, ArXiv e-prints
 *       Duffy L. D., & Sikivie, P. 2008, Phys. Rev. D, 61, 063508
 * 
 * The logarithmic halo calculation (variables vhalo, q, d, lpar, rcyl) is NOT used in this program
 */

#include <stdinc.h>
#include <potential_float.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>

local double omega = 0.0;		/* pattern speed */
local double miya_ascal = 0.0;
local double miya_bscal = 1.0;
local double miya_mass  = 1.0;
local double plu_rc   = 1.0;
local double plu_mass = 1.0;
local double vhalo = 1.0;
local double q = 1.0;
local double d = 1.0;
local const double G = 1.0;

local const int X = 0;
local const int Y = 1;
local const int Z = 2;

void inipotential (int *npar, double *par, string name) {
  int n;
  n = *npar;
  if (n>0) omega = par[0];
  if (n>1) miya_ascal = par[1];
  if (n>2) miya_bscal = par[2];
  if (n>3) miya_mass  = par[3];
  if (n>4) plu_rc   = par[4];
  if (n>5) plu_mass = par[5];
  if (n>6) vhalo = par[6];
  if (n>7) q = par[7];
  if (n>8) d = par[8];
  if (n>9) warning("mpc: only first 9 parameters recognized");
}

// galaxy disk calculations (acc[] and pot)
void apply_miyamoto_pot(double *pos, double *acc, double *pot) {
  double qpar, apar, spar;
  qpar = hypot(pos[Z], miya_bscal);
  apar = miya_ascal + qpar;
  spar = pos[X]*pos[X] + pos[Y]*pos[Y] + (miya_ascal + qpar)*(miya_ascal + qpar);

  *pot -= miya_mass / sqrt(spar);

  acc[X] -= miya_mass * pos[X] / pow(spar, 1.5);
  acc[Y] -= miya_mass * pos[Y] / pow(spar, 1.5);
  acc[Z] -= miya_mass * pos[Z] * apar / (qpar * pow(spar, 1.5));
}

// galaxy bulge calculations (acc[] and pot)
void apply_plummer_pot(double *pos, double *acc, double *pot) {
  double ppar, rpar;
  ppar = sqrt(pos[X]*pos[X] + pos[Y]*pos[Y] + pos[Z]*pos[Z]) + plu_rc;
  rpar = ppar - plu_rc;
  if (rpar < 0.000001) // mke sure rpar != zero (causes nan in accx accy accz at origin)
    rpar = 0.000001;

  *pot -= plu_mass / ppar;

  acc[X] -= plu_mass * pos[X] / (rpar * ppar * ppar);
  acc[Y] -= plu_mass * pos[Y] / (rpar * ppar * ppar);
  acc[Z] -= plu_mass * pos[Z] / (rpar * ppar * ppar);
}


void potential_double(int *ndim, double *pos, double *acc, double *pot, double *time) {

  *pot = 0;
  acc[X] = 0;
  acc[Y] = 0;
  acc[Z] = 0;
  apply_miyamoto_pot(pos, acc, pot);
  apply_plummer_pot(pos, acc, pot);
  apply_caustic_pot_double(pos, acc, pot);
}
/*
 * Unused stuff:
 * log halo calculations (NOT used in this program)
 * rcyl = hypot(pos[X],pos[Y]);
 * lpar = (rcyl*rcyl) + ((pos[Z]/q)*(pos[Z]/q)) + (d*d);
 */
