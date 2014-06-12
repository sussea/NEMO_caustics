/*
 * This file provides a caustic ring halo potential for use in stellar models. To add the
 * caustic ring potential to your own potential, simply '#include' this file and add a
 * call to apply_caustics_pot in your potential file. The apply_caustics_pot function
 * takes in your position and acceleration in x y z as well as your potential, and adds
 * the contribution due to caustic rings to your acceleration and potential.
 * An example of how to use this code should be provided in the "mpc.c" file, and more
 * details may be found in the README. A stand-alone caustics potential can be found in
 * caustics_standalone.c
 *
 * Refs: Tam, H. 2012, ArXiv e-prints
 *       Duffy L. D., & Sikivie, P. 2008, Phys. Rev. D, 61, 063508
 * 
 * htam wrote original code to calculate gfield_close (uses functions f1-f5, T1-T4)
 * jdumas combined gfield_close and gfield_far to calculate total caustic ring halo acceleration
 * 
 * 
 * Version 1
 * 30-may-14
 * jdumas added flags on minimum values for rho and rpar
 * adam susser cleaned bloat in gfield_close()
 */

#include <math.h>
#include <complex.h>
#include <stdio.h>

//local double G = 1.0;
const double G_caustics = 1.0;
const double a_n[] = {1.0,40.1,20.1,13.6,10.4,8.4,7.0,6.1,5.3,4.8,4.3,4.0,3.7,3.4,3.2,3.0,2.8,2.7,2.5,2.4,2.3};
const double V_n[] = {1.0,517,523,523,523,522,521,521,520,517,515,512,510,507,505,503,501,499,497,496,494};
const double rate_n[] = {1.0,53,23,14,10,7.8,6.3,5.3,4.5,3.9,3.4,3.1,2.8,2.5,2.3,2.1,2.0,1.8,1.7,1.6,1.5};
const double p_n[] = {1.0,0.3,0.3,1.0,0.3,0.15,0.12,0.6,0.23,0.41,0.25,0.19,0.17,0.11,0.09,0.09,0.09,0.09,0.09,0.09,0.09};

const double ONE_THIRD=1.0/3.0;
#define c1 (double complex) 1.0
#define c2 (double complex) 2.0
#define c3 (double complex) 3.0
#define c4 (double complex) 4.0
#define c12 (double complex) 12.0
#define c05 (double complex) 0.5
#define II (double complex) 1.0*I
const double TOLERANCE_caustics=0.00001;
const double TOLERANCE2_caustics=0.001;

#define X_caustics 0
#define Y_caustics 1
#define Z_caustics 2

// now let's compute the values of T at which the various poles cross the real axis
inline double complex f1(double complex x) {
  return (2.0 + 4.0*x) / 3.0;
}

inline double complex f2(double complex x, double complex z) {
  return 4.0 * ((1.0-x) * (1.0-x) - 3.0 * z*z);
} 

inline double complex f3(double complex x, double complex z) {
  return -128.0 * cpow(x-1.0,3.0) - 1728.0 * z*z - 1152.0 * (x-1.0) * z*z;
}

// for convenience define f4, which appears frequently in the roots
double complex f4(double complex x, double complex z) {
  return cpow( (f3(x,z) + csqrt( -256.0*f2(x,z)*f2(x,z)*f2(x,z) + f3(x,z)*f3(x,z) )) / 2.0, ONE_THIRD );
}

// f6 my suspicion is that the gravitational field in r and z constitute the real and complex components of the individual terms in f5; so let's do a test to see if this s right
inline double complex f5(double complex x, double complex z, double complex T) {
  return (csqrt(c2*T - c1 + x - II*z)) / c2;
}


inline double complex T1(double complex x, double complex z) {
  if (creal(x) >= -0.125 && cabs(f4(x,z)) < TOLERANCE2_caustics) {
    double complex temp = cpow((-1.0 - 6.0*x + 15.0*x*x -8.0*x*x*x), ONE_THIRD) / 3.0;
    return c05 * (c1 - csqrt(f1(x) / 2.0 + temp) - csqrt( f1(x) - temp - c2*x/csqrt(f1(x) / 2.0 + temp) ));
  } else {
    double complex temp = csqrt(f1(x) / 2.0 + f2(x,z) / (c3 * f4(x,z)) + f4(x,z) / c12);
    return c05 * (c1 - temp - csqrt( f1(x) - f2(x,z) / (c3 * f4(x,z)) - f4(x,z) / c12 - c2 * x / temp ));
  }
}

inline double complex T2(double complex x, double complex z) {
  if (creal(x) >= -0.125 && cabs(f4(x,z)) < TOLERANCE2_caustics) {
    double complex temp = cpow(( -1.0 - 6.0 * x + 15.0 * x*x - 8.0 * x*x*x), ONE_THIRD) / 3.0;
    return c05 * (c1 - csqrt(f1(x) / 2.0 + temp) + csqrt( f1(x) - temp - c2 * x / sqrt(f1(x) / 2.0 + temp) ));
  } else {
    double complex temp = csqrt(f1(x) / 2.0 + f2(x,z) / (c3 * f4(x,z)) + f4(x,z) / c12);
    return c05 * (c1 - temp + csqrt( f1(x) - f2(x,z) / (c3 * f4(x,z)) - f4(x,z) / c12 - c2 * x / temp ));
  }
}

inline double complex T3(double complex x, double complex z) {
  if (creal(x) >= -0.125 && cabs(f4(x,z)) < TOLERANCE2_caustics) {
    double complex temp = cpow((-1.0 - 6.0 * x + 15.0 * x*x - 8.0 * x*x*x), ONE_THIRD)/3.0;
    return c05 * (c1 + csqrt(f1(x) / 2.0 + temp) - csqrt( f1(x) - temp + c2 * x / csqrt(f1(x) / 2.0 + temp) ));
  } else {
    double complex temp = csqrt(f1(x) / 2.0 + f2(x,z) / (c3 * f4(x,z)) + f4(x,z) / c12);
    return c05*(c1 + temp - csqrt( f1(x) - f2(x,z) / (c3 * f4(x,z)) - f4(x,z) / c12 + c2 * x / temp ));
  }
}

inline double complex T4(double complex x, double complex z) {
  if (creal(x) >= -0.125 && cabs(f4(x,z))<TOLERANCE2_caustics) {
    double complex temp = cpow((-1.0 - 6.0 * x + 15.0 * x*x - 8.0 * x*x*x), ONE_THIRD)/3.0;
    return c05 * (c1 + csqrt(f1(x) / 2.0 + temp) + csqrt( f1(x) - temp + c2 * x / csqrt(f1(x) / 2.0 + temp) ));
  } else {
    double complex temp = csqrt(f1(x) / 2.0 + f2(x,z) / (c3 * f4(x,z)) + f4(x,z) / c12);
    return c05 * (c1 + temp + csqrt(f1(x) - f2(x,z) / (c3*f4(x,z)) - f4(x,z) / c12 + c2 * x / temp ));
  }
}

// calculate gfield close to the nth caustic ring flow
// two cases: within the caustic, all four roots are real, so we use T1 for d4 and T2,3,4 for d3
// outside the caustic there should be two real roots and two complex roots
// the vector intlimits hold the limits of integration; for the first case, intlimits is of length four; second, length two
void gfield_close(double rho, double z, int n, double *rfield, double *zfield) {


  z = (double complex)(z / p_n[n]);

  double complex x = (double complex)(( rho - a_n[n] ) / p_n[n]);

  double im[4] = {cabs(cimag(T1(x,z))), cabs(cimag(T2(x,z))), cabs(cimag(T3(x,z))), cabs(cimag(T4(x,z)))};
  double re[4] = {creal(T1(x,z))      , creal(T2(x,z))      , creal(T3(x,z))      , creal(T4(x,z))};

  int count = 0;
  int i, j;
  double temp;

  // get only the roots we want
  for (i = 0; i < 4; ++i) {
    if (im[i] < TOLERANCE_caustics) {
      re[count] = re[i];    
      ++count; // count the roots - there should be 2 or 4
    }
  }

  // sort all 2 or 4 roots in ascending order
  for (i = 0; i < count - 1; ++i) {
    for (j = i + 1; j < count; ++j)
      if (re[i] > re[j]) {
        temp  = re[i];
        re[i] = re[j];
        re[j] = temp;
      }
  }

  // roots stored in re[]
  double factor = -8.0 * M_PI * G_caustics * rate_n[n] * 4498.6589 / (rho * V_n[n] * 1.0226831);

  *rfield += factor * ( creal(f5(x,z,re[1])) - creal(f5(x,z,re[0])) - 0.5 ); // it seems that the answer is just off by 0.5, so we subtract it by hand
  *zfield += factor * ( cimag(f5(x,z,re[1])) - cimag(f5(x,z,re[0])) );

  if (count == 4){
    *rfield += factor * ( creal(f5(x,z,re[3])) - creal(f5(x,z,re[2])) );
    *zfield += factor * ( cimag(f5(x,z,re[3])) - cimag(f5(x,z,re[2])) );	
  }

}

// calculate gfield far away from the nth caustic ring flow
void gfield_far(double rho, double z, int n, double *rfield, double *zfield) {

  // simulation units are kpc, gyr, ms=222288.47*Ms
  double r_squared = rho*rho + z*z;
  double shift = a_n[n] + (p_n[n] / 4.0);
  double A_n = (8.0 * M_PI * G_caustics * rate_n[n] * 4498.6589) / (V_n[n] * 1.0226831);  //convert from M_sun/yr*sr to m_sim/gyr*sr with 4498 and km/s to kpc/gyr with 1.0226

  // caustic radius shifted by 0.25 so that g goes to zero beyond a_n[n]
  double s = hypot(r_squared - shift*shift, 2.0 * shift * z);
  double factor = -A_n / ( s * (2.0 * shift*shift + s) );

  // result in kpc/gyr^2
  *rfield += factor * (r_squared * rho - shift * shift * rho);
  *zfield += factor * (r_squared * z   + shift * shift * z);
}

// call this from your own potential file (double version)
void apply_caustic_pot_double(double *pos, double *acc, double *pot) {
  double rho, z, rfield, zfield;
  double r, l, tr, tl; // r (right), l (left), tr (top right), tl (top left)
  int n;

  rfield = 0.0;
  zfield = 0.0;
  rho = sqrt(pos[X_caustics]*pos[X_caustics] + pos[Y_caustics]*pos[Y_caustics]);
  z = pos[Z_caustics];

  // rho cannot be zero (causes nan in near field and ax and ay at origin)
  if (rho < 0.000001) {
    rho = 0.000001;
  }

  // calculate gfield at a position (x,y,z) by adding the contributions from all n caustic ring flows
  for (n = 1; n <= 20; ++n) {

    //the caustic ring has a tricusp boundary (see Tam 2012)
    r  = (3.0 - sqrt( 1.0 + (8.0 / p_n[n]) * (rho - a_n[n]) )) / 4.0;
    l  = (3.0 + sqrt( 1.0 + (8.0 / p_n[n]) * (rho - a_n[n]) )) / 4.0;
    tr = 2.0 * p_n[n] * sqrt(pow(r, 3.0) * (1.0 - r));
    tl = 2.0 * p_n[n] * sqrt(pow(l, 3.0) * (1.0 - l));

    //if position (rho,z) is inside ring use gfield_close, else use gfield_far
    if ( (z <= tr && z >= 0.0 && rho >= a_n[n] && rho <= a_n[n] + p_n[n])
     || (z >= tl  && z <= tr  && rho >= (a_n[n] - p_n[n] / 8.0) && rho <= a_n[n])
     || (z >= -tr && z <= 0.0 && rho >= a_n[n] && rho <= a_n[n] + p_n[n])
     || (z <= -tl && z >= -tr && rho >= (a_n[n] - p_n[n] / 8.0) && rho <= a_n[n]) )
    {
      gfield_close(rho, z, n, &rfield, &zfield);
    } else {
      gfield_far(rho, z, n, &rfield, &zfield);        
    }
  }

  // calculate potential at a position (x,y,z) by adding the contributions from all n caustic ring flows
  for (n = 1; n <= 20; ++n) {
    double r_squared = rho*rho + z*z;
    double shift = a_n[n] + p_n[n] / 4.0;  //caustic radius shifted by 0.25 so that g goes to zero beyond a_n[n]
    double A_n = (8.0 * M_PI * G_caustics * rate_n[n] * 4498.6589) / (V_n[n] * 1.0226831);  //convert from M_sun/yr*sr to m_sim/gyr*sr with 4498 and km/s to kpc/gyr with 1.0226
    double s = hypot(r_squared - shift*shift, 2.0 * shift *z);
    *pot += (A_n / 2.0) * log(1.0 + ( s / (2.0 * a_n[n]*a_n[n]) ));
  }

  acc[X_caustics] += rfield * pos[X_caustics] / rho;
  acc[Y_caustics] += rfield * pos[Y_caustics] / rho;
  acc[Z_caustics] += zfield;

}
// call this from your own potential file (float version)
void apply_caustic_pot_float(float *pos, float *acc, float *pot) {
  double pos_d[3];
  double acc_d[3];
  double pot_d = 0;

  pos_d[X_caustics] = (double) pos[X_caustics]; pos_d[Y_caustics] = (double) pos[Y_caustics]; pos_d[Z_caustics] = (double) pos[Z_caustics];
  acc_d[X_caustics] = 0;                        acc_d[Y_caustics] = 0;                        acc_d[Z_caustics] = 0;

  apply_caustic_pot_double(pos_d,acc_d,&pot_d);

  acc[X_caustics] += (float) acc_d[X_caustics];
  acc[Y_caustics] += (float) acc_d[Y_caustics];
  acc[Z_caustics] += (float) acc_d[Z_caustics];
  *pot += (float) pot_d;

}

