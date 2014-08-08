/*
Copyright (C) 2014 Julie Dumas
Copyright (C) 2014 Jake Weiss

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "nbody_caustic.h"

const real G_caustics = 1.0;
const real a_n[] = {1.0,40.1,20.1,13.6,10.4,8.4,7.0,6.1,5.3,4.8,4.3,4.0,3.7,3.4,3.2,3.0,2.8,2.7,2.5,2.4,2.3};
const real V_n[] = {1.0,517,523,523,523,522,521,521,520,517,515,512,510,507,505,503,501,499,497,496,494};
const real rate_n[] = {1.0,53,23,14,10,7.8,6.3,5.3,4.5,3.9,3.4,3.1,2.8,2.5,2.3,2.1,2.0,1.8,1.7,1.6,1.5};
const real p_n[] = {1.0,0.3,0.3,1.0,0.3,0.15,0.12,0.6,0.23,0.41,0.25,0.19,0.17,0.11,0.09,0.09,0.09,0.09,0.09,0.09,0.09};
real delta_V_max = 0;
real avg_delta_V_sum2 = 0;

#define ONE_THIRD ((double)1.0/(double)3.0)
#define c1 (double complex) 1.0
#define c2 (double complex) 2.0
#define c3 (double complex) 3.0
#define c4 (double complex) 4.0
#define c12 (double complex) 12.0
#define c05 (double complex) 0.5
#define II (double complex) 1.0*I
#define TOLERANCE_caustics ((double)0.00001)
#define TOLERANCE2_caustics ((double)0.001)

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
    double complex temp = cpow((-1.0 - 6.0*x + 15.0*x*x - 8.0*x*x*x), ONE_THIRD)/3.0;
    return c05 * (c1 - csqrt(f1(x) / 2.0 + temp) - csqrt( f1(x) - temp - c2*x / csqrt(f1(x)/2.0 + temp) ));
  } else {
    double complex temp = csqrt(f1(x) / 2.0 + f2(x,z) / (c3 * f4(x,z)) + f4(x,z) / c12);
    return c05 * (c1 - temp - csqrt( f1(x) - f2(x,z) / (c3 * f4(x,z)) - f4(x,z) / c12 - c2 * x / temp ));
  }
}

inline double complex T2(double complex x, double complex z) {
  if (creal(x) >= -0.125 && cabs(f4(x,z)) < TOLERANCE2_caustics) {
    double complex temp = cpow((-1.0 - 6.0*x + 15.0*x*x - 8.0*x*x*x), ONE_THIRD)/3.0;
    return c05 * (c1 - csqrt(f1(x) / 2.0 + temp) + csqrt( f1(x) - temp - c2*x / csqrt(f1(x)/2.0 + temp) ));
  } else {
    double complex temp = csqrt(f1(x) / 2.0 + f2(x,z) / (c3 * f4(x,z)) + f4(x,z) / c12);
    return c05 * (c1 - temp + csqrt( f1(x) - f2(x,z) / (c3 * f4(x,z)) - f4(x,z) / c12 - c2 * x / temp ));
  }
}

inline double complex T3(double complex x, double complex z) {
  if (creal(x) >= -0.125 && cabs(f4(x,z)) < TOLERANCE2_caustics) {
    double complex temp = cpow((-1.0 - 6.0*x + 15.0*x*x - 8.0*x*x*x), ONE_THIRD)/3.0;
    return c05 * (c1 + csqrt(f1(x) / 2.0 + temp) - csqrt( f1(x) - temp + c2*x / csqrt(f1(x)/2.0 + temp) ));
  } else {
    double complex temp = csqrt(f1(x) / 2.0 + f2(x,z) / (c3 * f4(x,z)) + f4(x,z) / c12);
    return c05 * (c1 + temp - csqrt( f1(x) - f2(x,z) / (c3 * f4(x,z)) - f4(x,z) / c12 + c2 * x / temp ));
  }
}

inline double complex T4(double complex x, double complex z) {
  if (creal(x) >= -0.125 && cabs(f4(x,z))<TOLERANCE2_caustics) {
    double complex temp = cpow((-1.0 - 6.0*x + 15.0*x*x - 8.0*x*x*x), ONE_THIRD)/3.0;
    return c05 * (c1 + csqrt(f1(x) / 2.0 + temp) + csqrt( f1(x) - temp + c2*x / csqrt(f1(x)/2.0 + temp) ));
  } else {
    double complex temp = csqrt(f1(x) / 2.0 + f2(x,z) / (c3 * f4(x,z)) + f4(x,z) / c12);
    return c05 * (c1 + temp + csqrt(f1(x) - f2(x,z) / (c3*f4(x,z)) - f4(x,z) / c12 + c2 * x / temp ));
  }
}

// returns number of useful values in T (number of flows)
int getTVals(real rho, real z, int n, real* T){

  double complex X = (double complex)(( rho - a_n[n] ) / p_n[n]);
  double complex Z = (double complex)(z / p_n[n]);

  real im[4] = {cabs(cimag(T1(X,Z))), cabs(cimag(T2(X,Z))), cabs(cimag(T3(X,Z))), cabs(cimag(T4(X,Z)))};
  real re[4] = {creal(T1(X,Z))      , creal(T2(X,Z))      , creal(T3(X,Z))      , creal(T4(X,Z))};

  int count = 0;
  int i, j;
  real temp;
  // get only the roots we want
  for (i = 0; i < 4; ++i) {
    if (im[i] < TOLERANCE_caustics) {
      T[count] = re[i];
      ++count; // count the roots - there should be 2 or 4
    }
  }

  // sort all 2 or 4 roots in ascending order
  for (i = 0; i < count - 1; ++i) {
    for (j = i + 1; j < count; ++j)
      if (T[i] > T[j]) {
        temp  = T[i];
        T[i] = T[j];
        T[j] = temp;
      }
  }
  return count;
}

// calculate gfield 	 to the nth caustic ring flow
// two cases: within the caustic, all four roots are real, so we use T1 for d4 and T2,3,4 for d3
// outside the caustic there should be two real roots and two complex roots
void gfield_close(real rho, real z, int n, real *rfield, real *zfield) {

  // roots stored in T[]
  real T[4] = {0, 0, 0, 0};
  int T_count = getTVals(rho, z, n, T);

  real factor = -8.0 * M_PI * G_caustics * rate_n[n] * 4498.6589 / (rho * V_n[n] * 1.0226831);

  double complex X = (double complex)(( rho - a_n[n] ) / p_n[n]);
  double complex Z = (double complex)(z / p_n[n]);

  *rfield += factor * ( creal(f5(X,Z,T[1])) - creal(f5(X,Z,T[0])) - 0.5 ); // it seems that the answer is just off by 0.5, so we subtract it by hand
  *zfield += factor * ( cimag(f5(X,Z,T[1])) - cimag(f5(X,Z,T[0])) );

  if (T_count == 4){
    *rfield += factor * ( creal(f5(X,Z,T[3])) - creal(f5(X,Z,T[2])) );
    *zfield += factor * ( cimag(f5(X,Z,T[3])) - cimag(f5(X,Z,T[2])) );	
  }
}

// calculate gfield far away from the nth caustic ring flow
void gfield_far(real rho, real z, int n, real *rfield, real *zfield) {

  // simulation units are kpc, gyr, ms=222288.47*Ms
  real r_squared = rho*rho + z*z;
  real shift = a_n[n] + (p_n[n] / 4.0);
  real A_n = (8.0 * M_PI * G_caustics * rate_n[n] * 4498.6589) / (V_n[n] * 1.0226831);  //convert from M_sun/yr*sr to m_sim/gyr*sr with 4498 and km/s to kpc/gyr with 1.0226

  // caustic radius shifted by 0.25 so that g goes to zero beyond a_n[n]
  real s = hypot(r_squared - shift*shift, 2.0 * shift * z);
  real factor = -A_n / ( s * (2.0 * shift*shift + s) );

  // result in kpc/gyr^2
  *rfield += factor * (r_squared * rho - shift * shift * rho);
  *zfield += factor * (r_squared * z   + shift * shift * z);
}

real get_density_close(real rho, real z, int n) {
  real density = 0;
  real T[4];
  int T_count = getTVals(rho,z,n,T);
  int flow;
  for (flow = 0; flow < T_count; ++flow) {
    real alpha_factor; // alpha_factor is (s*alpha^2)/2
    if (fabs(T[flow]) > 0.0000001) {
      alpha_factor = z*z/(T[flow]*T[flow])/(4*p_n[n]);
    } else {
      alpha_factor = a_n[n] - rho + p_n[n]*(T[flow]-1)*(T[flow]-1); // if negative, alpha*tao_0 is imaginary so we overcounted the flows/poles, ignore this one
      if ((fabs(alpha_factor) > 0.00001) && (alpha_factor < 0)) {
        continue;
      }
    }
    real D = 2*(V_n[n]*1.0226831)*(p_n[n]*T[flow]*(T[flow]-1) + alpha_factor);//convert from km/s to kpc/gyr with 1.0226
    density += fabs(4498.6589*rate_n[n]/(rho*D)); //convert from M_sun/yr*sr to m_sim/gyr*sr with 4498
  }
  return density;
}

// return 1 if inside the tricusp boundary/caustic ring envelope (see Tam 2012)
int in_caustic_envelope(mwvector pos, int n){

  real rho = hypot(X(pos),Y(pos));

  // r (right), l (left), tr (top right), tl (top left)
  real R=(3.0-mw_sqrt(1.0+(8.0/p_n[n])*(rho-a_n[n])))/4.0;
  real l=(3.0+mw_sqrt(1.0+(8.0/p_n[n])*(rho-a_n[n])))/4.0;
  real tr=2.0*p_n[n]*mw_sqrt(cube(R)*(1.0-R));
  real tl=2.0*p_n[n]*mw_sqrt(cube(l)*(1.0-l));


  return( (Z(pos)<=tr && Z(pos)>=0.0 && rho>=a_n[n] && rho<=a_n[n]+p_n[n])
   || (Z(pos)>=tl && Z(pos)<=tr && rho>=(a_n[n]-p_n[n]/8.0) && rho<=a_n[n])
   || (Z(pos)>=-tr && Z(pos)<=0.0 && rho>=a_n[n] && rho<=a_n[n]+p_n[n])
   || (Z(pos)<=-tl && Z(pos)>=-tr && rho>=(a_n[n]-p_n[n]/8.0) && rho<=a_n[n]) );

}

// return 1 if we crossed the two or bottom sheet when moving between pos1, pos2 (NOT USED RIGHT NOW)
int crossed_tb_sheet(mwvector pos1, mwvector pos2, int n) {
  real rho = hypot(X(pos1),Y(pos1));
  real R=(3.0-mw_sqrt(1.0+(8.0/p_n[n])*(rho-a_n[n])))/4.0;
  real tr=2.0*p_n[n]*mw_sqrt(cube(R)*(1.0-R));
  int in1 = (fabs(Z(pos1))<=tr && rho>=a_n[n] && rho<=a_n[n]+p_n[n]);

  rho = hypot(X(pos2),Y(pos2));
  R=(3.0-mw_sqrt(1.0+(8.0/p_n[n])*(rho-a_n[n])))/4.0;
  tr=2.0*p_n[n]*mw_sqrt(cube(R)*(1.0-R));
  int in2 = (fabs(Z(pos2))<=tr && rho>=a_n[n] && rho<=a_n[n]+p_n[n]);

  return (in1 != in2);

}

// returns number of collisions
// collision_info holds a max of 4 collisions, and each collision holds rho_s, z_s, and T_s for hitting the surface of the caustic ring
int try_caustic_collision(real rho1, real rho2, real z1, real z2, int n, real collision_info[4][3]) {
  real X1 = (rho1-a_n[n])/p_n[n];
  real X2 = (rho2-a_n[n])/p_n[n];
  real Z1 = z1/p_n[n];
  real Z2 = z2/p_n[n];
  real m = (Z2-Z1)/(X2-X1);
  real b = Z2-m*X2;

  real mpow4 = m*m*m*m;
  real mpow5 = mpow4*m;
  real mpow6 = mpow5*m;
  real mpow7 = mpow6*m;

  real bpow4 = b*b*b*b;
  real bpow5 = bpow4*b;
  real bpow6 = bpow5*b;


  real A = -(3*m*m+1)/(2*m*m+2);
  real B = -(13*m*m+4*b*m)/(3*m*m+3);
  real C = pow(2*mpow6-48*b*mpow5+384*b*b*mpow4-72*mpow4-1024*b*b*b*m*m*m+648*b*m*m*m-432*b*b*m*m+432*m*m-1152*b*b*b*m+864*b*m+432*b*b \
                 +sqrt(-6912*b*mpow7+89856*b*b*mpow6-6912*mpow6-131328*b*b*b*mpow5+269568*b*mpow5-1002240*bpow4*mpow4+518400*b*b*mpow4 \
                      +186624*mpow4-1216512*bpow5*m*m*m-573696*b*b*b*m*m*m+746496*b*m*m*m-442368*bpow6*m*m-2032128*bpow4*m*m \
                      +1119744*b*b*m*m-1658880*bpow5*m+746496*b*b*b*m-442368*bpow6+186624*bpow4),1.0/3.0);
  real D = 1.0/(12*pow(2,1.0/3.0)*(m*m+1));
  real E = (mpow4-16*b*m*m*m+64*b*b*m*m-24*m*m+24*b*m+48*b*b)/(6*pow(4,1.0/3.0)*(m*m+1));
  real F = -12*(m*m+b*m)/(m*m+1);

  real sqrt1 = sqrt(A*A+B/2+C*D+E/C);
  real T[6]; // if all T[0] through T[3] nan, we are on the line z=0 and T = 0 or 1, so we add two spots with 0 and 1 accessed only during such nans
  T[0] = -0.5*( A + sqrt1 + sqrt(2*A*A+B-C*D-E/C+(8*A*A*A+6*A*B+F)/(4*sqrt1)) );
  T[1] = -0.5*( A + sqrt1 - sqrt(2*A*A+B-C*D-E/C+(8*A*A*A+6*A*B+F)/(4*sqrt1)) );
  T[2] = -0.5*( A - sqrt1 + sqrt(2*A*A+B-C*D-E/C-(8*A*A*A+6*A*B+F)/(4*sqrt1)) );
  T[3] = -0.5*( A - sqrt1 - sqrt(2*A*A+B-C*D-E/C-(8*A*A*A+6*A*B+F)/(4*sqrt1)) );
  T[4] = 0;
  T[5] = 1;
  int i;
  real X_s, Z_s;
  int collnum = 0;
  for (i = 0; i < 6; ++i) {
    X_s = (T[i]-1)*(2*T[i]-1);
    Z_s = 2*sqrt(T[i]*T[i]*T[i]*(1-T[i]));
    if (( (X1 <= X_s)&&(X_s <= X2) )||( (X1 >= X_s)&&(X_s >= X2) )) {
      if (!(( (Z1 < Z_s)&&(Z_s < Z2) )||( (Z1 > Z_s)&&(Z_s > Z2) ))) Z_s = -Z_s;
      if (!(( (Z1 < Z_s)&&(Z_s < Z2) )||( (Z1 > Z_s)&&(Z_s > Z2) ))) continue;
      collision_info[collnum][0] = X_s*p_n[n]+a_n[n];
      collision_info[collnum][1] = Z_s*p_n[n];
      collision_info[collnum][2] = T[i];
      ++collnum;
    }
  }
  return collnum;
}

void apply_dynamical_friction_voluminous(NBodyState* st) {
  const int nbody = st->nbody;

  mwvector* accs = mw_assume_aligned(st->acctab, 16);
  mwvector dwarf_vel; X(dwarf_vel) = Y(dwarf_vel) = Z(dwarf_vel) = 0;
  int i;

  // get velocity of dwarf's center of mass
  for (i = 0; i < nbody; ++i) {
    Body* b = &st->bodytab[i];
    X(dwarf_vel) += X(Vel(b))/nbody;
    Y(dwarf_vel) += Y(Vel(b))/nbody;
    Z(dwarf_vel) += Z(Vel(b))/nbody;
  }
  real dwarf_mass = 10.0;
  real b_max = 1;
  int n;
  for (i = 0; i < nbody; ++i) {
    Body* b = &st->bodytab[i];
    real rho = hypot(X(Pos(b)),Y(Pos(b)));
    real z = Z(Pos(b));
    for (n = 1; n <= 20; ++n) {
      if (in_caustic_envelope(Pos(b),n) == 1) {
        mwvector V_flow;
        X(V_flow) = -(V_n[n]*1.0226831)*(Y(Pos(b))/rho);
        Y(V_flow) =  (V_n[n]*1.0226831)*(X(Pos(b))/rho);
        Z(V_flow) = 0;

        mwvector V0;
        X(V0) = X(V_flow) - X(dwarf_vel);
        Y(V0) = Y(V_flow) - Y(dwarf_vel);
        Z(V0) = Z(V_flow) - Z(dwarf_vel); // Z(V_flow) technically unneccessary
        real V0_mag = sqrt(X(V0)*X(V0)+Y(V0)*Y(V0)+Z(V0)*Z(V0));
        real lambda = b_max*V0_mag*V0_mag/(G_caustics*dwarf_mass);
        real factor = 4*M_PI*G_caustics*G_caustics;
        real acc_scalar = factor*log(lambda)*dwarf_mass*get_density_close(rho,z,n) / (V0_mag*V0_mag*V0_mag); // THIS IS NOT THE MAGNITUDE OF ACC
        printf("dwarf velocity: |(%11f,%11f,%9f)|=%11f --- flow velocity: |(%11f,%11f,%8f)|=%11f --- acc: |(%10f,%10f,%10f)|=%10f \n",
                                                       X(dwarf_vel),     Y(dwarf_vel),     Z(dwarf_vel),
                                                       sqrt(X(dwarf_vel)*X(dwarf_vel)+Y(dwarf_vel)*Y(dwarf_vel)+Z(dwarf_vel)*Z(dwarf_vel)),
                                                       X(V_flow),        Y(V_flow),        Z(V_flow),        V_n[n]*1.0226831,
                                                       acc_scalar*X(V0), acc_scalar*Y(V0), acc_scalar*Z(V0), acc_scalar*V0_mag);
        X(accs[i]) += acc_scalar*X(V0);
        Y(accs[i]) += acc_scalar*Y(V0);
        Z(accs[i]) += acc_scalar*Z(V0);
      }
    }
  }
}

void apply_dynamical_friction(NBodyState* st, const NBodyCtx* ctx) {
  real dwarf_mass = 10.0;
  real b_max = 1; // doesn't matter much, probably
  real factor = 2*M_PI*G_caustics*G_caustics*dwarf_mass;
  real envelope_density[21];
  real avg_delta_V_sum1 = 0;

  // area density of caustic envelope assuming ALL mass in the caustic is evenly distributed on the envelope
  int n;
  for (n = 1; n <= 20; ++n) {
    envelope_density[n] = 0.45*sqrt(3)*(rate_n[n]*4498.6589)/(V_n[n]*1.0226831)/(a_n[n]+p_n[n]/4);
  }

  // delta vel calculations
  const int nbody = st->nbody;
  int i;
  for (i = 0; i < nbody; ++i) {
    Body* b = &st->bodytab[i];
    real rho = hypot(X(Pos(b)),Y(Pos(b)));
    mwvector old_pos = mw_subv( Pos(b), mw_mulvs(Vel(b), (ctx->timestep)) );
    real collision_info[4][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}}; // coordinates of up to 4 collisions - rho, Z, and T in that order (T just for which sheet)
    int collnum; // number of collisions with the caustic envelope
    mwvector delta_vel = mw_vec(0,0,0);
    for (n = 1; n <= 20; ++n) {
      collnum = try_caustic_collision(rho,hypot(X(old_pos),Y(old_pos)),Z(Pos(b)),Z(old_pos),n,collision_info);
      int coll;
      for (coll = 0; coll < collnum; ++coll) {

        mwvector V_flow; // velocity vector of axion flow
        X(V_flow) =  (V_n[n]*1.0226831)*(Y(Pos(b))/rho);
        Y(V_flow) = -(V_n[n]*1.0226831)*(X(Pos(b))/rho);
        Z(V_flow) = 0;

        mwvector V0 = mw_subv(V_flow,Vel(b)); // relative velocity vector b/w axions and dwarf
        real V0_mag = mw_absv(V0);

        real lambda = b_max*V0_mag*V0_mag/(G_caustics*dwarf_mass);

        real rho_s = collision_info[coll][0];
        real z_s   = collision_info[coll][1];
        real T_s   = collision_info[coll][2];
        real phi   = (atan(Y(Pos(b))/X(Pos(b)))+atan(Y(old_pos)/X(old_pos)))/2;
        real theta;
        if (T_s < 0.75) {
          theta = (1-2*(z_s < 0))*atan(sqrt((1+sqrt(8*(rho_s-a_n[n])/p_n[n]+1))
                                           /(3-sqrt(8*(rho_s-a_n[n])/p_n[n]+1))));
          if (theta != theta) theta = (1-2*(z_s < 0))*M_PI/2; // theta nans when rho=a_n[n]+p_n[n], when this theta should be used.
        }
        else {
          theta = (1-2*(z_s < 0))*atan(sqrt((1-sqrt(8*(rho_s-a_n[n])/p_n[n]+1))
                                           /(3+sqrt(8*(rho_s-a_n[n])/p_n[n]+1))));
        }
        mwvector proj_vec = mw_vec(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
        real V_proj = mw_dotv(Vel(b),proj_vec);

        real dv_scalar = factor*log(1+lambda*lambda)*envelope_density[n] / (V0_mag*V0_mag*V0_mag) / V_proj; // THIS IS NOT THE MAGNITUDE OF DELTA V.
        delta_vel = mw_addv(delta_vel, mw_mulvs(V0,dv_scalar));
      }
      if (collnum != 0) {
        if (delta_V_max < mw_absv(delta_vel)) delta_V_max = mw_absv(delta_vel);
        avg_delta_V_sum1 += mw_absv(delta_vel);
        printf("colls: %d -- vel: (x=%7.2f,y=%7.2f,z=%7.2f)=%6f -- pos: (rho=%7.2f,theta=%7.2f,z=%7.2f) -- delta vel: %3.2f -- max: %5.2f -- avgsum: %7f\n",
                                                       collnum,
                                                       X(Vel(b)),        Y(Vel(b)),        Z(Vel(b)),
                                                       sqrt(X(Vel(b))*X(Vel(b))+Y(Vel(b))*Y(Vel(b))+Z(Vel(b))*Z(Vel(b))),
                                                       X(Pos(b)),     Y(Pos(b)),     Z(Pos(b)),
                                                       mw_absv(delta_vel), delta_V_max, avg_delta_V_sum2);
      }
      Vel(b) = mw_addv(Vel(b), delta_vel);
    }
  }
  avg_delta_V_sum2 += avg_delta_V_sum1/nbody;
}

 
mwvector causticHaloAccel(const Halo* h, mwvector pos, real r)
{
    mwvector accel;

    real rho=0.0, rfield, zfield;

    rfield = 0.0;
    zfield = 0.0;

    rho = mw_sqrt(sqr(X(pos))+sqr(Y(pos)));

    // rho cannot be zero (causes nan in near field and ax and ay at origin)
    if (rho < 0.000001) {
      rho = 0.000001;
    }

    int n;
    for (n = 1; n <= 20; n++)
    {

        if( in_caustic_envelope(pos,n) )
        {
            gfield_close(rho,Z(pos),n,&rfield,&zfield);
        }

        else
        {
            gfield_far(rho,Z(pos),n,&rfield,&zfield);        
        }
    }

    X(accel) = (rfield*X(pos))/rho;
    Y(accel) = (rfield*Y(pos))/rho;
    Z(accel) = zfield;
    
    return accel;
}
