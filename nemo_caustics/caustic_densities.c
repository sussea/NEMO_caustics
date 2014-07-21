#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

//local double G = 1.0;
const double G_caustics = 1.0;

// properties of n=1-20 caustic ring flows from tables in Duffy & Sikivie (2008)
// caustic flow number     1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15,   16,   17,   18,   19,   20    
double a_n[]    = {1.0, 40.1, 20.1, 13.6, 10.4,  8.4,  7.0,  6.1,  5.3,  4.8,  4.3,  4.0,  3.7,  3.4,  3.2,  3.0,  2.8,  2.7,  2.5,  2.4,  2.3}; // caustic ring radius (kpc)
double V_n[]    = {1.0,  517,  523,  523,  523,  522,  521,  521,  520,  517,  515,  512,  510,  507,  505,  503,  501,  499,  497,  496,  494}; // particle speed in caustic (km/s)
double p_n[]    = {1.0,  0.3,  0.3,  1.0,  0.3, 0.15, 0.12,  0.6, 0.23, 0.41, 0.25, 0.19, 0.17, 0.11, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09}; // rho width of caustic (kpc)
double rate_n[] = {1.0,   53,   23,   14,   10,  7.8,  6.3,  5.3,  4.5,  3.9,  3.4,  3.1,  2.8,  2.5,  2.3,  2.1,  2.0,  1.8,  1.7,  1.6,  1.5}; // dark matter mass infall per unit solid angle per time (M_sun/sterad*yr)

const double ONE_THIRD=1.0/3.0;
#define c1 (double complex) 1.0
#define c2 (double complex) 2.0
#define c3 (double complex) 3.0
#define c4 (double complex) 4.0
#define c12 (double complex) 12.0
#define c05 (double complex) 0.5
#define II (double complex) 1.0*I
const double TOLERANCE_caustics=0.00001;
const double TOLERANCE2_caustics=0.01;
double printed_density_cap = 1000;

#define X_caustics 0
#define Y_caustics 1
#define Z_caustics 2

// let's compute the values of T at which the various poles cross the real axis
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
    return c05 * (c1 - csqrt(f1(x) / 2.0 + temp) + csqrt( f1(x) - temp - c2 * x / csqrt(f1(x) / 2.0 + temp) ));
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

// returns number of useful values in T (number of flows)
int getTVals(double rho, double z, int n, double* T){

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

// If T[flow] == 0 then the z equation is degenerate, so get alpha_factor from the rho equation. Otherwise, use the z equation.
// (it seems you can entirely avoid the z equation and still get the right answer?)
double get_density_close(double rho, double z, int n) {
  double density = 0;
  double T[4];
  int T_count = getTVals(rho,z,n,T);
  int flow;
  for (flow = 0; flow < T_count; ++flow) {
    double alpha_factor; // alpha_factor is (s*alpha^2)/2
    if (fabs(T[flow]) > 0.0000001) {
      alpha_factor = pow(z/T[flow],2)/(4*p_n[n]);
    } else {
      alpha_factor = a_n[n] - rho + p_n[n]*(T[flow]-1)*(T[flow]-1); // if negative, alpha*tao_0 is imaginary so we overcounted the flows/poles, ignore this one
      if ((fabs(alpha_factor) > 0.00001) && (alpha_factor < 0)) {
        continue;
      }
    }
    double D = 2*(V_n[n]*1.0226831)*(p_n[n]*T[flow]*(T[flow]-1) + alpha_factor);//convert from km/s to kpc/gyr with 1.0226
    density += fabs(4498.6589*rate_n[n]/(a_n[n]*D)); //convert from M_sun/yr*sr to m_sim/gyr*sr with 4498
    // ^ more accurate to replace a_n[n] with rho. a_n[n] is used only to match estimate on page 5 of A. Natarajan, P. Sikivie, Phys. Rev.D76, 023505 (2007)
  } // in fact, final version should use rho, and you should not see this message.
  if (density > printed_density_cap) density = printed_density_cap;
  return density;
}

double get_density_far(double rho, double z, int n) {
  double r_squared = rho*rho + z*z;
  double shift = a_n[n]+p_n[n]/4;
  double s = hypot(r_squared - shift*shift, 2.0 * shift * z);
  double factor = (8.0 * rate_n[n] * 4498.6589) / (V_n[n] * 1.0226831);  //convert from M_sun/yr*sr to m_sim/gyr*sr with 4498 and km/s to kpc/gyr with 1.0226
  double density = factor*shift*shift*r_squared/(s*pow(2*shift*shift+s,2));
  if (density > printed_density_cap) density = printed_density_cap;
  return density;
}

// this is stupid - ignore it
double get_density_far2(double rho, double z, int n) {
  double density = 0;
  double step_velocity = p_n[n]/3000;
  double T;
  int steps = 0;
  for (T = 0; T < 1; T += step_velocity * (T*T*T)*(T-1) / ( (4*T-3)*sqrt(pow(T,6)*(T-1)*(T-1) - 1) )) {
    double shift = 2*a_n[n] + p_n[n]*(T-1)*(2*T-1);
    double r_squared = rho*rho + pow(z-2*p_n[n]*T*sqrt(T*(1-T)),2);
    double s = hypot(r_squared - shift*shift, 2.0 * shift * z);
  double factor = (8.0 * rate_n[n] * 4498.6589) / (V_n[n] * 1.0226831);  //convert from M_sun/yr*sr to m_sim/gyr*sr with 4498 and km/s to kpc/gyr with 1.0226
    density += factor*shift*shift*r_squared/(s*pow(2*shift*shift+s,2));
    ++steps;
  }
  density = density / steps;
  return density;
}

int in_caustic_envelope(double rho, double z, int n) {

  double r, l, tr, tl; // r (right), l (left), tr (top right), tl (top left)
  //the caustic ring has a tricusp boundary (see Tam 2012)
  r  = (3.0 - sqrt( 1.0 + (8.0 / p_n[n]) * (rho - a_n[n]) )) / 4.0;
  l  = (3.0 + sqrt( 1.0 + (8.0 / p_n[n]) * (rho - a_n[n]) )) / 4.0;
  tr = 2.0 * p_n[n] * sqrt(pow(r, 3.0) * (1.0 - r));
  tl = 2.0 * p_n[n] * sqrt(pow(l, 3.0) * (1.0 - l));

  //if position (rho,z) is inside ring use gfield_close, else use gfield_far
  return ( (z <= tr && z >= 0.0 && rho >= a_n[n] && rho <= a_n[n] + p_n[n])
   || (z >= tl  && z <= tr  && rho >= (a_n[n] - p_n[n] / 8.0) && rho <= a_n[n])
   || (z >= -tr && z <= 0.0 && rho >= a_n[n] && rho <= a_n[n] + p_n[n])
   || (z <= -tl && z >= -tr && rho >= (a_n[n] - p_n[n] / 8.0) && rho <= a_n[n]) );

}

double get_density(double rho, double z) {
  double density = 0;
  int n;

  // rho cannot be zero (causes nan in near field and ax and ay at origin)
  if (rho < 0.000001) {
    rho = 0.000001;
  }
  // calculate gfield at a position (x,y,z) by adding the contributions from all n caustic ring flows
  for (n = 1; n <= 20; ++n) {

    if (in_caustic_envelope(rho, z, n)) {
      density += get_density_close(rho, z, n);
    } else {
      density += get_density_far(rho, z, n);
    }
  }
  return density;
}

// incredibly rough
void apply_dynamical_friction(double* pos, double* vel, double* acc, double M, double density, int n) {
  double dwarf_mass = 10.0 / 100000; //m_sim
//  double line_density = 0.6*p_n[n]*p_n[n]*get_density_close(a_n[n]+p_n[n]/4,0,n);
//  double density = line_density / (p_n[n]*p_n[n]);
  density = 500; //m_sim/kpc^3 m_solar (Julie has 10^7 M_sun/kpc^3 = 45 m_sim/kpc^3)
  double relvel = 1.0226831*(V_n[n] - 400); //kpc/gyr
  double b_max; //kpc
  for (b_max = 0; b_max < 10; b_max += 0.1) { 
    double lambda = b_max * relvel * relvel / (G_caustics * dwarf_mass);
    double accel = -4*M_PI*log(lambda)*G_caustics*G_caustics*density*dwarf_mass/relvel/relvel;
    double velchange = accel/relvel;
    printf("b_max: %f | vel: %e kpc/gyr | acc: %e kpc/gyr^2 | velchange: %e (kpc/gyr)/kpc\n", b_max, relvel, accel, velchange);
  }
}

// you can't numerically integrate an infinite spike you dolt! (this function tries to do that)
double get_caustic_mass(double rho, double z, int n, double integration_length) {
    if (in_caustic_envelope(rho,z,n)) {
      double temp = printed_density_cap;
      printed_density_cap = 9999999999999;
      double integration_steps = 100;
      double partition_length = integration_length/integration_steps;
      double partition_volume = partition_length*partition_length*integration_length;

      double slice_mass = 0.6*p_n[n]*p_n[n]*get_density_close(a_n[n]+p_n[n]/4,0,n)*integration_length; // total mass of caustic slice, derived from eq 6 on page 5 of arXiv:0705.0001v1 [astro-ph]
      double density_cap = slice_mass/partition_volume; // if hit, cap should still overestimate avg density of the partition

      double densities_sum_rho_z = 0;

      double rho_i, z_i;

      // integrate through rho
      for (rho_i = rho - integration_length/2; rho_i < rho + integration_length/2; rho_i += partition_length) {

        // subtract half the top and bottom (endpoint) partitions for trapezoid rule (z integral)
        double partition_density_top = get_density_close(rho_i, z + integration_length/2, n);
        double partition_density_bot = get_density_close(rho_i, z - integration_length/2, n);
        double densities_sum_z = (partition_density_top < density_cap ? partition_density_top : density_cap);
        densities_sum_z       += (partition_density_bot < density_cap ? partition_density_bot : density_cap);
        densities_sum_z = -0.5*densities_sum_z;

        // integrate through z
        for (z_i = z - integration_length/2; z_i < z + integration_length/2; z_i += partition_length) {
          double partition_density = get_density_close(rho_i, z_i, n);
          if (partition_density > density_cap) partition_density = density_cap;
          densities_sum_z += partition_density;
        }

        // this should never happen, but just in case...
        if (densities_sum_z < 0) {
          printf("SOMETHING WENT WRONG: Function 'get_caustic_mass' in 'mpc.c' assigned negative density %f to a partition at rho=%f z=%f.\nTerminating program.\n",densities_sum_z,rho,z);
          exit(1);
        }
        double trapezoid_factor = (fabs(rho_i - integration_length/2) < partition_length ? 0.5 : 1.0); // divide end points by 2 (rho integral)
        densities_sum_rho_z += densities_sum_z*trapezoid_factor; // separating sums prevents small numbers from being lost to numerical errors

      }
      printed_density_cap = temp;
      return (densities_sum_rho_z * partition_volume);
    } else {
      return (get_density_far(rho, z, n)*pow(integration_length,3)); // no density jumps, so the density is approx the same throughout the volume
    }
}

void main() {
//  apply_dynamical_friction();
//  exit(0);
  int important_n = 3;
  double rho_min = a_n[important_n] - p_n[important_n]/10;//0;
  double rho_max = a_n[important_n] + 4*p_n[important_n] + p_n[important_n]/10;//a_n[1]+3*p_n[1];
  int steps = 5000;
  double integration_length = (rho_max-rho_min)/steps;
  double rho;
  double avgdensity = 0;
  double mass = 0;
  printf("rho = %f to %f   |   integration_length = %f\n",rho_min,rho_max,integration_length);
  for (rho = rho_min; rho < rho_max; rho += integration_length) {
    double z = 0.001;
    double density_full = get_density(rho,z); // uses close and far density equations depending on proximity to each caustic
    double density_close = 0; // only uses close equation, even where it is inaccurate
    double density_far = 0;   // only uses far equation, even where it is inaccurate
    double density_far2 = 0;  // only uses far2 equation, even where it is inaccurate
    int n;
    for (n = 1; n <= 20; ++n) {
//      density_close += get_density_close(rho,z,n);
//      density_far += get_density_far(rho,z,n);
//      density_far2 += get_density_far2(rho,z,n);
      mass += get_caustic_mass(rho,z,n,integration_length);
    }

    avgdensity = mass/pow(integration_length,3);
    printf("%f %f %f %f %f %f %f %f\n", rho, z, density_full* 1e9/(4498.6589), density_close, density_far, density_far2, mass, density_full - avgdensity);
  }
}


