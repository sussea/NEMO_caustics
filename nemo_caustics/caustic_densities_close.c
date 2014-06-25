#include <math.h>
#include <complex.h>
#include <stdio.h>

//local double G = 1.0;
const double G_caustics = 1.0;

// properties of n=1-20 caustic ring flows from tables in Duffy & Sikivie (2008)
// caustic flow number     1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15,   16,   17,   18,   19,   20    
double a_n[]    = {1.0, 40.1, 20.1, 13.6, 10.4,  8.4,  7.0,  6.1,  5.3,  4.8,  4.3,  4.0,  3.7,  3.4,  3.2,  3.0,  2.8,  2.7,  2.5,  2.4,  2.3}; // caustic ring radius
double V_n[]    = {1.0,  517,  523,  523,  523,  522,  521,  521,  520,  517,  515,  512,  510,  507,  505,  503,  501,  499,  497,  496,  494}; // particle speed in caustic
double p_n[]    = {1.0,  0.3,  0.3,  1.0,  0.3, 0.15, 0.12,  0.6, 0.23, 0.41, 0.25, 0.19, 0.17, 0.11, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09}; // longitudinal length of caustic (along rho)
double rate_n[] = {1.0,   53,   23,   14,   10,  7.8,  6.3,  5.3,  4.5,  3.9,  3.4,  3.1,  2.8,  2.5,  2.3,  2.1,  2.0,  1.8,  1.7,  1.6,  1.5}; // dark matter mass infall per unit solid angle per time

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
      if (T[count] < 0.0000001) T[count] = 0;
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

double get_density(double rho, double z, int n) {
  double density = 0;
  double T[4];
  double T_count = getTVals(rho,z,n,T);

  double factor = -rate_n[n]/(a_n[n]*V_n[n]);
  int flow;
  for (flow = 0; flow < T_count; ++flow) {
    double s  = 2*(a_n[n] - rho + p_n[n]*(T[flow]-1)*(T[flow]-1)) / (a_n[n]*a_n[n]);
    double u  = V_n[n]*V_n[n]/s;
    double t0_squared = 2*p_n[n]/u;
    double angle;
    if (T[flow] != 0)
      angle = z / (T[flow]*sqrt(t0_squared)*V_n[n]);
    else
      angle = sqrt((2*(a_n[n] - rho) + u*t0_squared*(T[n]-1)*(T[n]-1))/s);
    density += fabs(factor*cos(0)/(2*p_n[n]*T[flow]*(T[flow]-1) + s*angle*angle));
//    printf("s: %5g  |  u: %5g  |  t0_squared: %5g  |  angle %5g  |  angle_factor: sqrt(%5g)  |  T: %5g\n",s, u, t0_squared, angle, angle*angle*s/p_n[n], T[flow]);
  }
  return density;
}

void main() {
  int n = 1;
  double rho = a_n[n] + p_n[n]/4;
  double z   = 0;
  double density = get_density(rho, z, n);
  double correct_density = 4*rate_n[n]/(a_n[n]*V_n[n]*p_n[n]);
  printf("Generated: %.16f\nCorrect:   %.16f\n", density, correct_density);
}

