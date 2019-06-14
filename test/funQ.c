#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

//#define MODEL 1

static double FUN_Q (double Q, void *ptr);
static double D_FUN_Q (double Q, void *ptr);

struct FUN_Q_params
{
  double d;
  double M1;
  double M2;
  double M3;
  double E;
};

double Gamma = 1.66666;
double Gamma_1 = 0.66666;

int
main ()
{
   FILE *fptr1, *fptr;
#if (MODEL==0)
   fptr = fopen("type0.dat","w");
#else
   fptr1 = fopen("type1.dat","w");
#endif

 double d = 3.093660e+01;
 double M1 = 5.351120e+03;
 double M2 = -1.357348e-01;
 double M3 = 2.377042e-01;
 double E = 5.349582e+03;

 double Msqr=SQR(M1)+SQR(M2)+SQR(M3);
// double d =sqrt(Msqr)/Gamma_1;

/*
 if (sqrt(Msqr)>E){
   printf("sqrt(Msqr)>E\n");
   abort();
  }
*/
  struct FUN_Q_params params = { d, M1, M2, M3, E };

  double f;

 double Ql = -10000;
  double Qr = 10000;
  int N = 10000;
  double dQ = (Qr - Ql)/N;
  double Q = Ql;

  for (int i = 0; i <= N-1; i++)
    {
      Q = dQ + Q;

      f = FUN_Q (Q, &params);

#if (MODEL == 0)
      fprintf (fptr, "%e %e %e %e\n", Q, f, D_FUN_Q (Q, &params), f-Q/Gamma + E);
#else
      fprintf (fptr1, "%e %e %e %e\n", Q, f, D_FUN_Q (Q, &params), f-Q/Gamma + E);
#endif
    }

  double Q0 = 1.0e+32;

//  printf( "Q0=%e, D_FUN_Q=%e\n" ,Q0, D_FUN_Q (Q0, &params));

  return 0;
}


#if (MODEL==0)
static double
FUN_Q (double Q, void *ptr)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  double d = (params->d);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
  double E = (params->E);

  double U1 = M1 / Q;
  double U2 = M2 / Q;
  double U3 = M3 / Q;

  double rho = d / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3));

  double pres = (Gamma_1 / Gamma) * (Q / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - rho);

  double f = Q * sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - pres - E;

  return f;
}

static double
D_FUN_Q (double Q, void *ptr)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  double d = (params->d);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);

  double U1 = M1 / Q;
  double U2 = M2 / Q;
  double U3 = M3 / Q;

  double dU1 = -M1 / (Q * Q);
  double dU2 = -M2 / (Q * Q);
  double dU3 = -M3 / (Q * Q);

  double dd = -d * pow (1 + SQR (U1) + SQR (U2) + SQR(U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3);

  double dp = (Gamma_1 / Gamma) * (1 / sqrt (1 + SQR (U1) + SQR (U2) + SQR(U3))
		- Q * pow (1 + SQR (U1) + SQR (U2) + SQR(U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3) - dd);

  double df = sqrt (1 + SQR (U1) + SQR (U2) + SQR(U3)) + Q * (U1 * dU1 + U2 * dU2 + U3 * dU3) 
            / sqrt (1 + SQR (U1) + SQR (U2) + SQR(U3)) - dp;

  return df;
}


#else // MODEL == 1
static double
FUN_Q (double Q, void *ptr)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  double d = (params->d);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
  double E = (params->E);

  double Msqr = SQR(M1)+SQR(M2)+SQR(M3);
  double Qsqr = SQR(Q);

  //double f = (1.0/(Gamma*sqrt(Msqr+Qsqr)))*(Gamma*SIGN(Q)*Msqr+SIGN(Q)*Qsqr+d*fabs(Q)*Gamma_1)-E;
  double f = (1.0/(Gamma*sqrt(Msqr+Qsqr)))*(Gamma*Msqr+Qsqr+d*Q*Gamma_1)-E;

  return f;
}

static double
D_FUN_Q (double Q, void *ptr)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  double d = (params->d);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
 
  double Msqr = SQR(M1)+SQR(M2)+SQR(M3);
  double Qsqr = SQR(Q);  

  double a = (2-Gamma)*Msqr*Q;
  double b = d*(Gamma_1)*Msqr;

  double df = (Qsqr*Q+a+b)/(Gamma*pow(Msqr+Qsqr, 1.50));

  return df;
}
#endif
