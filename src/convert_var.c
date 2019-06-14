#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "defs.h"
#include "struct.h"
#include "globals.h"
#include "prototypes.h"

struct FUN_Q_params
{
  double d;
  double M1;
  double M2;
  double M3;
  double E;
};

static double FUN_Q (double, void *);
static double D_FUN_Q (double, void *);
static void FDF_FUN_Q (double, void *, double *, double *);

Prim1DS Cons1D_to_Prim1D (const Cons1DS *);
Cons1DS Prim1D_to_Cons1D (const Prim1DS *);

PrimS
Cons_to_Prim (const ConsS * pCons)
{
  Cons1DS U;
  Prim1DS W;
  PrimS Prim;

  U.d = pCons->d;
  U.Mx = pCons->M1;
  U.My = pCons->M2;
  U.Mz = pCons->M3;
  U.E = pCons->E;

  W = Cons1D_to_Prim1D (&U);

  Prim.d = W.d;
  Prim.U1 = W.Ux;
  Prim.U2 = W.Uy;
  Prim.U3 = W.Uz;
  Prim.P = W.P;

  return Prim;
}

ConsS
Prim_to_Cons (const PrimS * pW)
{
  Cons1DS U;
  Prim1DS W;
  ConsS Cons;

  W.d = pW->d;
  W.Ux = pW->U1;
  W.Uy = pW->U2;
  W.Uz = pW->U3;
  W.P = pW->P;

  U = Prim1D_to_Cons1D (&W);

  Cons.d = U.d;
  Cons.M1 = U.Mx;
  Cons.M2 = U.My;
  Cons.M3 = U.Mz;
  Cons.E = U.E;

  return Cons;
}

Prim1DS
Cons1D_to_Prim1D (const Cons1DS * cons)
{
#ifdef SR_DEBUG
  if ((cons->d < 0) || (cons->E < 0))
    {
      printf ("\n\nerror: D < 0 or E < 0!\n");
      printf ("file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
      printf ("line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, cons->d, cons->Mx, cons->My, cons->Mz, cons->E);
      abort ();
    }
#endif

  Prim1DS prim;

  double Msqr = SQR (cons->Mx) + SQR (cons->My) + SQR (cons->Mz);
  double M    = sqrt(Msqr);

#ifdef SR_DEBUG
  if (SQR(cons->E) <= Msqr + SQR(cons->d))
    {
      printf ("\n\nerror: E^2 <= |M|^2 + D^2!\n");
      printf ("file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
      printf ("line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, cons->d, cons->Mx, cons->My, cons->Mz, cons->E);
      printf ("E^2-|M|^2-D^2=%e\n\n", SQR( cons->E) - Msqr - SQR(cons->d));
      abort ();
    }
#endif

  if (fabs (M) > TINY_NUMBER)
    {
      int status;

      int iter = 0;
      int max_iter = 100;

      const gsl_root_fdfsolver_type *T;

      gsl_root_fdfsolver *s;

      double Q0;
      double Q;

/* initial guess Q  */
      if (cons->d > M / Gamma_1)
	{
	  Q = M * (cons->E - M) / ((1 - 1 / Gamma) * cons->d);
	}
      else
	{
	  Q = cons->E * Gamma;
	}

      gsl_function_fdf F;
      struct FUN_Q_params params = { cons->d, cons->Mx, cons->My, cons->Mz, cons->E };

      F.f = &FUN_Q;
      F.df = &D_FUN_Q;
      F.fdf = &FDF_FUN_Q;
      F.params = &params;

      T = gsl_root_fdfsolver_newton;

      s = gsl_root_fdfsolver_alloc (T);

      gsl_root_fdfsolver_set (s, &F, Q);

      //printf ("status = %s\n", gsl_strerror (status)); 
      do
	{
	  iter++;
	  status = gsl_root_fdfsolver_iterate (s);
	  //printf ("status = %s\n", gsl_strerror (status));
	  Q0 = Q;
	  //printf("Q=%20.16e\n",Q);
	  Q = gsl_root_fdfsolver_root (s);
	  status = gsl_root_test_delta (Q, Q0, 0, 1e-16);
	  //printf ("status = %s\n", gsl_strerror (status));
	}
      while (status == GSL_CONTINUE && iter < max_iter);
      //printf ("status = %s\n", gsl_strerror (status));
      prim.Ux = cons->Mx / Q;
      prim.Uy = cons->My / Q;
      prim.Uz = cons->Mz / Q;

      double Factor = sqrt (1 + SQR (prim.Ux) + SQR (prim.Uy) + SQR (prim.Uz));

      prim.d = cons->d / Factor;
      prim.P = (Gamma_1 / Gamma) * (Q / Factor - prim.d);

      gsl_root_fdfsolver_free (s);
#ifdef SR_DEBUG
      if (prim.P < 0)
	{
	  printf ("\n\nerror: P < 0!\n");
	  printf ("file: %s\nfunction: %s\nline:%d \n", __FILE__, __FUNCTION__, __LINE__);
	  printf ("Q = %e\n", Q);
	  printf ("D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", cons->d, cons->Mx, cons->My, cons->Mz, cons->E);
	  printf ("|M| > TINY_NUMBER!!\n");
	  abort ();
	}
#endif
    }
  else if (cons->E >= cons->d)
    {
      prim.Ux = 0.0;
      prim.Uy = 0.0;
      prim.Uz = 0.0;
      prim.d = cons->d;
      prim.P = Gamma_1 * (cons->E - cons->d);
    }
  else
    {
      printf ("\nToo critical to solve!\n");
      printf ("D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", cons->d, cons->Mx, cons->My, cons->Mz, cons->E);
      abort();
    }
  return prim;
}

Cons1DS
Prim1D_to_Cons1D (const Prim1DS * prim)
{
  double Ux = prim->Ux;
  double Uy = prim->Uy;
  double Uz = prim->Uz;
  double rho = prim->d;
  double P = prim->P;

  Cons1DS cons;

  double Usqr = 1 + SQR (Ux) + SQR (Uy) + SQR (Uz);
  double Usqrt = sqrt (Usqr);
  double a = (rho + (Gamma / Gamma_1) * P);
  double b = a * Usqrt;

  cons.d = rho * Usqrt;
  cons.Mx = b * Ux;
  cons.My = b * Uy;
  cons.Mz = b * Uz;
  cons.E = a * Usqr - P;

  return cons;
}

/*------ The following function have high stability in ultra-relativistic regime---------*/
/*
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
*/
/*------ The following function have high stability in ultra-relativistic regime---------*/
/*
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
*/
/*------ The following function have high stability in ultra-relativistic regime---------*/
/*
static void
FDF_FUN_Q (double Q, void *ptr, double *f, double *df)
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

  double rho = d / sqrt (1 + SQR (U1) + SQR (U2) + SQR(U3));
  double pres = (Gamma_1 / Gamma) * (Q / sqrt (1 + SQR (U1) + SQR (U2) + SQR(U3)) - rho);

  double dU1 = -M1 / (Q * Q);
  double dU2 = -M2 / (Q * Q);
  double dU3 = -M3 / (Q * Q);

  double drho = -d * pow (1 + SQR (U1) + SQR (U2) + SQR(U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3);

  double dp = (Gamma_1 / Gamma) * (1 / sqrt (1 + SQR (U1) + SQR (U2) + SQR(U3))
	- Q * pow (1 + SQR (U1) + SQR (U2) + SQR(U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3) - drho);

  *f = Q * sqrt (1 + SQR (U1) + SQR (U2) + SQR(U3)) - pres - E;
  *df = sqrt (1 + SQR (U1) + SQR (U2) + SQR(U3)) + Q * (U1 * dU1 + U2 * dU2 + U3 * dU3) 
         / sqrt (1 + SQR (U1) + SQR (U2) + SQR(U3)) - dp;
}

*/

static double
FUN_Q (double Q, void *ptr)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  double d = (params->d);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
  double E = (params->E);

  double Msqr = SQR (M1) + SQR (M2) + SQR (M3);
  double Qsqr = SQR (Q);
  double f = (1.0 / (Gamma * sqrt (Msqr + Qsqr))) * (Gamma * Msqr + Qsqr + d * Q * Gamma_1) - E;

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

  double Msqr = SQR (M1) + SQR (M2) + SQR (M3);
  double Qsqr = SQR (Q);

  double a = (2 - Gamma) * Msqr * Q;
  double b = d * (Gamma_1) * Msqr;

  double df = (Qsqr * Q + a + b) / (Gamma * pow (Msqr + Qsqr, 1.50));

  return df;
}
static void
FDF_FUN_Q (double Q, void *ptr, double *f, double *df)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  double d = (params->d);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
  double E = (params->E);

  double Msqr = SQR (M1) + SQR (M2) + SQR (M3);
  double Qsqr = SQR (Q);
  *f = (1.0 / (Gamma * sqrt (Msqr + Qsqr))) * (Gamma * Msqr + Qsqr + d * Q * Gamma_1) - E;

  double a = (2 - Gamma) * Msqr * Q;
  double b = d * (Gamma_1) * Msqr;
  *df = (Qsqr * Q + a + b) / (Gamma * pow (Msqr + Qsqr, 1.50));
}
