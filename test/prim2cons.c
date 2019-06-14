#define MAIN_C
#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "struct.h"
#include "globals.h"
#include "prototypes.h"

int
main ()
{
  Gamma = 1.66666;
  Gamma_1 = Gamma-1.0;
  Gamma_2 = Gamma-2.0;

  ConsS cons;
  PrimS prim;

  prim.d = 1.0e-130;
  prim.U1 = 1.0e-130;
  prim.U2 = 0.0;
  prim.U3 = 0.00;
  prim.P = 1.0e-130;

  printf ("Initial prim:\n\n");
  printf ("rho=%E, U1=%E, U2=%E U3=%E, P=%E\n\n", prim.d, prim.U1, prim.U2, prim.U3, prim.P);
  printf ("===========================================\n\n");

  cons = Prim_to_Cons (&prim);

  double d = cons.d;
  double M1 = cons.M1;
  double M2 = cons.M2;
  double M3 = cons.M3;
  double E = cons.E;

  printf ("after transform to cons:\n\n");
  printf ("d=%E, M1=%E, M2=%E M3=%E, E=%E\n\n", d, M1, M2, M3, E);
  printf ("===========================================\n\n");

  PrimS prim2 = Cons_to_Prim (&cons);

  if ((fabs (M1) > TINY_NUMBER) 
   && (fabs (M2) > TINY_NUMBER)
   && (fabs (M3) > TINY_NUMBER))
    {
      double err_d = (prim2.d - prim.d) / prim.d;
      double err_U1 = (prim2.U1 - prim.U1) / prim.U1;
      double err_U2 = (prim2.U2 - prim.U2) / prim.U2;
      double err_U3 = (prim2.U3 - prim.U3) / prim.U3;
      double err_P = (prim2.P - prim.P) / prim.P;

      printf ("recover to prim again and calculate relative error:\n\n");
      printf ("err_d=%E, err_U1=%E, err_U2=%E err_U3=%E, err_P=%E\n", err_d,
	      err_U1, err_U2, err_U3, err_P);
      printf ("===========================================\n\n");

    }
  else
    {
      printf ("recover to prim again:\n\n");
      printf ("rho=%E, U1=%E, U2=%E U3=%E, P=%E\n\n", prim2.d, prim2.U1,
	      prim2.U2, prim2.U3, prim2.P);
      printf ("===========================================\n\n");

    }

  return 0;
}
