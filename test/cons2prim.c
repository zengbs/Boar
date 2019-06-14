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
  Gamma = 1.4;
  Gamma_1 = Gamma - 1.0;
  Gamma_2 = Gamma - 2.0;

  PrimS prim;
  ConsS cons;
/*
 double M1 =2.5;
 double M2 =3.0;
 double M3 =0.0;
 double E = 21;

 double Msqr=SQR(M1)+SQR(M2)+SQR(M3);
 double d = -200+sqrt(Msqr)/Gamma_1;
*/
  cons.d =  5.918569e-01;
  cons.M1 = -1.750500e+07;
  cons.M2 = -1.750500e+07;
  cons.M3 = -1.570834e+08;
  cons.E =  1.535905e+08;
  double Msqr=SQR(cons.M1)+SQR(cons.M2)+SQR(cons.M3);

  printf ("Initial cons:\n\n");
  printf ("D=%E, M1=%E, M2=%E M3=%E, E=%E\n\n",
          cons.d, cons.M1, cons.M2, cons.M3, cons.E);
  printf ("===========================================\n\n");

  prim = Cons_to_Prim (&cons);
  double d = prim.d;
  double U1 = prim.U1;
  double U2 = prim.U2;
  double U3 = prim.U3;
  double P = prim.P;

  printf ("transform to prim:\n\n");
  printf ("d=%E, U1=%E, U2=%E U3=%E, P=%E\n\n", d, U1, U2, U3, P);
  printf ("===========================================\n\n");

  ConsS cons2 = Prim_to_Cons (&prim);

  printf ("transform back to cons:\n\n");
  printf ("D=%E, M1=%E, M2=%E M3=%E, E=%E\n\n", cons2.d, cons2.M1, cons2.M2, cons2.M3, cons2.E);
  printf ("===========================================\n\n");

  if ((fabs (cons2.M1) > TINY_NUMBER) 
   && (fabs (cons2.M2) > TINY_NUMBER)
   && (fabs (cons2.M3) > TINY_NUMBER))
    {
      double err_d = (cons2.d - cons.d) / cons.d;
      double err_M1 = (cons2.M1 - cons.M1) / cons.M1;
      double err_M2 = (cons2.M2 - cons.M2) / cons.M2;
      double err_M3 = (cons2.M3 - cons.M3) / cons.M3;
      double err_E = (cons2.E - cons.E) / cons.E;

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_M2=%E err_M3=%E, err_P=%E\n", err_d,  err_M1, err_M2, err_M3, err_E);
      printf ("===========================================\n\n");

    }
  else if ((fabs (cons2.M1) < TINY_NUMBER) 
	&& (fabs (cons2.M2) < TINY_NUMBER)
	&& (fabs (cons2.M3) > TINY_NUMBER))
    {
      double err_d = (cons2.d - cons.d) / cons.d;
      double err_E = (cons2.E - cons.E) / cons.E;
      double err_M3 = (cons2.M3 - cons.M3) / cons.M3;

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M3=%E, err_P=%E\n", err_d, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (cons2.M1) < TINY_NUMBER) 
	&& (fabs (cons2.M2) > TINY_NUMBER)
	&& (fabs (cons2.M3) < TINY_NUMBER))
    {
      double err_d = (cons2.d - cons.d) / cons.d;
      double err_E = (cons2.E - cons.E) / cons.E;
      double err_M2 = (cons2.M2 - cons.M2) / cons.M2;

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M2=%E, err_P=%E\n", err_d, err_M2, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (cons2.M1) > TINY_NUMBER) 
	&& (fabs (cons2.M2) < TINY_NUMBER)
	&& (fabs (cons2.M3) < TINY_NUMBER))
    {
      double err_d = (cons2.d - cons.d) / cons.d;
      double err_E = (cons2.E - cons.E) / cons.E;
      double err_M1 = (cons2.M1 - cons.M1) / cons.M1;

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_P=%E\n", err_d, err_M1, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (cons2.M1) > TINY_NUMBER) 
	&& (fabs (cons2.M2) > TINY_NUMBER)
	&& (fabs (cons2.M3) < TINY_NUMBER))
    {
      double err_d = (cons2.d - cons.d) / cons.d;
      double err_E = (cons2.E - cons.E) / cons.E;
      double err_M1 = (cons2.M1 - cons.M1) / cons.M1;
      double err_M2 = (cons2.M2 - cons.M2) / cons.M2;

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_M2=%E, err_P=%E\n", err_d, err_M1, err_M2, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (cons2.M1) > TINY_NUMBER) 
	&& (fabs (cons2.M2) < TINY_NUMBER)
	&& (fabs (cons2.M3) > TINY_NUMBER))
    {
      double err_d = (cons2.d - cons.d) / cons.d;
      double err_E = (cons2.E - cons.E) / cons.E;
      double err_M1 = (cons2.M1 - cons.M1) / cons.M1;
      double err_M3 = (cons2.M3 - cons.M3) / cons.M3;

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_M3=%E, err_P=%E\n", err_d, err_M1, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (cons2.M1) < TINY_NUMBER) 
	&& (fabs (cons2.M2) > TINY_NUMBER)
	&& (fabs (cons2.M3) > TINY_NUMBER))
    {
      double err_d = (cons2.d - cons.d) / cons.d;
      double err_E = (cons2.E - cons.E) / cons.E;
      double err_M2 = (cons2.M2 - cons.M2) / cons.M2;
      double err_M3 = (cons2.M3 - cons.M3) / cons.M3;

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M2=%E, err_M3=%E, err_P=%E\n", err_d, err_M2, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (cons2.M1) < TINY_NUMBER) 
	&& (fabs (cons2.M2) < TINY_NUMBER)
	&& (fabs (cons2.M3) < TINY_NUMBER))
    {
      double err_d = (cons2.d - cons.d) / cons.d;
      double err_E = (cons2.E - cons.E) / cons.E;

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_P=%E\n", err_d, err_E);
      printf ("===========================================\n\n");
    }

  return 0;
}
