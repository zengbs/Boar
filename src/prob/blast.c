/*============================================================================*/
/*! \file blast.c
 *  \brief Problem generator for spherical blast wave problem.
 *
 * PURPOSE: Problem generator for spherical blast wave problem.
 *
 * REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
 *   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.     */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "struct.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  Prim1DS W;
  Cons1DS U1d;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  double drat,prat,rad,pa,da,x1,x2,x3,xc1,xc2,xc3;
  double rin;
  int centerX,centerY,centerZ; /*blast center*/

  rin = 0.01; /*Radius of the inner sphere*/
  pa  = 1.0e-3; /*ambient pressure*/
  da  = 0.8; /*ambient density*/
  drat = 2.0; /*density ratio*/
  prat = 1.0e+10; /*Pressure ratio initially*/

/* setup uniform ambient medium with spherical over-pressured region */

  W.d = da;
  W.Ux = 0.0;
  W.Uy = 0.0;
  W.Uz = 0.0;

       centerX = ceil((double)(pDomain->Grid->Nx[0])/2.0);
       centerY = ceil((double)(pDomain->Grid->Nx[1])/2.0);
       centerZ = ceil((double)(pDomain->Grid->Nx[2])/2.0);
       cc_pos(pGrid,centerX,centerY,centerZ,&xc1,&xc2,&xc3);

       printf("%d %d %d\n", centerX, centerY, centerZ);
       printf("%f %f %f\n", xc1, xc2, xc3);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	rad = sqrt(SQR(x1-xc1) + SQR(x2-xc2) + SQR(x3-xc3));
        W.P = pa;
	if (rad < rin) W.P = prat*pa;
        W.d = da;
	if (rad < rin) W.d = drat*da;

        U1d = Prim1D_to_Cons1D(&(W));

	pGrid->U[k][j][i].d  = U1d.d;
	pGrid->U[k][j][i].M1 = U1d.Mx;
	pGrid->U[k][j][i].M2 = U1d.My;
	pGrid->U[k][j][i].M3 = U1d.Mz;
	pGrid->U[k][j][i].E  = U1d.E;
      }
    }
  }
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/
/*
void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}
*/
