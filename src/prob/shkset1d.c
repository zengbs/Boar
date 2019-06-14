/*============================================================================*/
/*! \file shkset1d.c 
 *  \brief Problem generator for 1-D Riemann problems.  
 *
 * PURPOSE: Problem generator for 1-D Riemann problems.  Initial discontinuity
 *   is located so there are equal numbers of cells to the left and right (at
 *   center of grid based on integer index).  Initializes plane-parallel shock
 *   along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
 *
 * If error_test=1 in the <problem> block, then the L1 error in the final
 * solution will be computed for the Sod shocktube (hydrodynamics) or the RJ2a
 * test (MHD).  This is useful for regression tests.
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "globals.h"
#include "prototypes.h"

void
problem (DomainS * pDomain)
{
  GridS *pGrid = (pDomain->Grid);
  int i, il, iu, j, jl, ju, k, kl, ku;
//  int is,ie,js,je,ks,ke;

  int shk_dir;			/* Shock direction: {1,2,3} -> {x1,x2,x3} */

  double lVx, rVx, lVy, rVy, lVz, rVz;
  double factor;
  double x1, x2, x3;
  Prim1DS Wl, Wr;
  Cons1DS U1d, Ul, Ur;
/*
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
*/
/* Parse left state read from input file: dl,pl,ul,vl,wl,bxl,byl,bzl */

  Wl.d = 0.9;
  Wl.P = 5.0e+10;
  lVx = 0.1;
  lVy = 0.0;
  lVz = 0.0;

  if (SQR (lVx) + SQR (lVy) + SQR (lVz) >= 1.0)
    {
      printf ("\n|V| >= speed of light!\n");
      abort ();
    }

  factor = 1.0 / sqrt (1.0 - SQR (lVx) - SQR (lVy) - SQR (lVz));

  Wl.Ux = factor * lVx;
  Wl.Uy = factor * lVy;
  Wl.Uz = factor * lVz;

/* Parse right state read from input file: dr,pr,ur,vr,wr,bxr,byr,bzr */

  Wr.d = 1.0;
  Wr.P = 1e-5;
  rVx = 0.1;
  rVy = 0.0;
  rVz = 0.0;

  if (SQR (rVx) + SQR (rVy) + SQR (rVz) >= 1.0)
    {
      printf ("\n|V| >= speed of light!\n");
      abort ();
    }

  factor = 1.0 / sqrt (1.0 - SQR (rVx) - SQR (rVy) - SQR (rVz));

  Wr.Ux = factor * rVx;
  Wr.Uy = factor * rVy;
  Wr.Uz = factor * rVz;

  Ul = Prim1D_to_Cons1D (&Wl);
  Ur = Prim1D_to_Cons1D (&Wr);

/* Parse shock direction */
  shk_dir = 1;			/*1, 2 or 3 */

/* Set up the index bounds for initializing the grid */
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1)
    {
      ju = pGrid->je + nghost;
      jl = pGrid->js - nghost;
    }
  else
    {
      ju = pGrid->je;
      jl = pGrid->js;
    }

  if (pGrid->Nx[2] > 1)
    {
      ku = pGrid->ke + nghost;
      kl = pGrid->ks - nghost;
    }
  else
    {
      ku = pGrid->ke;
      kl = pGrid->ks;
    }

/* Initialize the grid including the ghost cells.  Discontinuity is always
 * located at x=0, so xmin/xmax in input file must be set appropriately. */

  switch (shk_dir)
    {
/*--- shock in 1-direction ---------------------------------------------------*/
    case 1:			/* shock in 1-direction  */

      for (k = kl; k <= ku; k++)
	{
	  for (j = jl; j <= ju; j++)
	    {
	      for (i = il; i <= iu; i++)
		{
		  cc_pos (pGrid, i, j, k, &x1, &x2, &x3);

/* set primitive and conserved variables to be L or R state */
		  if (x1 <= 0.5 * (pGrid->MaxX[0] + pGrid->MinX[0]))
		    {
		      U1d = Ul;
		    }
		  else
		    {
		      U1d = Ur;
		    }

/* Initialize conserved (and with SR the primitive) variables in Grid */
		  pGrid->U[k][j][i].d = U1d.d;
		  pGrid->U[k][j][i].M1 = U1d.Mx;
		  pGrid->U[k][j][i].M2 = U1d.My;
		  pGrid->U[k][j][i].M3 = U1d.Mz;
		  pGrid->U[k][j][i].E = U1d.E;
		}
	    }
	}
      break;

/*--- shock in 2-direction ---------------------------------------------------*/
    case 2:			/* shock in 2-direction  */
      for (k = kl; k <= ku; k++)
	{
	  for (j = jl; j <= ju; j++)
	    {
	      for (i = il; i <= iu; i++)
		{
		  cc_pos (pGrid, i, j, k, &x1, &x2, &x3);

/* set primitive variables to be L or R state */
		  if (x2 <= 0.5 * (pGrid->MaxX[0] + pGrid->MinX[0]))
		    {
		      U1d = Ul;
		    }
		  else
		    {
		      U1d = Ur;
		    }

/* Initialize conserved (and with SR the primitive) variables in Grid */
		  pGrid->U[k][j][i].d = U1d.d;
		  pGrid->U[k][j][i].M1 = -U1d.My;
		  pGrid->U[k][j][i].M2 = U1d.Mx;
		  pGrid->U[k][j][i].M3 = U1d.Mz;
		  pGrid->U[k][j][i].E = U1d.E;
		}
	    }
	}
      break;

/*--- shock in 3-direction ---------------------------------------------------*/
    case 3:			/* shock in 3-direction  */
      for (k = kl; k <= ku; k++)
	{
	  for (j = jl; j <= ju; j++)
	    {
	      for (i = il; i <= iu; i++)
		{
		  cc_pos (pGrid, i, j, k, &x1, &x2, &x3);

/* set primitive variables to be L or R state */
		  if (x3 <= 0.5 * (pGrid->MaxX[0] + pGrid->MinX[0]))
		    {
		      U1d = Ul;
		    }
		  else
		    {
		      U1d = Ur;
		    }

/* Initialize conserved (and with SR the primitive) variables in Grid */
		  pGrid->U[k][j][i].d = U1d.d;
		  pGrid->U[k][j][i].M1 = -U1d.Mz;
		  pGrid->U[k][j][i].M2 = U1d.My;
		  pGrid->U[k][j][i].M3 = U1d.Mx;
		  pGrid->U[k][j][i].E = U1d.E;
		}
	    }
	}
      break;
    }
  return;
}
