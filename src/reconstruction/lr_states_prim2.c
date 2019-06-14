/*============================================================================*/
/*! \file lr_states_prim2.c
 *  \brief Second order (piecewise linear) spatial reconstruction in the
 *   primitive variables. 
 *
 * PURPOSE: Second order (piecewise linear) spatial reconstruction in the
 *   primitive variables. With the CTU integrator, a time-evolution
 *   (characteristic tracing) step is used to interpolate interface values
 *   to the half time level {n+1/2}.
 *
 *   Limiting is performed in the primitive (rather than characteristic)
 *   variables.  When used with the VL integrator, an eigenvalue decomposition
 *   is NOT needed.
 *
 * NOTATION: 
 * - W_{L,i-1/2} is reconstructed value on the left-side of interface at i-1/2
 * - W_{R,i-1/2} is reconstructed value on the right-side of interface at i-1/2
 *
 *   The L- and R-states at the left-interface in each cell are indexed i.
 * - W_{L,i-1/2} is denoted by Wl[i  ];   W_{R,i-1/2} is denoted by Wr[i  ]
 * - W_{L,i+1/2} is denoted by Wl[i+1];   W_{R,i+1/2} is denoted by Wr[i+1]
 *
 *   Internally, in this routine, Wlv and Wrv are the reconstructed values on
 *   the left-and right-side of cell center.  Thus (see Step 8),
 * -   W_{L,i-1/2} = Wrv(i-1);  W_{R,i-1/2} = Wlv(i)
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - lr_states()          - computes L/R states
 * - lr_states_init()     - initializes memory for static global arrays
 * - lr_states_destruct() - frees memory for static global arrays	      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../struct.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

static double **pW=NULL;
/*static double **vel=NULL;*/

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states(const GridS *pG, const Prim1DS W[], const double Bxc[],
 *               const double dt, const double dx, const int il, const int iu,
 *               Prim1DS Wl[], Prim1DS Wr[], const int dir)
 *  \brief Computes L/R states
 * Input Arguments:
 * - W = PRIMITIVE variables at cell centers along 1-D slice
 * - Bxc = B in direction of slice at cell center
 * - dtodx = dt/dx
 * - il,iu = lower and upper indices of zone centers in slice
 * W and Bxc must be initialized over [il-2:iu+2]
 *
 * Output Arguments:
 * - Wl,Wr = L/R-states of PRIMITIVE variables at interfaces over [il:iu+1]
 */

void lr_states(const GridS *pG __attribute((unused)),
               const Prim1DS W[], const double Bxc[],
               const double dt, const double dx, const int il, const int iu,
               Prim1DS Wl[], Prim1DS Wr[], 
               const int dir __attribute__((unused)))
{
  int i,n;
  double lim_slope1,lim_slope2;
  double dWc[5+0],dWl[5+0];
  double dWr[5+0],dWg[5+0];
  double Wlv[5+0],Wrv[5+0];
  double dWm[5+0];
  double *pWl, *pWr;
  double minmod_coeff =2.0;
/* Set pointer to primitive variables */

  for (i=il-2; i<=iu+2; i++) pW[i] = (double*)&(W[i]);
/*========================== START BIG LOOP OVER i =======================*/
  for (i=il-1; i<=iu+1; i++) {

/*--- Step 1. ------------------------------------------------------------------
 * Compute centered, L/R, and van Leer differences of primitive variables
 * Note we access contiguous array elements by indexing pointers for speed */
/* W1d = (d, V1, V2, V3, P, B2c, B3c, s[n]) */

    for (n=0; n<(5+0); n++) {
      dWc[n] = pW[i+1][n] - pW[i-1][n];
      dWl[n] = pW[i][n]   - pW[i-1][n];
      dWr[n] = pW[i+1][n] - pW[i][n];
      if (dWl[n]*dWr[n] > 0.0) {
        dWg[n] = 2.0*dWl[n]*dWr[n]/(dWl[n]+dWr[n]);
      } else {
        dWg[n] = 0.0;
      }
    }

/*--- Step 2. ------------------------------------------------------------------
 * Apply monotonicity constraints to differences in primitive vars. */
    for (n=0; n<(5+0); n++) {
      dWm[n] = 0.0;
      if (dWl[n]*dWr[n] > 0.0) {
        lim_slope1 = MIN(    fabs(dWl[n]),fabs(dWr[n]));
        lim_slope2 = MIN(0.5*fabs(dWc[n]),fabs(dWg[n]));
        dWm[n] = SIGN(dWc[n])*MIN(minmod_coeff*lim_slope1,lim_slope2);
      }
    }

/*--- Step 3. ------------------------------------------------------------------
 * Compute L/R values, ensure they lie between neighboring cell-centered vals */

    for (n=0; n<(5+0); n++) {
      Wlv[n] = pW[i][n] - 0.5*dWm[n];
      Wrv[n] = pW[i][n] + 0.5*dWm[n];
    }

    for (n=0; n<(5+0); n++) {
      Wlv[n] = MAX(MIN(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wlv[n] = MIN(MAX(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wrv[n] = MAX(MIN(pW[i][n],pW[i+1][n]),Wrv[n]);
      Wrv[n] = MIN(MAX(pW[i][n],pW[i+1][n]),Wrv[n]);
    }

/*--- Step 4. ------------------------------------------------------------------
 * Set L/R values */

    pWl = (double *) &(Wl[i+1]);
    pWr = (double *) &(Wr[i]);

    for (n=0; n<(5+0); n++) {
      pWl[n] = Wrv[n];
      pWr[n] = Wlv[n];
    }

  } /*===================== END BIG LOOP OVER i ===========================*/

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states_init(MeshS *pM)
 *  \brief Allocate enough memory for work arrays */

void lr_states_init(DomainS *pD)
{
  int nmax,size1=0,size2=0,size3=0;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
        if (pD->Grid->Nx[0] > size1){
          size1 = pD->Grid->Nx[0];
        }
        if (pD->Grid->Nx[1] > size2){
          size2 = pD->Grid->Nx[1];
        }
        if (pD->Grid->Nx[2] > size3){
          size3 = pD->Grid->Nx[2];
        }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

  if ((pW = (double**)malloc(nmax*sizeof(double*))) == NULL) goto on_error;
  return;
  on_error:
    lr_states_destruct();
    printf ("[lr_states_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states_destruct(void)
 *  \brief Free memory used by work arrays */

void lr_states_destruct(void)
{
  if (pW != NULL) free(pW);
  //if (vel != NULL) free_2d_array(vel);
  return;
}
