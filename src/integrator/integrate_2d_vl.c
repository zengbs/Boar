/*============================================================================*/
/*! \file integrate_2d_vl_sr.c
 *  \brief Integrate MHD equations using 2D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator.
 *
 * PURPOSE: Integrate MHD equations using 2D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator.  The variables updated are:
 *   -  U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *   -  B1i, B2i  -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, and the H-correction
 *   of Sanders et al.
 *
 * REFERENCE: 
 * - J.M Stone & T.A. Gardiner, "A simple, unsplit Godunov method
 *   for multidimensional MHD", NewA 14, 139 (2009)
 *
 * - R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensional dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_2d_vl()
 * - integrate_destruct_2d()
 * - integrate_init_2d() */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../struct.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

/* The L/R states of primitive variables and fluxes at each cell face */
static Prim1DS **Wl_x1Face=NULL, **Wr_x1Face=NULL;
static Prim1DS **Wl_x2Face=NULL, **Wr_x2Face=NULL;
static Cons1DS **x1Flux=NULL, **x2Flux=NULL;
static Cons1DS **x1FluxP=NULL, **x2FluxP=NULL; /* first order flux */

/* 1D scratch vectors used by lr_states and flux functions */
static double *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W1d=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* primitive variables at t^{n} */
static PrimS **W=NULL;

/* conserved and primitive variables at t^{n+1/2} computed in predict step */
static ConsS **Uhalf=NULL;
static PrimS **Whalf=NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf3_corner() - the upwind CT method of Gardiner & Stone (2005) 
 *   FixCell() - apply first-order correction to one cell
 *============================================================================*/
static void FixCell(GridS *pG, Int3Vect);
/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_2d_vl(DomainS *pD)
 *  \brief van Leer unsplit integrator in 2D. 
 *
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 2D.
 */
void integrate_2d_vl(DomainS *pD)
{
  GridS *pG = pD->Grid;
  pG->dt = pD->dt;
  pG->time = pD->time;
  double dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2;
  double hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  int cart_x1 = 1, cart_x2 = 2;//, cart_x3 = 3;
  int flag_cell=0,negd=0,negP=0,NaNFlux=0;
//  ConsS Ucheck;
  PrimS Wcheck;
  Int3Vect BadCell;
  int il=is-(nghost-1), iu=ie+(nghost-1);
  int jl=js-(nghost-1), ju=je+(nghost-1);

/* Uhalf = conserved variables at t = n
 * W     = primitive variables at t = n 
 * */
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      Uhalf[j][i] = pG->U[ks][j][i];
      W[j][i] = Cons_to_Prim(&(pG->U[ks][j][i]));
      if ((W[j][i].d<0.0)||(W[j][i].P)<0.0)
      {
        printf("W[%i][%i].d=%f, W[%i][%i].P=%f\n", j,i,W[j][i].d,j,i,W[j][i].P);
      }
    }
  }
/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * W1d = (d, U1, U2, U3, P, B2c, B3c, s[n])
 */

  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      W1d[i].d  = W[j][i].d;
      W1d[i].Ux = W[j][i].U1;
      W1d[i].Uy = W[j][i].U2;
      W1d[i].Uz = W[j][i].U3;
      W1d[i].P  = W[j][i].P;
    }

/*--- Step 1b ------------------------------------------------------------------
 * Compute first-order L/R states */
    for (i=il; i<=ie+nghost; i++) {
      Wl[i] = W1d[i-1];
      Wr[i] = W1d[i  ];

/* Compute U from W in case Pfloor used in Cons1D_to_Prim1D */
      Ul[i] = Prim1D_to_Cons1D(&Wl[i]);
      Ur[i] = Prim1D_to_Cons1D(&Wr[i]);
    }

/*--- Step 1c ------------------------------------------------------------------
 * No source terms needed */

/*--- Step 1d ------------------------------------------------------------------
 * Compute flux in x1-direction */

    for (i=il; i<=ie+nghost; i++) {
      fluxes(Ul[i],Ur[i],Wl[i],Wr[i],&x1Flux[j][i]);
    }
  }

/*=== STEP 2: Compute first-order fluxes at t^{n} in x2-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, U2, U3, U1, P, B3c, B1c, s[n])
 */

  for (i=is-nghost; i<=ie+nghost; i++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      W1d[j].d  = W[j][i].d;
      W1d[j].Ux = W[j][i].U2;
      W1d[j].Uy = W[j][i].U3;
      W1d[j].Uz = W[j][i].U1;
      W1d[j].P  = W[j][i].P;
    }

/*--- Step 2b ------------------------------------------------------------------
 * Compute first-order L/R states */
    for (j=jl; j<=je+nghost; j++) {
      Wl[j] = W1d[j-1];
      Wr[j] = W1d[j  ];

/* Compute U from W in case Pfloor used in Cons1D_to_Prim1D */
      Ul[j] = Prim1D_to_Cons1D(&Wl[j]);
      Ur[j] = Prim1D_to_Cons1D(&Wr[j]);
    }

/*--- Step 2c ------------------------------------------------------------------
 * No source terms needed */

/*--- Step 2d ------------------------------------------------------------------
 * Compute flux in x2-direction */

    for (j=jl; j<=je+nghost; j++) {
      fluxes(Ul[j],Ur[j],Wl[j],Wr[j],&x2Flux[j][i]);
    }
  }

/*=== STEP 3: Not needed in 2D ===*/

/*=== STEP 4:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 4a ------------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at t^{n} and integrate
 * to corner.
 */


/*=== STEP 5: Update cell-centered variables to half-timestep ================*/

/*--- Step 5a ------------------------------------------------------------------
 * Update cell-centered variables (including B3c) to half-timestep with x1Flux
 */

  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      Uhalf[j][i].d   -= hdtodx1*(x1Flux[j][i+1].d  - x1Flux[j][i].d );
      Uhalf[j][i].M1  -= hdtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx);
      Uhalf[j][i].M2  -= hdtodx1*(x1Flux[j][i+1].My - x1Flux[j][i].My);
      Uhalf[j][i].M3  -= hdtodx1*(x1Flux[j][i+1].Mz - x1Flux[j][i].Mz);
      Uhalf[j][i].E   -= hdtodx1*(x1Flux[j][i+1].E  - x1Flux[j][i].E );
    }
  }

/*--- Step 5b ------------------------------------------------------------------
 * Update cell-centered variables (including B3c) to half-timestep with x2Flux
 */

  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      Uhalf[j][i].d   -= hdtodx2*(x2Flux[j+1][i].d  - x2Flux[j][i].d );
      Uhalf[j][i].M1  -= hdtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);
      Uhalf[j][i].M2  -= hdtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);
      Uhalf[j][i].M3  -= hdtodx2*(x2Flux[j+1][i].My - x2Flux[j][i].My);
      Uhalf[j][i].E   -= hdtodx2*(x2Flux[j+1][i].E  - x2Flux[j][i].E );
    }
  }

/*--- Step 5d ------------------------------------------------------------------
 * With first-order flux correction, check predict fluxes & emf3 for NaN's
 * and save predict fluxes & emf3 for first-order flux correction
 */

  NaNFlux = 0;
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      if ((x1Flux[j][i].d  != x1Flux[j][i].d)  ||
          (x1Flux[j][i].E  != x1Flux[j][i].E)  ||
          (x1Flux[j][i].Mx != x1Flux[j][i].Mx) ||
          (x1Flux[j][i].My != x1Flux[j][i].My) ||
          (x1Flux[j][i].Mz != x1Flux[j][i].Mz)) {
        NaNFlux++;
      }          
      if ((x2Flux[j][i].d  != x2Flux[j][i].d)  ||
          (x2Flux[j][i].E  != x2Flux[j][i].E)  ||
          (x2Flux[j][i].Mx != x2Flux[j][i].Mx) ||
          (x2Flux[j][i].My != x2Flux[j][i].My) ||
          (x2Flux[j][i].Mz != x2Flux[j][i].Mz)) {
        NaNFlux++;
      }
      x1FluxP[j][i] = x1Flux[j][i];
      x2FluxP[j][i] = x2Flux[j][i];
    }
  }
  if (NaNFlux != 0) {
    printf ("[Step5d ]: %3i first-order fluxes are NaN! nstep=%5i\n",NaNFlux, pD->nstep);
    exit(EXIT_FAILURE);
  }
/*=== STEP 7: Conserved->Primitive variable inversion at t^{n+1/2} ===========*/
        
/* Invert conserved variables at t^{n+1/2} to primitive variables. With FOFC, 
 * check if cell-centered d < 0, P< 0, or v^2 > 1. With Entropy fix, correct
 * by computing new primitive state using the entropy variable, otherwise
 * correct by switching back to values at beginning of step, rendering update
 * first order in time for that cell.
 */

  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      Whalf[j][i] = Cons_to_Prim(&(Uhalf[j][i]));
    }
  }

  negd = 0; /* # of cell that posses minus density */
  negP = 0; /* # of cell that posses minus pressure */
  flag_cell = 0;
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      Whalf[j][i] = Cons_to_Prim(&(Uhalf[j][i]));
      if (Whalf[j][i].d < 0.0) {
        flag_cell = 1;
        negd++;
      }
      if (Whalf[j][i].P < 0.0) {
        flag_cell = 1;
        negP++;
      }
/* Replace W^{n+0.5} with W^{n}
 * */
      if (flag_cell != 0) {
          printf("Uhalf[%i][%i].d  =%f\n",j,i,Uhalf[j][i].d);
          printf("Uhalf[%i][%i].M1=%f\n",j,i,Uhalf[j][i].M1);
          printf("Uhalf[%i][%i].M2=%f\n",j,i,Uhalf[j][i].M2);
          printf("Uhalf[%i][%i].M3=%f\n",j,i,Uhalf[j][i].M3);
          printf("Uhalf[%i][%i].E  =%f\n",j,i,Uhalf[j][i].E);
          printf("Whalf[%i][%i].d  =%f\n",j,i,Whalf[j][i].d);
          printf("Whalf[%i][%i].U1  =%f\n",j,i,Whalf[j][i].U1);
          printf("Whalf[%i][%i].U2  =%f\n",j,i,Whalf[j][i].U2);
          printf("Whalf[%i][%i].U3  =%f\n",j,i,Whalf[j][i].U3);
          printf("Whalf[%i][%i].P  =%f\n",j,i,Whalf[j][i].P);
	  Whalf[j][i].d = W[j][i].d;
	  Whalf[j][i].P = W[j][i].P;
	  flag_cell=0;
      }
    }
  }

  if (negd > 0 || negP > 0 ){
    printf("[Step7  ]: %3i cells had d<0; %3i cells had P<0; nstep=%5i\n",negd,negP,pD->nstep);
  }
/*=== STEP 8: Compute second-order L/R x1-interface states ===================*/

/*--- Step 8a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, U1, U2, U3, P, B2c, B3c, s[n])
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=il; i<=iu; i++) {
      W1d[i].d  = Whalf[j][i].d;
      W1d[i].Ux = Whalf[j][i].U1;
      W1d[i].Uy = Whalf[j][i].U2;
      W1d[i].Uz = Whalf[j][i].U3;
      W1d[i].P  = Whalf[j][i].P;
    }

/*--- Step 8b ------------------------------------------------------------------
 * Compute L/R states on x1-interfaces, store into arrays
 */

    lr_states(pG,W1d,Bxc,pG->dt,pG->dx1,is,ie,Wl,Wr,cart_x1);

    for (i=is; i<=ie+1; i++) {
      Wl_x1Face[j][i] = Wl[i];
      Wr_x1Face[j][i] = Wr[i];

       if ((Wr[i].P < 0.0) || (Wr[i].d < 0.0) || (Wl[i].P < 0.0) || (Wl[i].d < 0.0))
        {
          printf("Wr[%d].d=%f, Wr[%d].P=%f\n", i, Wr[i].d, i, Wr[i].P);
          printf("Wl[%d].d=%f, Wl[%d].P=%f\n", i, Wl[i].d, i, Wl[i].P);
          abort();
        }

    }
  }


/*=== STEP 9: Compute second-order L/R x2-interface states ===================*/

/*--- Step 9a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, U2, U3, U1, P, B3c, B1c, s[n])
 */

  for (i=is-1; i<=ie+1; i++) {
    for (j=jl; j<=ju; j++) {
      W1d[j].d  = Whalf[j][i].d;
      W1d[j].Ux = Whalf[j][i].U2;
      W1d[j].Uy = Whalf[j][i].U3;
      W1d[j].Uz = Whalf[j][i].U1;
      W1d[j].P  = Whalf[j][i].P;
    }

/*--- Step 9b ------------------------------------------------------------------
 * Compute L/R states on x2-interfaces, store into arrays
 */

    lr_states(pG,W1d,Bxc,pG->dt,pG->dx2,js,je,Wl,Wr,cart_x2);

    for (j=js; j<=je+1; j++) {
      Wl_x2Face[j][i] = Wl[j];
      Wr_x2Face[j][i] = Wr[j];

       if ((Wr[i].P < 0.0) || (Wr[i].d < 0.0) || (Wl[i].P < 0.0) || (Wl[i].d < 0.0))
        {
          printf("Wr[%d].d=%f, Wr[%d].P=%f\n", i, Wr[i].d, i, Wr[i].P);
          printf("Wl[%d].d=%f, Wl[%d].P=%f\n", i, Wl[i].d, i, Wl[i].P);
          abort();
        }

    }
  }

/*=== STEP 10: Not needed in 2D ===*/

/*=== STEP 11: Compute 2D x1-Flux, x2-Flux, ==================================*/

/*--- Step 11b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      Ul[i] = Prim1D_to_Cons1D(&Wl_x1Face[j][i]);
      Ur[i] = Prim1D_to_Cons1D(&Wr_x1Face[j][i]);

      fluxes(Ul[i],Ur[i],Wl_x1Face[j][i],Wr_x1Face[j][i],&x1Flux[j][i]);

/* revert to predictor flux if this flux Nan'ed */
      if ((x1Flux[j][i].d  != x1Flux[j][i].d)  ||
          (x1Flux[j][i].E  != x1Flux[j][i].E)  ||
          (x1Flux[j][i].Mx != x1Flux[j][i].Mx) ||
          (x1Flux[j][i].My != x1Flux[j][i].My) ||
          (x1Flux[j][i].Mz != x1Flux[j][i].Mz)) {
        x1Flux[j][i] = x1FluxP[j][i];
        NaNFlux++;
      }
    }
  }

/*--- Step 11c -----------------------------------------------------------------
 * Compute second-order fluxes in x2-direction
 */

  for (j=js; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      Ul[i] = Prim1D_to_Cons1D(&Wl_x2Face[j][i]);
      Ur[i] = Prim1D_to_Cons1D(&Wr_x2Face[j][i]);

      fluxes(Ul[i],Ur[i],Wl_x2Face[j][i],Wr_x2Face[j][i],&x2Flux[j][i]);
/* revert to predictor flux if this flux NaN'ed */
      if ((x2Flux[j][i].d  != x2Flux[j][i].d)  ||
          (x2Flux[j][i].E  != x2Flux[j][i].E)  ||
          (x2Flux[j][i].Mx != x2Flux[j][i].Mx) ||
          (x2Flux[j][i].My != x2Flux[j][i].My) ||
          (x2Flux[j][i].Mz != x2Flux[j][i].Mz)) {
        x2Flux[j][i] = x2FluxP[j][i];
        NaNFlux++;
      }
    }
  }

  if (NaNFlux != 0) {
    printf("[Step11 ]: %3i second-order fluxes is(are) replaced; nstep=%5i\n",NaNFlux,pD->nstep);
    NaNFlux=0;
  }

/*=== STEP 14: Update cell-centered values for a full timestep ===============*/

/*--- Step 14a -----------------------------------------------------------------
 * Update cell-centered variables in pG (including B3c) using 2D x1-Fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].d   -= dtodx1*(x1Flux[j][i+1].d  - x1Flux[j][i].d );
      pG->U[ks][j][i].M1  -= dtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx);
      pG->U[ks][j][i].M2  -= dtodx1*(x1Flux[j][i+1].My - x1Flux[j][i].My);
      pG->U[ks][j][i].M3  -= dtodx1*(x1Flux[j][i+1].Mz - x1Flux[j][i].Mz);
      pG->U[ks][j][i].E   -= dtodx1*(x1Flux[j][i+1].E  - x1Flux[j][i].E );
    }
  }

/*--- Step 14b -----------------------------------------------------------------
 * Update cell-centered variables in pG (including B3c) using 2D x2-Fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].d   -= dtodx2*(x2Flux[j+1][i].d  - x2Flux[j][i].d );
      pG->U[ks][j][i].M1  -= dtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);
      pG->U[ks][j][i].M2  -= dtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);
      pG->U[ks][j][i].M3  -= dtodx2*(x2Flux[j+1][i].My - x2Flux[j][i].My);
      pG->U[ks][j][i].E   -= dtodx2*(x2Flux[j+1][i].E  - x2Flux[j][i].E );
    }
  }  

/*=== STEP 15: First-order flux correction ===================================*/

/*--- Step 15a -----------------------------------------------------------------
 * If cell-centered d or P have gone negative, or if v^2 > 1 in SR, correct
 * by using 1st order predictor fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      Wcheck = Cons_to_Prim(&(pG->U[ks][j][i]));
      if (Wcheck.d < 0.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = j;
        BadCell.k = ks;
        negd++;
        printf("b4: i=%i, j=%i, d=%f\n",i,j,Wcheck.d);
        printf("b4: i=%i, j=%i, P=%f\n",i,j,Wcheck.P);
      }
      if (Wcheck.P < 0.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = j;
        BadCell.k = ks;
        negP++;
        printf("b4: i=%i, j=%i, d=%f\n",i,j,Wcheck.d);
        printf("b4: i=%i, j=%i, P=%f\n",i,j,Wcheck.P);
        printf("d=%f,Mx=%f,My=%f,Mz=%f,E=%f\n",pG->U[ks][j][i].d,pG->U[ks][j][i].M1,pG->U[ks][j][i].M2,pG->U[ks][j][i].M3,pG->U[ks][j][i].E);
      }
      if (flag_cell != 0) {
        FixCell(pG, BadCell);
        flag_cell=0;
        Wcheck = Cons_to_Prim(&(pG->U[ks][j][i]));
        printf("i=%i, j=%i, d=%f\n",i,j,Wcheck.d);
        printf("i=%i, j=%i, P=%f\n",i,j,Wcheck.P);
      }
    }
  }

  if (negd > 0 || negP > 0){
    printf("[Step15a]: %3i cells had d<0; %3i cells had P<0; nstep=%5i\n",negd,negP,pD->nstep);
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_2d(MeshS *pM)
 *  \brief Allocate temporary integration arrays */
void integrate_init_2d(DomainS *pD)
{
  int nmax,size1=0,size2=0;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2 */
      if (pD->Grid != NULL) {
        if (pD->Grid->Nx[0] > size1){
          size1 = pD->Grid->Nx[0];
        }
        if (pD->Grid->Nx[1] > size2){
          size2 = pD->Grid->Nx[1];
        }
      }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  nmax = MAX(size1,size2);

  if ((Bxc = (double*)malloc(nmax*sizeof(double))) == NULL) goto on_error;
  if ((Bxi = (double*)malloc(nmax*sizeof(double))) == NULL) goto on_error;
  if ((U1d= (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W1d= (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Wl_x1Face = (Prim1DS**)calloc_2d_array(size2,size1, sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x1Face = (Prim1DS**)calloc_2d_array(size2,size1, sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wl_x2Face = (Prim1DS**)calloc_2d_array(size2,size1, sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x2Face = (Prim1DS**)calloc_2d_array(size2,size1, sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((x1Flux    = (Cons1DS**)calloc_2d_array(size2,size1, sizeof(Cons1DS))) 
    == NULL) goto on_error;
  if ((x2Flux    = (Cons1DS**)calloc_2d_array(size2,size1, sizeof(Cons1DS))) 
    == NULL) goto on_error;
  if ((x1FluxP = (Cons1DS**)calloc_2d_array(size2,size1, sizeof(Cons1DS))) 
    == NULL) goto on_error;
  if ((x2FluxP = (Cons1DS**)calloc_2d_array(size2,size1, sizeof(Cons1DS))) 
    == NULL) goto on_error;

  if ((Uhalf = (ConsS**)calloc_2d_array(size2,size1,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((Whalf = (PrimS**)calloc_2d_array(size2,size1,sizeof(PrimS)))==NULL)
    goto on_error;
  if ((W     = (PrimS**)calloc_2d_array(size2,size1,sizeof(PrimS)))==NULL)
    goto on_error;

  return;

  on_error:
//    integrate_destruct();
    printf("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_2d(void)
 *  \brief Free temporary integration arrays */
void integrate_destruct_2d(void)
{
  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
  if (U1d != NULL) free(U1d);
  if (Ul  != NULL) free(Ul);
  if (Ur  != NULL) free(Ur);
  if (W1d != NULL) free(W1d);
  if (Wl  != NULL) free(Wl);
  if (Wr  != NULL) free(Wr);

  if (Wl_x1Face != NULL) free_2d_array(Wl_x1Face);
  if (Wr_x1Face != NULL) free_2d_array(Wr_x1Face);
  if (Wl_x2Face != NULL) free_2d_array(Wl_x2Face);
  if (Wr_x2Face != NULL) free_2d_array(Wr_x2Face);
  if (x1Flux    != NULL) free_2d_array(x1Flux);
  if (x2Flux    != NULL) free_2d_array(x2Flux);
  if (x1FluxP   != NULL) free_2d_array(x1FluxP);
  if (x2FluxP   != NULL) free_2d_array(x2FluxP);
  if (Uhalf    != NULL) free_2d_array(Uhalf);
  if (Whalf    != NULL) free_2d_array(Whalf);
  if (W        != NULL) free_2d_array(W);

  return;
}
/*----------------------------------------------------------------------------*/
/*! \fn static void FixCell(GridS *pG, Int3Vect ix)
 *  \brief Uses first order fluxes to fix negative d,P or superluminal v
 */ 
static void FixCell(GridS *pG, Int3Vect ix)
{
  int ks=pG->ks;
  double dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2;
  Cons1DS x1FD_i, x1FD_ip1, x2FD_j, x2FD_jp1;

/* Compute difference of predictor and corrector fluxes at cell faces */

  x1FD_i.d = x1Flux[ix.j][ix.i].d - x1FluxP[ix.j][ix.i].d;
  x2FD_j.d = x2Flux[ix.j][ix.i].d - x2FluxP[ix.j][ix.i].d;

  x1FD_ip1.d = x1Flux[ix.j][ix.i+1].d - x1FluxP[ix.j][ix.i+1].d;
  x2FD_jp1.d = x2Flux[ix.j+1][ix.i].d - x2FluxP[ix.j+1][ix.i].d;

  x1FD_i.Mx = x1Flux[ix.j][ix.i].Mx - x1FluxP[ix.j][ix.i].Mx;
  x2FD_j.Mx = x2Flux[ix.j][ix.i].Mx - x2FluxP[ix.j][ix.i].Mx;

  x1FD_ip1.Mx = x1Flux[ix.j][ix.i+1].Mx - x1FluxP[ix.j][ix.i+1].Mx;
  x2FD_jp1.Mx = x2Flux[ix.j+1][ix.i].Mx - x2FluxP[ix.j+1][ix.i].Mx;

  x1FD_i.My = x1Flux[ix.j][ix.i].My - x1FluxP[ix.j][ix.i].My;
  x2FD_j.My = x2Flux[ix.j][ix.i].My - x2FluxP[ix.j][ix.i].My;

  x1FD_ip1.My = x1Flux[ix.j][ix.i+1].My - x1FluxP[ix.j][ix.i+1].My;
  x2FD_jp1.My = x2Flux[ix.j+1][ix.i].My - x2FluxP[ix.j+1][ix.i].My;

  x1FD_i.Mz = x1Flux[ix.j][ix.i].Mz - x1FluxP[ix.j][ix.i].Mz;
  x2FD_j.Mz = x2Flux[ix.j][ix.i].Mz - x2FluxP[ix.j][ix.i].Mz;

  x1FD_ip1.Mz = x1Flux[ix.j][ix.i+1].Mz - x1FluxP[ix.j][ix.i+1].Mz;
  x2FD_jp1.Mz = x2Flux[ix.j+1][ix.i].Mz - x2FluxP[ix.j+1][ix.i].Mz;

  x1FD_i.E = x1Flux[ix.j][ix.i].E - x1FluxP[ix.j][ix.i].E;
  x2FD_j.E = x2Flux[ix.j][ix.i].E - x2FluxP[ix.j][ix.i].E;

  x1FD_ip1.E = x1Flux[ix.j][ix.i+1].E - x1FluxP[ix.j][ix.i+1].E;
  x2FD_jp1.E = x2Flux[ix.j+1][ix.i].E - x2FluxP[ix.j+1][ix.i].E;

/* Use flux differences to correct bad cell-centered quantities */

  pG->U[ks][ix.j][ix.i].d  += dtodx1*(x1FD_ip1.d  - x1FD_i.d );
  pG->U[ks][ix.j][ix.i].M1 += dtodx1*(x1FD_ip1.Mx - x1FD_i.Mx);
  pG->U[ks][ix.j][ix.i].M2 += dtodx1*(x1FD_ip1.My - x1FD_i.My);
  pG->U[ks][ix.j][ix.i].M3 += dtodx1*(x1FD_ip1.Mz - x1FD_i.Mz);
  pG->U[ks][ix.j][ix.i].E  += dtodx1*(x1FD_ip1.E  - x1FD_i.E );

  pG->U[ks][ix.j][ix.i].d  += dtodx2*(x2FD_jp1.d  - x2FD_j.d );
  pG->U[ks][ix.j][ix.i].M1 += dtodx2*(x2FD_jp1.Mz - x2FD_j.Mz);
  pG->U[ks][ix.j][ix.i].M2 += dtodx2*(x2FD_jp1.Mx - x2FD_j.Mx);
  pG->U[ks][ix.j][ix.i].M3 += dtodx2*(x2FD_jp1.My - x2FD_j.My);
  pG->U[ks][ix.j][ix.i].E  += dtodx2*(x2FD_jp1.E  - x2FD_j.E );
/* Use flux differences to correct cell-centered values at i-1 and i+1 */

  if (ix.i > pG->is) {
    pG->U[ks][ix.j][ix.i-1].d  += dtodx1*(x1FD_i.d );
    pG->U[ks][ix.j][ix.i-1].M1 += dtodx1*(x1FD_i.Mx);
    pG->U[ks][ix.j][ix.i-1].M2 += dtodx1*(x1FD_i.My);
    pG->U[ks][ix.j][ix.i-1].M3 += dtodx1*(x1FD_i.Mz);
    pG->U[ks][ix.j][ix.i-1].E  += dtodx1*(x1FD_i.E );
  }

  if (ix.i < pG->ie) {
    pG->U[ks][ix.j][ix.i+1].d  -= dtodx1*(x1FD_ip1.d );
    pG->U[ks][ix.j][ix.i+1].M1 -= dtodx1*(x1FD_ip1.Mx);
    pG->U[ks][ix.j][ix.i+1].M2 -= dtodx1*(x1FD_ip1.My);
    pG->U[ks][ix.j][ix.i+1].M3 -= dtodx1*(x1FD_ip1.Mz);
    pG->U[ks][ix.j][ix.i+1].E  -= dtodx1*(x1FD_ip1.E );
  }

/* Use flux differences to correct cell-centered values at j-1 and j+1 */

  if (ix.j > pG->js) {
    pG->U[ks][ix.j-1][ix.i].d  += dtodx2*(x2FD_j.d );
    pG->U[ks][ix.j-1][ix.i].M1 += dtodx2*(x2FD_j.Mz);
    pG->U[ks][ix.j-1][ix.i].M2 += dtodx2*(x2FD_j.Mx);
    pG->U[ks][ix.j-1][ix.i].M3 += dtodx2*(x2FD_j.My);
    pG->U[ks][ix.j-1][ix.i].E  += dtodx2*(x2FD_j.E );
  }

  if (ix.j < pG->je) {
    pG->U[ks][ix.j+1][ix.i].d  -= dtodx2*(x2FD_jp1.d );
    pG->U[ks][ix.j+1][ix.i].M1 -= dtodx2*(x2FD_jp1.Mz);
    pG->U[ks][ix.j+1][ix.i].M2 -= dtodx2*(x2FD_jp1.Mx);
    pG->U[ks][ix.j+1][ix.i].M3 -= dtodx2*(x2FD_jp1.My);
    pG->U[ks][ix.j+1][ix.i].E  -= dtodx2*(x2FD_jp1.E );
  }
}
