/*============================================================================*/
/*! \file integrate_3d_vl_sr.c
 *  \brief Integrate MHD equations using 3D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator. 
 *
 * PURPOSE: Integrate MHD equations using 3D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator.  The variables updated are:
 *    - U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *    - B1i, B2i, B3i  -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, and the H-correction
 *   of Sanders et al.
 *   - For adb hydro, requires (9*Cons1DS + 3*double + 1*ConsS) = 53 3D arrays
 *   - For adb mhd, requires   (9*Cons1DS + 9*double + 1*ConsS) = 80 3D arrays
 *
 * REFERENCE: 
 * - J.M Stone & T.A. Gardiner, "A simple, unsplit Godunov method
 *   for multidimensional MHD", NewA 14, 139 (2009)
 *
 * - R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_3d_vl()
 * - integrate_destruct_3d()
 * - integrate_init_3d() */
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
static Prim1DS ***Wl_x1Face=NULL, ***Wr_x1Face=NULL;
static Prim1DS ***Wl_x2Face=NULL, ***Wr_x2Face=NULL;
static Prim1DS ***Wl_x3Face=NULL, ***Wr_x3Face=NULL;
static Cons1DS ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;


/* 1D scratch vectors used by lr_states and flux functions */
static double *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W1d=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* primitive variables at t^{n} computed in predict step */
static PrimS ***W=NULL;

/* conserved & primitive variables at t^{n+1/2} computed in predict step */
static ConsS ***Uhalf=NULL;
static PrimS ***Whalf=NULL;


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_3d_vl(DomainS *pD) 
 *  \brief 3D van Leer unsplit integrator for MHD. 
 */
void integrate_3d_vl(DomainS *pD)
{
  GridS *pG=pD->Grid;
  pG->dt=pD->dt;
  pG->time=pD->time;
  double dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  double q1 = 0.5*dtodx1, q2 = 0.5*dtodx2, q3 = 0.5*dtodx3;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int cart_x1 = 1, cart_x2 = 2, cart_x3 = 3;

  int il=is-(nghost-1), iu=ie+(nghost-1);
  int jl=js-(nghost-1), ju=je+(nghost-1);
  int kl=ks-(nghost-1), ku=ke+(nghost-1);


  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        Uhalf[k][j][i] = pG->U[k][j][i];
        W[k][j][i] = Cons_to_Prim(&(pG->U[k][j][i]));
      }
    }
  }

/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, U1, U2, U3, P, B2c, B3c, s[n])
 */

  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
	W1d[i].d  = W[k][j][i].d;
	W1d[i].Ux = W[k][j][i].U1;
	W1d[i].Uy = W[k][j][i].U2;
	W1d[i].Uz = W[k][j][i].U3;
	W1d[i].P  = W[k][j][i].P;
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

/*--- Step 1d ------------------------------------------------------------------
 * Compute flux in x1-direction */

      for (i=il; i<=ie+nghost; i++) {
        fluxes(Ul[i],Ur[i],Wl[i],Wr[i],&x1Flux[k][j][i]);
      }
    }
  }

/*=== STEP 2: Compute first-order fluxes at t^{n} in x2-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, U2, U3, U1, P, B3c, B1c, s[n])
 */

  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      for (j=js-nghost; j<=je+nghost; j++) {
	W1d[j].d  = W[k][j][i].d;
	W1d[j].Ux = W[k][j][i].U2;
	W1d[j].Uy = W[k][j][i].U3;
	W1d[j].Uz = W[k][j][i].U1;
	W1d[j].P  = W[k][j][i].P;
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

/*--- Step 2d ------------------------------------------------------------------
 * Compute flux in x2-direction */

      for (j=jl; j<=je+nghost; j++) {
        fluxes(Ul[j],Ur[j],Wl[j],Wr[j],&x2Flux[k][j][i]);
      }
    }
  }


/*=== STEP 3: Compute first-order fluxes at t^{n} in x3-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 3a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, U3, U1, U2, P, B1c, B2c, s[n])
 */

  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      for (k=ks-nghost; k<=ke+nghost; k++) {
	W1d[k].d  = W[k][j][i].d;
	W1d[k].Ux = W[k][j][i].U3;
	W1d[k].Uy = W[k][j][i].U1;
	W1d[k].Uz = W[k][j][i].U2;
	W1d[k].P  = W[k][j][i].P;
      }

/*--- Step 3b ------------------------------------------------------------------
 * Compute first-order L/R states */      
      for (k=kl; k<=ke+nghost; k++) { 
        Wl[k] = W1d[k-1];
        Wr[k] = W1d[k  ]; 

/* Compute U from W in case Pfloor used in Cons1D_to_Prim1D */
        Ul[k] = Prim1D_to_Cons1D(&Wl[k]);
        Ur[k] = Prim1D_to_Cons1D(&Wr[k]);
      }

/*--- Step 3d ------------------------------------------------------------------
 * Compute flux in x1-direction */

      for (k=kl; k<=ke+nghost; k++) {
        fluxes(Ul[k],Ur[k],Wl[k],Wr[k],&x3Flux[k][j][i]);
      }
    }
  }


/*=== STEP 5: Update cell-centered variables to half-timestep ================*/

/*--- Step 5a ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x1-fluxes
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Uhalf[k][j][i].d   -= q1*(x1Flux[k][j][i+1].d  - x1Flux[k][j][i].d );
        Uhalf[k][j][i].M1  -= q1*(x1Flux[k][j][i+1].Mx - x1Flux[k][j][i].Mx);
        Uhalf[k][j][i].M2  -= q1*(x1Flux[k][j][i+1].My - x1Flux[k][j][i].My);
        Uhalf[k][j][i].M3  -= q1*(x1Flux[k][j][i+1].Mz - x1Flux[k][j][i].Mz);
        Uhalf[k][j][i].E   -= q1*(x1Flux[k][j][i+1].E  - x1Flux[k][j][i].E );
      }
    }
  }

/*--- Step 5b ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x2-fluxes
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Uhalf[k][j][i].d   -= q2*(x2Flux[k][j+1][i].d  - x2Flux[k][j][i].d );
        Uhalf[k][j][i].M1  -= q2*(x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz);
        Uhalf[k][j][i].M2  -= q2*(x2Flux[k][j+1][i].Mx - x2Flux[k][j][i].Mx);
        Uhalf[k][j][i].M3  -= q2*(x2Flux[k][j+1][i].My - x2Flux[k][j][i].My);
        Uhalf[k][j][i].E   -= q2*(x2Flux[k][j+1][i].E  - x2Flux[k][j][i].E );
      }
    }
  }

/*--- Step 5c ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x3-fluxes
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Uhalf[k][j][i].d   -= q3*(x3Flux[k+1][j][i].d  - x3Flux[k][j][i].d );
        Uhalf[k][j][i].M1  -= q3*(x3Flux[k+1][j][i].My - x3Flux[k][j][i].My);
        Uhalf[k][j][i].M2  -= q3*(x3Flux[k+1][j][i].Mz - x3Flux[k][j][i].Mz);
        Uhalf[k][j][i].M3  -= q3*(x3Flux[k+1][j][i].Mx - x3Flux[k][j][i].Mx);
        Uhalf[k][j][i].E   -= q3*(x3Flux[k+1][j][i].E  - x3Flux[k][j][i].E );
      }
    }
  }

        
  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      for (j=js-nghost; j<=je+nghost; j++) {
	//Whalf[k][j][i] = check_Prim(&(Uhalf[k][j][i]));
	Whalf[k][j][i] = Cons_to_Prim(&(Uhalf[k][j][i]));
      }
    }
  }


/*=== STEP 8: Compute second-order L/R x1-interface states ===================*/

/*--- Step 8a ------------------------------------------------------------------
 * Load 1D vector of primitve variables;
 * W1d = (d, U1, U2, U3, P, B2c, B3c, s[n])
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=il; i<=iu; i++) {
        W1d[i].d  = Whalf[k][j][i].d;
        W1d[i].Ux = Whalf[k][j][i].U1;
        W1d[i].Uy = Whalf[k][j][i].U2;
        W1d[i].Uz = Whalf[k][j][i].U3;
        W1d[i].P  = Whalf[k][j][i].P;
      }

/*--- Step 8b ------------------------------------------------------------------
 * Compute L and R states at x1-interfaces, store in 3D array
 */

      lr_states(pG,W1d,Bxc,pG->dt,pG->dx1,is,ie,Wl,Wr,cart_x1);


      for (i=is; i<=ie+1; i++) {
        Wl_x1Face[k][j][i] = Wl[i];
        Wr_x1Face[k][j][i] = Wr[i];
      }
    }
  }

/*=== STEP 9: Compute second-order L/R x2-interface states ===================*/

/*--- Step 9a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * W1d = (d, U2, U3, U1, P, B3c, B1c, s[n])
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (i=is-1; i<=ie+1; i++) {
      for (j=jl; j<=ju; j++) {
        W1d[j].d  = Whalf[k][j][i].d;
        W1d[j].Ux = Whalf[k][j][i].U2;
        W1d[j].Uy = Whalf[k][j][i].U3;
        W1d[j].Uz = Whalf[k][j][i].U1;
        W1d[j].P  = Whalf[k][j][i].P;
      }

/*--- Step 9b ------------------------------------------------------------------
 * Compute L and R states at x2-interfaces, store in 3D array
 */

      lr_states(pG,W1d,Bxc,pG->dt,pG->dx2,js,je,Wl,Wr,cart_x2);

      for (j=js; j<=je+1; j++) {
        Wl_x2Face[k][j][i] = Wl[j];
        Wr_x2Face[k][j][i] = Wr[j];
      }
    }
  }

/*=== STEP 10: Compute second-order L/R x3-interface states ==================*/

/*--- Step 9a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * W1d = (d, U3, U1, U2, P, B1c, B2c, s[n])
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      for (k=kl; k<=ku; k++) {
        W1d[k].d  = Whalf[k][j][i].d;
        W1d[k].Ux = Whalf[k][j][i].U3;
        W1d[k].Uy = Whalf[k][j][i].U1;
        W1d[k].Uz = Whalf[k][j][i].U2;
        W1d[k].P  = Whalf[k][j][i].P;
      }

/*--- Step 9b ------------------------------------------------------------------
 * Compute L and R states at x3-interfaces, store in 3D array
 */
      lr_states(pG,W1d,Bxc,pG->dt,pG->dx3,ks,ke,Wl,Wr,cart_x3);


      for (k=ks; k<=ke+1; k++) {
        Wl_x3Face[k][j][i] = Wl[k];
        Wr_x3Face[k][j][i] = Wr[k];
      }
    }
  }


/*--- Step 11b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        Ul[i] = Prim1D_to_Cons1D(&Wl_x1Face[k][j][i]);
        Ur[i] = Prim1D_to_Cons1D(&Wr_x1Face[k][j][i]);
        fluxes(Ul[i],Ur[i],Wl_x1Face[k][j][i],Wr_x1Face[k][j][i],
               &x1Flux[k][j][i]);
      }
    }
  }

/*--- Step 11c -----------------------------------------------------------------
 * Compute second-order fluxes in x2-direction
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        Ul[i] = Prim1D_to_Cons1D(&Wl_x2Face[k][j][i]);
        Ur[i] = Prim1D_to_Cons1D(&Wr_x2Face[k][j][i]);

        fluxes(Ul[i],Ur[i],Wl_x2Face[k][j][i],Wr_x2Face[k][j][i],
               &x2Flux[k][j][i]);
      }
    }
  }

/*--- Step 11d -----------------------------------------------------------------
 * Compute second-order fluxes in x3-direction
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        Ul[i] = Prim1D_to_Cons1D(&Wl_x3Face[k][j][i]);
        Ur[i] = Prim1D_to_Cons1D(&Wr_x3Face[k][j][i]);

        fluxes(Ul[i],Ur[i],Wl_x3Face[k][j][i],Wr_x3Face[k][j][i],
               &x3Flux[k][j][i]);
      }
    }
  }

/*=== STEP 12: Update face-centered B for a full timestep ====================*/
        
/*--- Step 12a -----------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at the half-time-step.
 */



  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].d -=dtodx1*(x1Flux[k][j][i+1].d -x1Flux[k][j][i].d );
        pG->U[k][j][i].M1-=dtodx1*(x1Flux[k][j][i+1].Mx-x1Flux[k][j][i].Mx);
        pG->U[k][j][i].M2-=dtodx1*(x1Flux[k][j][i+1].My-x1Flux[k][j][i].My);
        pG->U[k][j][i].M3-=dtodx1*(x1Flux[k][j][i+1].Mz-x1Flux[k][j][i].Mz);
        pG->U[k][j][i].E -=dtodx1*(x1Flux[k][j][i+1].E -x1Flux[k][j][i].E );
      }
    }
  }

/*--- Step 14b -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x2-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].d -=dtodx2*(x2Flux[k][j+1][i].d -x2Flux[k][j][i].d );
        pG->U[k][j][i].M1-=dtodx2*(x2Flux[k][j+1][i].Mz-x2Flux[k][j][i].Mz);
        pG->U[k][j][i].M2-=dtodx2*(x2Flux[k][j+1][i].Mx-x2Flux[k][j][i].Mx);
        pG->U[k][j][i].M3-=dtodx2*(x2Flux[k][j+1][i].My-x2Flux[k][j][i].My);
        pG->U[k][j][i].E -=dtodx2*(x2Flux[k][j+1][i].E -x2Flux[k][j][i].E );
      }
    }
  }

/*--- Step 14c -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x3-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].d -=dtodx3*(x3Flux[k+1][j][i].d -x3Flux[k][j][i].d );
        pG->U[k][j][i].M1-=dtodx3*(x3Flux[k+1][j][i].My-x3Flux[k][j][i].My);
        pG->U[k][j][i].M2-=dtodx3*(x3Flux[k+1][j][i].Mz-x3Flux[k][j][i].Mz);
        pG->U[k][j][i].M3-=dtodx3*(x3Flux[k+1][j][i].Mx-x3Flux[k][j][i].Mx);
        pG->U[k][j][i].E -=dtodx3*(x3Flux[k+1][j][i].E -x3Flux[k][j][i].E );
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_3d(MeshS *pM)
 *  \brief Allocate temporary integration arrays */
void integrate_init_3d(DomainS *pD)
{
  int nmax,size1=0,size2=0,size3=0;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
//  for (nl=0; nl<(pD->NLevels); nl++){
//    for (nd=0; nd<(pD->DomainsPerLevel[nl]); nd++){
      if (pD->Grid != NULL) {
        if (pD->Grid->Nx[0] > size1){
          size1 = pD->Grid->Nx[0];
        }
        if (pD->Grid->Nx[1] > size2){
          size2 = pD->Grid->Nx[1];
        }
        if (pD->Grid->Nx[2] > size3){
          size3 = pD->Grid->Nx[2];
        }
      }
//    }
//  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

  if ((Wl_x1Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x1Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wl_x2Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x2Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wl_x3Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x3Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;

  if ((Bxc = (double*)malloc(nmax*sizeof(double))) == NULL) goto on_error;
  if ((Bxi = (double*)malloc(nmax*sizeof(double))) == NULL) goto on_error;

  if ((U1d = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul  = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur  = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W1d = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((x1Flux = (Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2Flux = (Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3Flux = (Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;

  if ((Uhalf = (ConsS***)calloc_3d_array(size3,size2,size1,sizeof(ConsS)))
    == NULL) goto on_error;
  if ((Whalf = (PrimS***)calloc_3d_array(size3,size2,size1,sizeof(PrimS)))
    == NULL) goto on_error;
  if ((W     = (PrimS***)calloc_3d_array(size3,size2,size1,sizeof(PrimS)))
    == NULL) goto on_error;

  return;

  on_error:
    integrate_destruct();
    printf("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_3d(void)
 *  \brief Free temporary integration arrays */
void integrate_destruct_3d(void)
{
  if (Wl_x1Face != NULL) free_3d_array(Wl_x1Face);
  if (Wr_x1Face != NULL) free_3d_array(Wr_x1Face);
  if (Wl_x2Face != NULL) free_3d_array(Wl_x2Face);
  if (Wr_x2Face != NULL) free_3d_array(Wr_x2Face);
  if (Wl_x3Face != NULL) free_3d_array(Wl_x3Face);
  if (Wr_x3Face != NULL) free_3d_array(Wr_x3Face);

  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);

  if (U1d != NULL) free(U1d);
  if (Ul  != NULL) free(Ul);
  if (Ur  != NULL) free(Ur);
  if (W1d != NULL) free(W1d);
  if (Wl  != NULL) free(Wl);
  if (Wr  != NULL) free(Wr);

  if (x1Flux  != NULL) free_3d_array(x1Flux);
  if (x2Flux  != NULL) free_3d_array(x2Flux);
  if (x3Flux  != NULL) free_3d_array(x3Flux);

  if (Uhalf  != NULL) free_3d_array(Uhalf);
  if (Whalf  != NULL) free_3d_array(Whalf);
  if (W      != NULL) free_3d_array(W);

  return;
}
