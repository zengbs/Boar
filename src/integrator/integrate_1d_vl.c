/*============================================================================*/
/*! \file integrate_1d_vl_sr.c
 *  \brief Integrate SRMHD equations using 1D version of MUSCL-Hancock (VL)
 *   integrator.
 *
 * PURPOSE: Integrate SRMHD equations using 1D version of MUSCL-Hancock (VL)
 *   integrator.  Updates U.[d,M1,M2,M3,E,B2c,B3c,s] in Grid structure.
 *   Adds gravitational source terms, self-gravity.
 *        
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_1d_vl()
 * - integrate_init_1d()
 * - integrate_destruct_1d() */
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
static Prim1DS *Wl_x1Face=NULL, *Wr_x1Face=NULL;
static Cons1DS *x1Flux=NULL;

/* 1D scratch vectors used by lr_states and flux functions */
static double *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W1d=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;
static Cons1DS *x1FluxP=NULL;

/* conserved & primitive variables at t^{n} */
static PrimS *W=NULL;

/* conserved variables at t^{n+1/2} computed in predict step */
static ConsS *Uhalf=NULL;
static PrimS *Whalf=NULL;

static void FixCell(GridS *pG, Int3Vect);


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_1d_vl(DomainS *pD) 
 *  \brief 1D version of van Leer unsplit integrator for MHD. 
 *
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 1D.
 */
void integrate_1d_vl(DomainS *pD)
{
  GridS *pG=pD->Grid;
  pG->dt=pD->dt;
  pG->time=pD->time;
//  ConsS U;
  double dtodx1=pG->dt/pG->dx1, hdtodx1=0.5*pG->dt/pG->dx1;
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  int cart_x1 = 1;//, cart_x2 = 2, cart_x3 = 3;
  //double x1,x2,x3,phicl,phicr,phifc,phil,phir,phic;
  int flag_cell=0,negd=0,negP=0;//,NaNFlux=0;
//  int fail=0,final=0;
//  double Bx;
  PrimS Wcheck;
  Int3Vect BadCell;

  int il=is-(nghost-1), iu=ie+(nghost-1);

  for (i=is-nghost; i<=ie+nghost; i++) {
    Uhalf[i] = pG->U[ks][js][i];
    W[i] = Cons_to_Prim(&(pG->U[ks][js][i]));
  }

/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;  
 * W1d = (d, U1, U2, U3, P, B2c, B3c, s[n])
 */

  for (i=is-nghost; i<=ie+nghost; i++) {
    W1d[i].d  = W[i].d;
    W1d[i].Ux = W[i].U1;
    W1d[i].Uy = W[i].U2;
    W1d[i].Uz = W[i].U3;
    W1d[i].P  = W[i].P;
  }

/*--- Step 1b ------------------------------------------------------------------
 * Compute first-order L/R states */

/* Ensure that W & U are consistent */
  for (i=is-nghost; i<=ie+nghost; i++) {
    U1d[i] = Prim1D_to_Cons1D(&W1d[i]);
  }

  for (i=il; i<=ie+nghost; i++) {
    Wl[i] = W1d[i-1];
    Wr[i] = W1d[i  ];

    Ul[i] = U1d[i-1];
    Ur[i] = U1d[i  ];
  }

/*--- Step 1c ------------------------------------------------------------------
 * No source terms needed */

/*--- Step 1d ------------------------------------------------------------------
 * Compute flux in x1-direction */

  for (i=il; i<=ie+nghost; i++) {
    fluxes(Ul[i],Ur[i],Wl[i],Wr[i],&x1Flux[i]);
  }

/*=== STEPS 2-4: Not needed in 1D ===*/

/*=== STEP 5: Update cell-centered variables to half-timestep ================*/

/*--- Step 5a ------------------------------------------------------------------
 * Update cell-centered variables (including B2c and B3c) to half-timestep
 */

  for (i=il; i<=iu; i++) {
    Uhalf[i].d   -= hdtodx1*(x1Flux[i+1].d  - x1Flux[i].d );
    Uhalf[i].M1  -= hdtodx1*(x1Flux[i+1].Mx - x1Flux[i].Mx);
    Uhalf[i].M2  -= hdtodx1*(x1Flux[i+1].My - x1Flux[i].My);
    Uhalf[i].M3  -= hdtodx1*(x1Flux[i+1].Mz - x1Flux[i].Mz);
    Uhalf[i].E   -= hdtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
          x1FluxP[i] = x1Flux[i];
  }

/*=== STEP 6: Add source terms to predict values at half-timestep ============*/

/*--- Step 6a ------------------------------------------------------------------
 * Add source terms from a static gravitational potential for 0.5*dt to predict
 * step.  To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

/*=== STEP 7: Conserved->Primitive variable inversion at t^{n+1/2} ===========*/
        
/* Invert conserved variables at t^{n+1/2} to primitive variables. With FOFC, 
 * if cell-centered d < 0, P< 0, or v^2 > 1, correct by switching back to 
 * values at beginning of step, rendering update first order in time for that
 * cell.
 */

  negd = 0;
  negP = 0;
  flag_cell = 0;
  for (i=il; i<=iu; i++) {
    Whalf[i] = Cons_to_Prim(&Uhalf[i]);
    if (Whalf[i].d < 0.0) {
      flag_cell = 1;
      negd++;
    }
    if (Whalf[i].P < 0.0) {
      flag_cell = 1;
      negP++;
    }
    if (flag_cell != 0) {
      Whalf[i].d = W[i].d;
      Whalf[i].U1 = W[i].U1;
      Whalf[i].U2 = W[i].U2;
      Whalf[i].U3 = W[i].U3;
      Whalf[i].P = W[i].P;
      flag_cell=0;
    }
  }

  if (negd > 0 || negP > 0)
    printf("[Step7]: %i cells had d<0; %i cells had P<0\n"
                                 ,negd,negP);


/*=== STEP 8: Compute second-order L/R x1-interface states ===================*/

/*--- Step 8a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W = (d, U1, U2, U3, P, B2c, B3c, s[n])
 */

  for (i=il; i<=iu; i++) {
    W1d[i].d  = Whalf[i].d;
    W1d[i].Ux = Whalf[i].U1;
    W1d[i].Uy = Whalf[i].U2;
    W1d[i].Uz = Whalf[i].U3;
    W1d[i].P  = Whalf[i].P;
  }

/*--- Step 8b ------------------------------------------------------------------
 * Compute L/R states on x1-interfaces, store into arrays
 */

  lr_states(pG,W1d,Bxc,pG->dt,pG->dx1,is,ie,Wl,Wr,cart_x1);

  for (i=is; i<=ie+1; i++) {
    Wl_x1Face[i] = Wl[i];
    Wr_x1Face[i] = Wr[i];
  }

/*=== STEPS 9-10: Not needed in 1D ===*/

/*=== STEP 11: Compute x1-Flux ===============================================*/

/*--- Step 11b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

  for (i=is; i<=ie+1; i++) {
    Ul[i] = Prim1D_to_Cons1D(&Wl_x1Face[i]);
    Ur[i] = Prim1D_to_Cons1D(&Wr_x1Face[i]);

    fluxes(Ul[i],Ur[i],Wl_x1Face[i],Wr_x1Face[i],&x1Flux[i]);
  }

/*=== STEP 12: Not needed in 1D ===*/
        
/*=== STEP 13: Add source terms for a full timestep using n+1/2 states =======*/
       
/*--- Step 13a -----------------------------------------------------------------
 * Add gravitational source terms due to a Static Potential
 * To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

/*=== STEP 14: Update cell-centered values for a full timestep ===============*/

/*--- Step 14a -----------------------------------------------------------------
 * Update cell-centered variables in pG (including B2c and B3c) using x1-Fluxes
 */

  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].d  -= dtodx1*(x1Flux[i+1].d  - x1Flux[i].d );
    pG->U[ks][js][i].M1 -= dtodx1*(x1Flux[i+1].Mx - x1Flux[i].Mx);
    pG->U[ks][js][i].M2 -= dtodx1*(x1Flux[i+1].My - x1Flux[i].My);
    pG->U[ks][js][i].M3 -= dtodx1*(x1Flux[i+1].Mz - x1Flux[i].Mz);
    pG->U[ks][js][i].E  -= dtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
  }

/*=== STEP 15: First-order flux correction ===================================*/

/*--- Step 15a -----------------------------------------------------------------
 * If cell-centered d or P have gone negative, or if v^2 > 1, correct
 * by using 1st order predictor fluxes */
        
  for (i=is; i<=ie; i++) {
      Wcheck = Cons_to_Prim(&(pG->U[ks][js][i]));
      if (Wcheck.d < 0.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = js;
        BadCell.k = ks;
        negd++;
      }
      if (Wcheck.P < 0.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = js;
        BadCell.k = ks;
        negP++;
      }
      if (flag_cell != 0) {
        FixCell(pG, BadCell);
        flag_cell=0;
      }

  }

  if (negd > 0 || negP > 0){
    printf("[Step15a]: %i cells had d<0; %i cells had P<0;\n",negd,negP);
  }

/*--- Step 15b -----------------------------------------------------------------
 * In SR the first-order flux correction can fail to fix an unphysical state.
 * We must fix these cells in order to avoid NaN's at the next timestep,
 * particuarly if v^2 > 1. We have 2 approaches; firstly, we use the entropy
 * equation (which we have not applied a 1st order flux correction to) to
 * calculate the pressure and the Lorentz factor of the gas. If this produces
 * and unphysical state, then we floor the pressure and iterate on v^2 until
 * v^2 < 1. Possibly could improved by averaging density and pressure from
 * adjacent cells and then calculating pressure.
 */
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_1d(MeshS *pM)
 *  \brief Allocate temporary integration arrays */
void integrate_init_1d(DomainS *pD)
{
  int size1=0;//,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1 */
      if (pD->Grid != NULL) {
        if (pD->Grid->Nx[0] > size1){
          size1 = pD->Grid->Nx[0];
        }
      }

  size1 = size1 + 2*nghost;

  if ((Wl_x1Face=(Prim1DS*)malloc(size1*sizeof(Prim1DS))) ==NULL) goto on_error;
  if ((Wr_x1Face=(Prim1DS*)malloc(size1*sizeof(Prim1DS))) ==NULL) goto on_error;
  if ((x1Flux   =(Cons1DS*)malloc(size1*sizeof(Cons1DS))) ==NULL) goto on_error;
  if ((x1FluxP  =(Cons1DS*)malloc(size1*sizeof(Cons1DS))) ==NULL) goto on_error;

  if ((Bxc = (double*)malloc(size1*sizeof(double))) == NULL) goto on_error;
  if ((Bxi = (double*)malloc(size1*sizeof(double))) == NULL) goto on_error;

  if ((U1d= (Cons1DS*)malloc(size1*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul = (Cons1DS*)malloc(size1*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur = (Cons1DS*)malloc(size1*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W1d= (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Uhalf = (ConsS*)malloc(size1*sizeof(ConsS)))==NULL) goto on_error;
  if ((Whalf = (PrimS*)malloc(size1*sizeof(PrimS)))==NULL) goto on_error;
  if ((W     = (PrimS*)malloc(size1*sizeof(PrimS)))==NULL) goto on_error;

  return;

  on_error:
    integrate_destruct();
    printf("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_1d(void)
 *  \brief Free temporary integration arrays */
void integrate_destruct_1d(void)
{
  if (Wl_x1Face != NULL) free(Wl_x1Face);
  if (Wr_x1Face != NULL) free(Wr_x1Face);
  if (x1Flux    != NULL) free(x1Flux);
  if (x1FluxP   != NULL) free(x1FluxP);


  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W1d      != NULL) free(W1d);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Uhalf    != NULL) free(Uhalf);
  if (Whalf    != NULL) free(Whalf);
  if (W        != NULL) free(W);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void FixCell(GridS *pG, Int3Vect indx)
 *  \brief Uses first order fluxes to fix negative d,P or superluminal v
 */ 
static void FixCell(GridS *pG, Int3Vect indx)
{
  int ks=pG->ks,js=pG->js;
  double dtodx1=pG->dt/pG->dx1;
  Cons1DS x1FD_i, x1FD_ip1;//, x2FD_j, x2FD_jp1;
  
  /* Compute difference of predictor and corrector fluxes at cell faces */
  
  x1FD_i.d    = x1Flux[indx.i  ].d - x1FluxP[indx.i  ].d;
  x1FD_ip1.d  = x1Flux[indx.i+1].d - x1FluxP[indx.i+1].d;
        
  x1FD_i.Mx    = x1Flux[indx.i  ].Mx - x1FluxP[indx.i  ].Mx;
  x1FD_ip1.Mx  = x1Flux[indx.i+1].Mx - x1FluxP[indx.i+1].Mx;
        
  x1FD_i.My    = x1Flux[indx.i  ].My - x1FluxP[indx.i  ].My;
  x1FD_ip1.My  = x1Flux[indx.i+1].My - x1FluxP[indx.i+1].My;
        
  x1FD_i.Mz    = x1Flux[indx.i  ].Mz - x1FluxP[indx.i  ].Mz;
  x1FD_ip1.Mz  = x1Flux[indx.i+1].Mz - x1FluxP[indx.i+1].Mz;
        
  x1FD_i.E    = x1Flux[indx.i  ].E - x1FluxP[indx.i  ].E;
  x1FD_ip1.E  = x1Flux[indx.i+1].E - x1FluxP[indx.i+1].E;
        
  /* Use flux differences to correct bad cell */
  pG->U[ks][js][indx.i].d  += dtodx1*(x1FD_ip1.d  - x1FD_i.d );
  pG->U[ks][js][indx.i].M1 += dtodx1*(x1FD_ip1.Mx - x1FD_i.Mx);
  pG->U[ks][js][indx.i].M2 += dtodx1*(x1FD_ip1.My - x1FD_i.My);
  pG->U[ks][js][indx.i].M3 += dtodx1*(x1FD_ip1.Mz - x1FD_i.Mz);
  pG->U[ks][js][indx.i].E  += dtodx1*(x1FD_ip1.E  - x1FD_i.E );
        
        
  /* Use flux differences to correct bad cell neighbors at i-1 and i+1 */      
  if (indx.i > pG->is) {
    pG->U[ks][js][indx.i-1].d  += dtodx1*(x1FD_i.d );
    pG->U[ks][js][indx.i-1].M1 += dtodx1*(x1FD_i.Mx);
    pG->U[ks][js][indx.i-1].M2 += dtodx1*(x1FD_i.My);
    pG->U[ks][js][indx.i-1].M3 += dtodx1*(x1FD_i.Mz);
    pG->U[ks][js][indx.i-1].E  += dtodx1*(x1FD_i.E );
  }
        
  if (indx.i < pG->ie) {
    pG->U[ks][js][indx.i+1].d  -= dtodx1*(x1FD_ip1.d );
    pG->U[ks][js][indx.i+1].M1 -= dtodx1*(x1FD_ip1.Mx);
    pG->U[ks][js][indx.i+1].M2 -= dtodx1*(x1FD_ip1.My);
    pG->U[ks][js][indx.i+1].M3 -= dtodx1*(x1FD_ip1.Mz);
    pG->U[ks][js][indx.i+1].E  -= dtodx1*(x1FD_ip1.E );
  }
        
}
