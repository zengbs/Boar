/*============================================================================*/
/*! \file init_grid.c 
 *  \brief Initializes most variables in the Grid structure.
 *
 * PURPOSE: Initializes most variables in the Grid structure.  Allocates memory
 *   for 3D arrays of Cons, interface B, etc.  With SMR, finds all overlaps
 *   between child and parent Grids, and initializes data needed for restriction
 *   flux-correction, and prolongation steps.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - init_grid()
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - checkOverlap() - checks for overlap of cubes, and returns overlap coords
 * - checkOverlapTouch() - same as above, but checks for overlap and/or touch */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "struct.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *  checkOverlap() - checks for overlap of cubes, and returns overlap coords
 *  checkOverlapTouch() - same as above, but checks for overlap and/or touch
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn void init_grid(MeshS *pM)
 *  \brief Initializes most variables in the Grid structure.
 */

void init_grid(DomainS *pD, int nz[])
{
  GridS *pG = NULL;
  pG = pD->Grid;

  int n1z, n2z, n3z;
/* number of dimensions in Grid. */

/* ---------------------  Intialize grid in 1-direction --------------------- */
/* Initialize is,ie,dx1
 * Compute Disp, MinX[0], and MaxX[0] using displacement of Domain and Grid
 * location within Domain */
      pG->Nx[0] = nz[0];

      if(pG->Nx[0] > 1) {
        pG->is = nghost;
        pG->ie = pG->Nx[0] + nghost - 1;
      }
      else
        pG->is = pG->ie = 0;
    
      pG->MaxX[0] = +1.0;
      pG->MinX[0] = +0.0;

      pG->dx1 = (pG->MaxX[0]-pG->MinX[0])/pG->Nx[0];
/* ---------------------  Intialize grid in 2-direction --------------------- */
/* Initialize js,je,dx2
 * Compute Disp, MinX[1], and MaxX[1] using displacement of Domain and Grid
 * location within Domain */

      pG->Nx[1] = nz[1];

      if(pG->Nx[1] > 1) {
        pG->js = nghost;
        pG->je = pG->Nx[1] + nghost - 1;
      }
      else
        pG->js = pG->je = 0;
    
      pG->MaxX[1] = +1.0;
      pG->MinX[1] = +0.0;

      pG->dx2 = (pG->MaxX[1]-pG->MinX[1])/pG->Nx[1];

/* ---------------------  Intialize grid in 3-direction --------------------- */
/* Initialize ks,ke,dx3
 * Compute Disp, MinX[2], and MaxX[2] using displacement of Domain and Grid
 * location within Domain */

      pG->Nx[2] = nz[2];

      if(pG->Nx[2] > 1) {
        pG->ks = nghost;
        pG->ke = pG->Nx[2] + nghost - 1;
      }
      else
        pG->ks = pG->ke = 0;
    
      pG->MaxX[2] = +1.0;
      pG->MinX[2] = +0.0;

      pG->dx3 = (pG->MaxX[2]-pG->MinX[2])/pG->Nx[2];
/* ---------  Allocate 3D arrays to hold Cons based on size of grid --------- */

      if (pG->Nx[0] > 1)
        n1z = pG->Nx[0] + 2*nghost;
      else
        n1z = 1;

      if (pG->Nx[1] > 1)
        n2z = pG->Nx[1] + 2*nghost;
      else
        n2z = 1;

      if (pG->Nx[2] > 1)
        n3z = pG->Nx[2] + 2*nghost;
      else
        n3z = 1;

/* Build a 3D array of type ConsS */

      pG->U = (ConsS***)calloc_3d_array(n3z, n2z, n1z, sizeof(ConsS));
      if (pG->U == NULL) goto on_error1;

return;

/*--- Error messages ---------------------------------------------------------*/

on_error1:
    free_3d_array(pG->U);
    printf ("[init_grid]: Error allocating memory\n");
}
