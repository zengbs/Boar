/*============================================================================*/
/*! \file integrate.c
 *  \brief Contains public functions to set integrator.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_init()        - set pointer to integrate function based on dim
 * - integrate_destruct()    - call destruct integrate function based on dim */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../struct.h"
#include "prototypes.h"
#include "../prototypes.h"

/* dimension of calculation (determined at runtime) */
static int dim=0;

/*----------------------------------------------------------------------------*/
/*! \fn VDFun_t integrate_init(MeshS *pD)
 *  \brief Initialize integrator; VGDFun_t defined in athena.h   */
VDFun_t integrate_init(DomainS *pD)
{
  int i;
/* Calculate the dimensions (using root Domain)  */
  dim = 0;
  for (i=0; i<3; i++) if(pD->Grid->Nx[i] > 1) dim++;

/* set function pointer to appropriate integrator based on dimensions */
  switch(dim){

  case 1:
    if(pD->Grid->Nx[0] <= 1) break;
    integrate_init_1d(pD);
    return integrate_1d_vl;

  case 2:
    if(pD->Grid->Nx[2] > 1) break;
    integrate_init_2d(pD);
    return integrate_2d_vl;

  case 3:
    integrate_init_3d(pD);
    return integrate_3d_vl;
  }

  if (dim == 1)
     printf("[integrate_init]: 1D problem must have Nx1 > 1: Nx1=%d, Nx2=%d, Nx3=%d\n",
            pD->Grid->Nx[0],pD->Grid->Nx[1],pD->Grid->Nx[2]);
  if (dim == 2)
     printf("[integrate_init]: 2D problem must have Nx1 and Nx2 > 1: Nx1=%d, Nx2=%d, Nx3=%d\n",
            pD->Grid->Nx[0],pD->Grid->Nx[1],pD->Grid->Nx[2]);

/* This is never executed, but generates a warning on some compilers. */
  return NULL;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct()
 *  \brief Free memory */
void integrate_destruct()
{
  switch(dim){
  case 1:
    integrate_destruct_1d();
    return;
  case 2:
    integrate_destruct_2d();
    return;
  case 3:
    integrate_destruct_3d();
    return;
  }

  printf("[integrate_destruct]: Grid dimension = %d\n",dim);
}
