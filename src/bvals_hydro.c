#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "struct.h"
#include "globals.h"
#include "prototypes.h"

static void reflect_ix1(GridS *pG);
static void reflect_ox1(GridS *pG);
static void reflect_ix2(GridS *pG);
static void reflect_ox2(GridS *pG);
static void reflect_ix3(GridS *pG);
static void reflect_ox3(GridS *pG);

static void outflow_ix1(GridS *pG);
static void outflow_ox1(GridS *pG);
static void outflow_ix2(GridS *pG);
static void outflow_ox2(GridS *pG);
static void outflow_ix3(GridS *pG);
static void outflow_ox3(GridS *pG);

static void periodic_ix1(GridS *pG);
static void periodic_ox1(GridS *pG);
static void periodic_ix2(GridS *pG);
static void periodic_ox2(GridS *pG);
static void periodic_ix3(GridS *pG);
static void periodic_ox3(GridS *pG);

static void conduct_ix1(GridS *pG);
static void conduct_ox1(GridS *pG);
static void conduct_ix2(GridS *pG);
static void conduct_ox2(GridS *pG);
static void conduct_ix3(GridS *pG);
static void conduct_ox3(GridS *pG);

void
bvals_hydro (DomainS * pD)
{
  GridS *pG = NULL;
  pG = pD->Grid;


if (pG->Nx[0] > 1){
/*--- Step 1.
 * Boundary Conditions in x1-direction */
  (*(pD->ix1_BCFun)) (pG);
  (*(pD->ox1_BCFun)) (pG);
}


if (pG->Nx[1] > 1){
/*--- Step 2.
 * Boundary Conditions in x2-direction */
  (*(pD->ix2_BCFun)) (pG);
  (*(pD->ox2_BCFun)) (pG);
}

if (pG->Nx[2] > 1){
/*--- Step 3.
 * Boundary Conditions in x3-direction */
  (*(pD->ix3_BCFun)) (pG);
  (*(pD->ox3_BCFun)) (pG);
  }
}


void
bvals_hydro_init (DomainS * pD)
{
  GridS *pG = NULL;
  pG = pD->Grid;

  int BCFlag_ix1, BCFlag_ox1;	/*!< BC flag on grid for inner/outer x1 */
  int BCFlag_ix2, BCFlag_ox2;	/*!< BC flag on grid for inner/outer x2 */
  int BCFlag_ix3, BCFlag_ox3;	/*!< BC flag on grid for inner/outer x3 */

  BCFlag_ix1 = 2;
  BCFlag_ox1 = 2;
  BCFlag_ix2 = 2;
  BCFlag_ox2 = 2;
  BCFlag_ix3 = 2;
  BCFlag_ox3 = 2;

/*---- ix1 boundary ----------------------------------------------------------*/

if (pG->Nx[0] > 1){

  switch (BCFlag_ix1)
    {
    case 1:
      pD->ix1_BCFun = reflect_ix1;
      break;

    case 2:			/* Outflow */
      pD->ix1_BCFun = outflow_ix1;
      break;

    case 4:			/* Periodic. Handle with MPI calls for parallel jobs. */
      pD->ix1_BCFun = periodic_ix1;
      break;

    case 5:			/* Reflecting, B_normal!=0 */
      pD->ix1_BCFun = conduct_ix1;
      break;
/*
    default:
      printf (-1, "[bvals_init]:bc_ix1=%d unknown\n", pD->BCFlag_ix1);
      exit (EXIT_FAILURE);
*/
    }

/*---- ox1 boundary ----------------------------------------------------------*/

  switch (BCFlag_ox1)
    {
    case 1:
      pD->ox1_BCFun = reflect_ox1;
      break;

    case 2:			/* Outflow */
      pD->ox1_BCFun = outflow_ox1;
      break;

    case 4:			/* Periodic. Handle with MPI calls for parallel jobs. */
      pD->ox1_BCFun = periodic_ox1;
      break;

    case 5:			/* Reflecting, B_normal!=0 */
      pD->ox1_BCFun = conduct_ox1;
      break;
/*
    default:
      printf (-1, "[bvals_init]:bc_ox1=%d unknown\n", pD->BCFlag_ox1);
      exit (EXIT_FAILURE);
*/
    }
}
/*---- ix2 boundary ----------------------------------------------------------*/
if (pG->Nx[1] > 1){

  switch (BCFlag_ix2)
    {
    case 1:
      pD->ix2_BCFun = reflect_ix2;
      break;

    case 2:			/* Outflow */
      pD->ix2_BCFun = outflow_ix2;
      break;

    case 4:			/* Periodic. Handle with MPI calls for parallel jobs. */
      pD->ix2_BCFun = periodic_ix2;
      break;

    case 5:			/* Reflecting, B_normal!=0 */
      pD->ix2_BCFun = conduct_ix2;
      break;
/*
    default:
      printf (-1, "[bvals_init]:bc_ix2=%d unknown\n", pD->BCFlag_ix2);
      exit (EXIT_FAILURE);
*/
   }

/*---- ox2 boundary ----------------------------------------------------------*/

  switch (BCFlag_ox2)
    {
    case 1:
      pD->ox2_BCFun = reflect_ox2;
      break;

    case 2:			/* Outflow */
      pD->ox2_BCFun = outflow_ox2;
      break;

    case 4:			/* Periodic. Handle with MPI calls for parallel jobs. */
      pD->ox2_BCFun = periodic_ox2;
      break;

    case 5:			/* Reflecting, B_normal!=0 */
      pD->ox2_BCFun = conduct_ox2;
      break;
/*
    default:
      printf (-1, "[bvals_init]:bc_ox2=%d unknown\n", pD->BCFlag_ox2);
      exit (EXIT_FAILURE);
*/
    }
  }
/*---- ix3 boundary ----------------------------------------------------------*/

if (pG->Nx[2] > 1){

  switch (BCFlag_ix3)
    {
    case 1:
      pD->ix3_BCFun = reflect_ix3;
      break;

    case 2:			/* Outflow */
      pD->ix3_BCFun = outflow_ix3;
      break;

    case 4:			/* Periodic. Handle with MPI calls for parallel jobs. */
      pD->ix3_BCFun = periodic_ix3;
      break;

    case 5:			/* Reflecting, B_normal!=0 */
      pD->ix3_BCFun = conduct_ix3;
      break;
/*
    default:
      printf (-1, "[bvals_init]:bc_ix3=%d unknown\n", pD->BCFlag_ix3);
      exit (EXIT_FAILURE);
*/
    }

/*---- ox3 boundary ----------------------------------------------------------*/

  switch (BCFlag_ox3)
    {
    case 1:
      pD->ox3_BCFun = reflect_ox3;
      break;

    case 2:			/* Outflow */
      pD->ox3_BCFun = outflow_ox3;
      break;

    case 4:			/* Periodic. Handle with MPI calls for parallel jobs. */
      pD->ox3_BCFun = periodic_ox3;
      break;

    case 5:			/* Reflecting, B_normal!=0 */
      pD->ox3_BCFun = conduct_ox3;
      break;
/*
    default:
      printf (-1, "[bvals_init]:bc_ox3=%d unknown\n", pD->BCFlag_ox3);
      exit (EXIT_FAILURE);
*/
    }
  } 
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   reflecting_???:   where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 *   outflow_???
 *   periodic_???
 *   conduct_???
 *   pack_???
 *   unpack_???
 */

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix1(GridS *pG)
 *  \brief  REFLECTING boundary conditions, Inner x1 boundary (bc_ix1=1) */

static void reflect_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][is-i]    =  pG->U[k][j][is+(i-1)];
        pG->U[k][j][is-i].M1 = -pG->U[k][j][is-i].M1; /* reflect 1-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox1(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x1 boundary (bc_ox1=1). */

static void reflect_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][ie+i]    =  pG->U[k][j][ie-(i-1)];
        pG->U[k][j][ie+i].M1 = -pG->U[k][j][ie+i].M1; /* reflect 1-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix2(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x2 boundary (bc_ix2=1) */

static void reflect_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][js-j][i]    =  pG->U[k][js+(j-1)][i];
        pG->U[k][js-j][i].M2 = -pG->U[k][js-j][i].M2; /* reflect 2-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox2(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x2 boundary (bc_ox2=1) */

static void reflect_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][je+j][i]    =  pG->U[k][je-(j-1)][i];
        pG->U[k][je+j][i].M2 = -pG->U[k][je+j][i].M2; /* reflect 2-mom. */
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix3(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x3 boundary (bc_ix3=1) */

static void reflect_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ks-k][j][i]    =  pG->U[ks+(k-1)][j][i];
        pG->U[ks-k][j][i].M3 = -pG->U[ks-k][j][i].M3; /* reflect 3-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox3(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x3 boundary (bc_ox3=1) */

static void reflect_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ke = pG->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ke+k][j][i]    =  pG->U[ke-(k-1)][j][i];
        pG->U[ke+k][j][i].M3 = -pG->U[ke+k][j][i].M3; /* reflect 3-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix1(GridS *pG)
 *  \brief OUTFLOW boundary condition, Inner x1 boundary (bc_ix1=2) */

static void outflow_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][is-i] = pG->U[k][j][is];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox1(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Outer x1 boundary (bc_ox1=2) */

static void outflow_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][ie+i] = pG->U[k][j][ie];
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix2(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Inner x2 boundary (bc_ix2=2) */

static void outflow_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][js-j][i] = pG->U[k][js][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox2(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Outer x2 boundary (bc_ox2=2) */

static void outflow_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][je+j][i] = pG->U[k][je][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix3(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Inner x3 boundary (bc_ix3=2) */

static void outflow_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ks-k][j][i] = pG->U[ks][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox3(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Outer x3 boundary (bc_ox3=2) */

static void outflow_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ke = pG->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ke+k][j][i] = pG->U[ke][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix1(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x1 boundary (bc_ix1=4) */

static void periodic_ix1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][is-i] = pG->U[k][j][ie-(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox1(GridS *pG)
 *  \brief PERIODIC boundary conditions (cont), Outer x1 boundary (bc_ox1=4) */

static void periodic_ox1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][ie+i] = pG->U[k][j][is+(i-1)];
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix2(GridS *pG)
 *  \brief PERIODIC boundary conditions (cont), Inner x2 boundary (bc_ix2=4) */

static void periodic_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][js-j][i] = pG->U[k][je-(j-1)][i];
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox2(GridS *pG)
 *  \brief PERIODIC boundary conditions (cont), Outer x2 boundary (bc_ox2=4) */

static void periodic_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][je+j][i] = pG->U[k][js+(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix3(GridS *pG)
 *  \brief PERIODIC boundary conditions (cont), Inner x3 boundary (bc_ix3=4) */

static void periodic_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ks-k][j][i] = pG->U[ke-(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox3(GridS *pG)
 *  \brief PERIODIC boundary conditions (cont), Outer x3 boundary (bc_ox3=4) */

static void periodic_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ke+k][j][i] = pG->U[ks+(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ix1(GridS *pG)
 *  \brief CONDUCTOR boundary conditions, Inner x1 boundary (bc_ix1=5) */

static void conduct_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][is-i]    =  pG->U[k][j][is+(i-1)];
        pG->U[k][j][is-i].M1 = -pG->U[k][j][is-i].M1; /* reflect 1-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ox1(GridS *pG)
 *  \brief CONDUCTOR boundary conditions, Outer x1 boundary (bc_ox1=5) */

static void conduct_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][ie+i]    =  pG->U[k][j][ie-(i-1)];
        pG->U[k][j][ie+i].M1 = -pG->U[k][j][ie+i].M1; /* reflect 1-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ix2(GridS *pG)
 *  \brief CONDUCTOR boundary conditions, Inner x2 boundary (bc_ix2=5) */

static void conduct_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][js-j][i]    =  pG->U[k][js+(j-1)][i];
        pG->U[k][js-j][i].M2 = -pG->U[k][js-j][i].M2; /* reflect 2-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ox2(GridS *pG)
 *  \brief CONDUCTOR boundary conditions, Outer x2 boundary (bc_ox2=5) */

static void conduct_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][je+j][i]    =  pG->U[k][je-(j-1)][i];
        pG->U[k][je+j][i].M2 = -pG->U[k][je+j][i].M2; /* reflect 2-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ix3(GridS *pG)
 *  \brief CONDUCTOR boundary conditions, Inner x3 boundary (bc_ix3=5) */

static void conduct_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ks-k][j][i]    =  pG->U[ks+(k-1)][j][i];
        pG->U[ks-k][j][i].M3 = -pG->U[ks-k][j][i].M3; /* reflect 3-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ox3(GridS *pG)
 *  \brief CONDUCTOR boundary conditions, Outer x3 boundary (bc_ox3=5) */

static void conduct_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ke = pG->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ke+k][j][i]    =  pG->U[ke-(k-1)][j][i];
        pG->U[ke+k][j][i].M3 = -pG->U[ke+k][j][i].M3; /* reflect 3-mom. */
      }
    }
  }

  return;
}
