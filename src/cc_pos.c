/*============================================================================*/
/*! \file cc_pos.c
 *  \brief Functions to compute (x1,x2,x3) positions of cells i,j,k.  
 *
 * PURPOSE: Functions to compute (x1,x2,x3) positions of cells i,j,k.  
 *   In a nested grid, each Grid structure is a patch in a larger computational
 *   domain (with the exception of the level0 grid).  The displacement of the
 *   origin of the Grid from the origin of the computational domain (level0 
 *   grid) is x1_{disp} = idisp*dx1.  Furthermore, the origin of the level0
 *   grid can be displaced by a distance x1_{0} from the origin of x1.  Thus,
 *   the x1 position of the center of cell i (x1_{cc,i}) in any level Grid is
 *            x1_{cc,i} = x1_{0} + ((i + idisp) + 0.5)*dx1
 *   Similarly for x2 and x3.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - cc_pos() - given i,j,k returns cell-centered x1,x2,x3
 * - fc_pos() - given i,j,k returns face-centered x1,x2,x3
 * - x1cc() - given i, returns cell-centered x1.
 * - x2cc() - given j, returns cell-centered x2.
 * - x3cc() - given k, returns cell-centered x3.
 * - celli() - given x, returns containing cell first index. 
 * - cellj() - given y, returns containing cell first index.  
 * - cellk() - given y, returns containing cell first index.		      */
/*============================================================================*/

#include "struct.h"
#include "defs.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void cc_pos(const GridS *pG, const int i, const int j,const int k,
 *                  double *px1, double *px2, double *px3)
 *  \brief given i,j,k returns cell-centered x1,x2,x3
 */
void cc_pos(const GridS *pG, const int i, const int j,const int k,
	    double *px1, double *px2, double *px3)
{
  *px1 = pG->MinX[0] + ((double)(i - pG->is) + 0.5)*pG->dx1;
  *px2 = pG->MinX[1] + ((double)(j - pG->js) + 0.5)*pG->dx2;
  *px3 = pG->MinX[2] + ((double)(k - pG->ks) + 0.5)*pG->dx3;
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void fc_pos(const GridS *pG, const int i, const int j,const int k,
 *                  double *px1, double *px2, double *px3)
 *  \brief given i,j,k returns face-centered x1,x2,x3
 */

void fc_pos(const GridS *pG, const int i, const int j,const int k,
	    double *px1, double *px2, double *px3)
{
  *px1 = pG->MinX[0] + ((double)(i - pG->is))*pG->dx1;
  *px2 = pG->MinX[1] + ((double)(j - pG->js))*pG->dx2;
  *px3 = pG->MinX[2] + ((double)(k - pG->ks))*pG->dx3;
  return;
}
