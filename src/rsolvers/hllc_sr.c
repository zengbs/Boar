/*============================================================================*/
/*! \file hllc_sr.c
 *  \brief Computes 1D fluxes using the relativistic HLLC Riemann solver. 
 *
 * PURPOSE: Computes 1D fluxes using the relativistic HLLC Riemann solver, 
 *   an extension of the HLLE fluxes to include the contact wave.  Currently 
 *   only works for hydrodynamics.  For an extension to MHD, see hlld_sr.c
 *
 * REFERENCES:
 * - A. Mignone and G. Bodo, "An HLLC Riemann solver for relativistic flows",
 *   Mon. Not. R. Astron. Soc. 364, 126-136 (2005)
 *
 * - A. Mignone and G. Bodo, "An HLLC Solver for Relativistic Flows - II:
 *   Magnetohydrodynamics", arxiv:astro-ph/0601640v1 (2006)
 *
 * HISTORY: Written by Jonathan FUlton, February 2009
 *          Extended to MHD by Kris Beckwith, Spring 2010
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - fluxes() - all Riemann solvers in Athena must have this function name and
 *              use the same argument list as defined in rsolvers/prototypes.h*/
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../struct.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *      const Prim1DS Wl, const Prim1DS Wr, const double Bxi, Cons1DS *pFlux)
 *  \brief Computes 1D fluxes
 *   Input Arguments:
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface 
 *   - Wl,Wr = L/R-states of PRIMITIVE variables at cell interface 
 *   Output Arguments:
 *   - pFlux = pointer to fluxes of CONSERVED variables at cell interface 
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, Cons1DS *pFlux)
{
  Cons1DS Fl,Fr,Fhll,Uhll,Usl,Usr;
  double rhl, rhr, csl, csr, cslsq, csrsq, vsql, vsqr, gammasql, gammasqr;
  double ssl, ssr, radl, radr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
  double lmdal,lmdar; /* Left and Right wave speeds */
  double lmdas; /* Contact wave speed */
  double ovlrmll;
  double a,b,c,quad;
  double den,ps; /* PressUre in inner region */
  double lVx, lVy, lVz, rVx, rVy, rVz;
  double lFactor,rFactor;
 
/*--- Step 0. ------------------------------------------------------------------
 * Transform 4-velocity to 3-velocity
 */
  lFactor = 1.0/sqrt(1+Wl.Ux*Wl.Ux+Wl.Uy*Wl.Uy+Wl.Uz*Wl.Uz);
  rFactor = 1.0/sqrt(1+Wr.Ux*Wr.Ux+Wr.Uy*Wr.Uy+Wr.Uz*Wr.Uz);

  lVx = Wl.Ux * lFactor;
  lVy = Wl.Uy * lFactor;
  lVz = Wl.Uz * lFactor;

  rVx = Wr.Ux * rFactor;
  rVy = Wr.Uy * rFactor;
  rVz = Wr.Uz * rFactor;

/*--- Step 1. ------------------------------------------------------------------
 * Compute the max and min wave speeds used in Mignone 
 */
  rhl = Wl.d + Wl.P * Gamma / Gamma_1; /* Mignone Eq 3.5 */
  rhr = Wr.d + Wr.P * Gamma / Gamma_1;

  csl = sqrt(Gamma * Wl.P / rhl); /* Mignone Eq 4 */
  csr = sqrt(Gamma * Wr.P / rhr);

  cslsq = SQR(csl);
  csrsq = SQR(csr);

  vsql = SQR(lVx) + SQR(lVy) + SQR(lVz);
  vsqr = SQR(rVx) + SQR(rVy) + SQR(rVz);

  gammasql = 1.0 / (1.0 - vsql);
  gammasqr = 1.0 / (1.0 - vsqr);

  ssl = cslsq / ( gammasql * (1.0 - cslsq) ); /* Mignone Eq 22.5 */
  ssr = csrsq / ( gammasqr * (1.0 - csrsq) );

  radl = sqrt( ssl*(1.0-SQR(lVx)+ssl) ); /* Mignone Eq 23 (radical part) */
  radr = sqrt( ssr*(1.0-SQR(rVx)+ssr) );


  lmdapl = (lVx + radl) / (1.0 + ssl); /* Mignone Eq 23 */
  lmdapr = (rVx + radr) / (1.0 + ssr);
  lmdaml = (lVx - radl) / (1.0 + ssl);
  lmdamr = (rVx - radr) / (1.0 + ssr);


  lmdal = MIN(lmdaml, lmdamr); /* Mignone Eq 21 */
  lmdar = MAX(lmdapl, lmdapr);
  

/*--- Step 2. ------------------------------------------------------------------
 * Compute L/R fluxes according to Mignone 2
 */

  Fl.d  = Ul.d * lVx;
  Fl.Mx = Ul.Mx * lVx + Wl.P;
  Fl.My = Ul.My * lVx;
  Fl.Mz = Ul.Mz * lVx;
  Fl.E  = Ul.Mx;

  Fr.d  = Ur.d * rVx;
  Fr.Mx = Ur.Mx * rVx + Wr.P;
  Fr.My = Ur.My * rVx;
  Fr.Mz = Ur.Mz * rVx;
  Fr.E  = Ur.Mx;

/*--- Step 3. ------------------------------------------------------------------
 * Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
 * Compute HLL conserved quantities using Mignone eq 9
 */

  ovlrmll = 1.0 / ( lmdar - lmdal );
  lmdatlmda = lmdal*lmdar;

  Fhll.d  = (lmdar*Fl.d  - lmdal*Fr.d  + lmdatlmda * (Ur.d  - Ul.d) ) * ovlrmll;
  Fhll.Mx = (lmdar*Fl.Mx - lmdal*Fr.Mx + lmdatlmda * (Ur.Mx - Ul.Mx)) * ovlrmll;
  Fhll.My = (lmdar*Fl.My - lmdal*Fr.My + lmdatlmda * (Ur.My - Ul.My)) * ovlrmll;
  Fhll.Mz = (lmdar*Fl.Mz - lmdal*Fr.Mz + lmdatlmda * (Ur.Mz - Ul.Mz)) * ovlrmll;
  Fhll.E  = (lmdar*Fl.E  - lmdal*Fr.E  + lmdatlmda * (Ur.E  - Ul.E )) * ovlrmll;

  Uhll.d  = (lmdar * Ur.d  - lmdal * Ul.d  + Fl.d  - Fr.d ) * ovlrmll;
  Uhll.Mx = (lmdar * Ur.Mx - lmdal * Ul.Mx + Fl.Mx - Fr.Mx) * ovlrmll;
  Uhll.My = (lmdar * Ur.My - lmdal * Ul.My + Fl.My - Fr.My) * ovlrmll;
  Uhll.Mz = (lmdar * Ur.Mz - lmdal * Ul.Mz + Fl.Mz - Fr.Mz) * ovlrmll;
  Uhll.E  = (lmdar * Ur.E  - lmdal * Ul.E  + Fl.E  - Fr.E ) * ovlrmll;

/*--- Step 4. ------------------------------------------------------------------
 * Compute contact wave speed using larger root from Mignone Eq 18
 * Physical root is the root with the minus sign
 */

  /* quadratic formUla calcUlation */

  a = Fhll.E;
  b = -(Uhll.E + Fhll.Mx);
  c = Uhll.Mx;


  quad = -0.5*(b + SIGN(b)*sqrt(b*b - 4.0*a*c));
  lmdas = c/quad;

/*--- Step 5. ------------------------------------------------------------------
 * Determine intercell flux according to Mignone 13
 */
  if( lmdal >= 0.0){ /* Fl */
    /* intercell flux is left flux */
    pFlux->d  = Fl.d;
    pFlux->Mx = Fl.Mx;
    pFlux->My = Fl.My;
    pFlux->Mz = Fl.Mz;
    pFlux->E  = Fl.E;

   return;
  }
  else if( lmdas >= 0.0){ /* Fls */

    /* Mignone 2006 Eq 48 */
    ps = -Fhll.E*lmdas + Fhll.Mx;

    /* now calcUlate Usl with Mignone Eq 16 */
    den = 1.0 / (lmdal - lmdas);

    Usl.d  =  Ul.d  * (lmdal - lVx) * den;
    Usl.Mx = (Ul.Mx * (lmdal - lVx) + ps - Wl.P) * den;
    Usl.My =  Ul.My * (lmdal - lVx) * den;
    Usl.Mz =  Ul.Mz * (lmdal - lVx) * den;
    Usl.E  = (Ul.E  * (lmdal - lVx) + ps * lmdas - Wl.P * lVx) * den;

    /* now calcUlate Fsr using Mignone Eq 14 */

    pFlux->d  = lmdal*(Usl.d  - Ul.d ) + Fl.d;
    pFlux->Mx = lmdal*(Usl.Mx - Ul.Mx) + Fl.Mx;
    pFlux->My = lmdal*(Usl.My - Ul.My) + Fl.My;
    pFlux->Mz = lmdal*(Usl.Mz - Ul.Mz) + Fl.Mz;
    pFlux->E  = lmdal*(Usl.E  - Ul.E ) + Fl.E;

    return;
  }
  else if( lmdar >= 0.0){ /* Frs */

    /* Mignone 2006 Eq 48 */
    ps = -Fhll.E*lmdas + Fhll.Mx;

    /* now calcUlate Usr with Mignone Eq 16 */
    den = 1.0 / (lmdar - lmdas);

    Usr.d  =  Ur.d *  (lmdar - rVx) * den;
    Usr.Mx = (Ur.Mx * (lmdar - rVx) + ps - Wr.P) * den;
    Usr.My =  Ur.My * (lmdar - rVx) * den;
    Usr.Mz =  Ur.Mz * (lmdar - rVx) * den;
    Usr.E  = (Ur.E *  (lmdar - rVx) + ps * lmdas - Wr.P * rVx) * den;

    /* now calcUlate Fsr using Mignone Eq 14 */

    pFlux->d  = lmdar*(Usr.d  - Ur.d ) + Fr.d;
    pFlux->Mx = lmdar*(Usr.Mx - Ur.Mx) + Fr.Mx;
    pFlux->My = lmdar*(Usr.My - Ur.My) + Fr.My;
    pFlux->Mz = lmdar*(Usr.Mz - Ur.Mz) + Fr.Mz;
    pFlux->E  = lmdar*(Usr.E  - Ur.E ) + Fr.E;

    return;
  }
  else{ /* Fr */
    /* intercell flux is right flux */
    pFlux->d  = Fr.d;
    pFlux->Mx = Fr.Mx;
    pFlux->My = Fr.My;
    pFlux->Mz = Fr.Mz;
    pFlux->E  = Fr.E;

    return;
  }
}
