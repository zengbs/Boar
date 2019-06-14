#ifndef GLOBALS_H
#define GLOBALS_H  
/*============================================================================*/
/*! \file globals.h
 *  \brief Contains global variables.
 *
 * PURPOSE: Contains global variables:
 *   The first occurence in this file is included in main.c and defines the
 *   variables.  The second is included everywhere else.		      */
/*============================================================================*/

#ifdef MAIN_C

double CourNo;                 /*!< Courant, Friedrichs, & Lewy (CFL) number */
double Gamma;                  /*!< adiabatic index (ratio of specific heats) */
double Gamma_1, Gamma_2;       /*!< (Gamma)-1 and (Gamma)-2 */

/*----------------------------------------------------------------------------*/
/* definitions included everywhere except main.c  */

#else /* MAIN_C */

extern double CourNo;
extern double Gamma, Gamma_1, Gamma_2;
extern double d_MIN;
extern double etah;
#endif /* MAIN_C */
#endif /* GLOBALS_H */
