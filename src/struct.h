#ifndef STRUCT_H
#define STRUCT_H

/*! \struct Int3Vect
 *  \brief General 3-vectors of ints.
 */
typedef struct Int3Vect_s{
  int i, j, k;
}Int3Vect;

/*conservative variables*/
typedef struct Cons_s
{
  double d;
  double M1;
  double M2;
  double M3;
  double E;
} ConsS;

/*primitive variables*/
typedef struct Prim_s
{
  double d;
  double U1;
  double U2;
  double U3;
  double P;
} PrimS;

typedef struct Cons1D_s
{
  double d;
  double Mx;
  double My;
  double Mz;
  double E;
} Cons1DS;


typedef struct Prim1D_s
{
  double d;
  double Ux;
  double Uy;
  double Uz;
  double P;
} Prim1DS;


typedef struct Grid_s
{
  ConsS ***U;			/*conserved variables */
  double MinX[3];		/*man(x) in each dir on this Grid [0,1,2]=[x1,x2,x3] */
  double MaxX[3];		/*min(x) in each dir on this Grid [0,1,2]=[x1,x2,x3] */
  double dx1, dx2, dx3;		/*cell size on this Grid */
  int is, ie;			/*star/end cell index in x1 direction */
  int js, je;			/*star/end cell index in x2 direction */
  int ks, ke;			/*star/end cell index in x3 direction */
  int Nx[3];			/*# of zones in each dir on Grid */
  double time, dt;		/*current time and time step */
} GridS;

typedef void (*VGFun_t) (GridS * pG);	/* generic void function of Grid */

typedef struct Domain_s
{
  GridS *Grid;
  VGFun_t ix1_BCFun, ox1_BCFun;	/*!< ix1/ox1 BC function pointers for grid */
  VGFun_t ix2_BCFun, ox2_BCFun;	/*!< ix1/ox1 BC function pointers for grid */
  VGFun_t ix3_BCFun, ox3_BCFun;	/*!< ix1/ox1 BC function pointers for grid */
  int BCFlag_ix1, BCFlag_ox1;  /*!< BC flag on root domain for inner/outer x1 */
  int BCFlag_ix2, BCFlag_ox2;  /*!< BC flag on root domain for inner/outer x2 */
  int BCFlag_ix3, BCFlag_ox3;  /*!< BC flag on root domain for inner/outer x3 */
  double time, dt;		/*current time and time step */
  int nstep;
}DomainS;

typedef void (*VDFun_t) (DomainS * pD);	/* generic void function of Grid */

#endif
