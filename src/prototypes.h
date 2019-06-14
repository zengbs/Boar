#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include "struct.h"
#include "defs.h"
#include "integrator/prototypes.h"
#include "reconstruction/prototypes.h"
#include "rsolvers/prototypes.h"

/*----------------------------------------------------------------------------*/
/* array.c */
void*   calloc_1d_array(                      size_t nc, size_t size);
void**  calloc_2d_array(           size_t nr, size_t nc, size_t size);
void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size);
void free_1d_array(void *array);
void free_2d_array(void *array);
void free_3d_array(void *array);

/*----------------------------------------------------------------------------*/
/* bvals_hydro.c  */
void bvals_hydro_init(DomainS *pD);
void bvals_hydro(DomainS *pD);

/*----------------------------------------------------------------------------*/
/* cc_pos.c */
void cc_pos(const GridS *pG, const int i, const int j,const int k,
            double *px1, double *px2, double *px3);

/*----------------------------------------------------------------------------*/
/* convert_var.c */
PrimS Cons_to_Prim(const ConsS *pU);
ConsS Prim_to_Cons(const PrimS *pW);
Prim1DS Cons1D_to_Prim1D(const Cons1DS *pU);
Cons1DS Prim1D_to_Cons1D(const Prim1DS *pW);
PrimS check_Prim(const ConsS *pU);
Prim1DS check_Prim1D (const Cons1DS *pU);

/*----------------------------------------------------------------------------*/
/* init_grid.c */
void init_grid(DomainS *pD, int *);

/*----------------------------------------------------------------------------*/
/* prob/PROBLEM.c ; linked to problem.c */
void problem(DomainS *pD);
/*
void Userwork_in_loop(MeshS *pM);
void Userwork_after_loop(MeshS *pM);
void problem_read_restart(MeshS *pM, FILE *fp);
void problem_write_restart(MeshS *pM, FILE *fp);
ConsFun_t get_usr_expr(const char *expr);
VOutFun_t get_usr_out_fun(const char *name);
*/

/*----------------------------------------------------------------------------*/
/* output.c */
void output(DomainS *pD);
void init_output(DomainS *pD);

#endif
