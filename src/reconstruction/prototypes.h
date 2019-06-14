#ifndef RECONSTRUCTION_PROTOTYPES_H
#define RECONSTRUCTION_PROTOTYPES_H

#include <stdio.h>
#include "../struct.h"
#include "../defs.h"


/* lr_states_prim2.c */
void lr_states_destruct(void);
void lr_states_init(DomainS *);
void lr_states(const GridS *, const Prim1DS *, const double *, const double, 
               const double,  const int, const int, Prim1DS *, Prim1DS *, const int);

#endif
