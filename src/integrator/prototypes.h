#ifndef INTEGRATORS_PROTOTYPES_H
#define INTEGRATORS_PROTOTYPES_H

#include "stdio.h"
#include "../struct.h"
#include "../defs.h"

/*integrate.c*/
VDFun_t integrate_init(DomainS *pD);
void integrate_destruct(void);

void integrate_destruct_1d(void);
void integrate_destruct_2d(void);
void integrate_destruct_3d(void);
void integrate_init_1d(DomainS *pD);
void integrate_init_2d(DomainS *pD);
void integrate_init_3d(DomainS *pD);
void integrate_1d_vl(DomainS *pD);
void integrate_2d_vl(DomainS *pD);
void integrate_3d_vl(DomainS *pD);

#endif
