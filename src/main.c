#define MAIN_C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include "defs.h"
#include "struct.h"
#include "globals.h"
#include "prototypes.h"

int
main (void)
{
  int nz[3];
  double tlim = 0.5;
  VDFun_t Integrate;
  DomainS Domain;
  GridS Grid;
  Domain.Grid = &Grid;

  Domain.time = 0.0; // initial time
  Domain.dt = 2.3437500e-03; // time step
  Domain.nstep = 0;  // initial step

// number of cells at base level along x/y/z
  nz[0] = 128;
  nz[1] = 128;
  nz[2] = 128;
 
  Gamma = 1.4;
  Gamma_1 = Gamma - 1;
  Gamma_2 = Gamma - 2;

  init_grid (&Domain, nz);
  problem (&Domain);
  bvals_hydro_init (&Domain);
  bvals_hydro (&Domain);
  Integrate = integrate_init(&Domain);
  lr_states_init (&Domain);
  init_output (&Domain);

  while (Domain.time <= tlim)
    {
      Domain.nstep++;
      Domain.time += Domain.dt;
      (*Integrate)(&Domain);
      bvals_hydro (&Domain);
      output(&Domain);
    }

  lr_states_destruct();
  integrate_destruct();
  free_3d_array(Domain.Grid->U);

  return 0;
}
