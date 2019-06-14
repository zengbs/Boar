#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "struct.h"
#include "defs.h"
#include "globals.h"
#include "prototypes.h"

static char filePath[100];
static void output_1d (DomainS *, FILE *);
static void output_2d (DomainS *, FILE *);
static void output_3d (DomainS *, FILE *);
static int dim = 0;

void
init_output (DomainS * pD)
{
  FILE *fptr[1];
  char fileName[100] = "";
  char file[200] = "";

  sprintf (filePath, "./");
  sprintf (fileName, "%06d.dat", pD->nstep);

  strcat (file, filePath);
  strcat (file, fileName);

  fptr[0] = fopen (file, "w");
  if (fptr[0] == NULL)
    {
      printf ("fail to create data file!\n");
      abort ();
    }
  dim = 0;
  for (int i = 0; i < 3; i++)
    if (pD->Grid->Nx[i] > 1)
      dim++;

  switch (dim)
    {
    case 1:
      if (pD->Grid->Nx[0] <= 1)
	break;
      output_1d (pD, fptr[0]);
      break;

    case 2:
      if (pD->Grid->Nx[2] > 1)
	break;
      output_2d (pD, fptr[0]);
      break;

    case 3:
      output_3d (pD, fptr[0]);
      break;
    }
  fclose (fptr[0]);
}

void
output (DomainS * pD)
{
  FILE *fptr[100000];
  char fileName[100] = "";
  char file[200] = "";

  sprintf (fileName, "%06d.dat", pD->nstep);

  strcat (file, filePath);
  strcat (file, fileName);

  fptr[pD->nstep] = fopen (file, "w");
  if (fptr[pD->nstep] == NULL)
    {
      printf ("fail to create data file!\n");
      abort ();
    }

  dim = 0;
  for (int i = 0; i < 3; i++)
    if (pD->Grid->Nx[i] > 1)
      dim++;
  switch (dim)
    {
    case 1:
      if (pD->Grid->Nx[0] <= 1)
	break;
      output_1d (pD, fptr[pD->nstep]);
      break;
    case 2:
      if (pD->Grid->Nx[2] > 1)
	break;
      output_2d (pD, fptr[pD->nstep]);
      break;

    case 3:
      output_3d (pD, fptr[pD->nstep]);
      break;
    }

  fclose (fptr[pD->nstep]);
}

static void
output_1d (DomainS * pD, FILE * fptr)
{
  GridS *pG = pD->Grid;

  Prim1DS W;
  Cons1DS U1d;

  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  double x1, x2, x3;
  double fac;
  double d, Vx, P, Lfac;

  fprintf (fptr, "#%10s %10s %10s %20s %20s %20s", "i[1]", "j[2]", "k[3]", "x[4]", "y[5]", "z[6]");
  fprintf (fptr, "%14s%14s%14s%14s%14s", "Dens[7]", "MomX[8]", "MomY[9]", "MomZ[10]", "Engy[11]");
  fprintf (fptr, "%14s%14s%14s%14s%14s%19s","PrimDens[12]", "Vx[13]", "Vy[14]", "Vz[15]", "Pressure[16]", "Lorentz Fac[17]\n");

  for (k = ks; k <= ke; k++)
    {
      for (j = js; j <= je; j++)
	{
	  for (i = is; i <= ie; i++)
	    {
	      cc_pos (pG, i, j, k, &x1, &x2, &x3);

	      U1d.d = pG->U[k][j][i].d;
	      U1d.Mx = pG->U[k][j][i].M1;
	      U1d.My = pG->U[k][j][i].M2;
	      U1d.Mz = pG->U[k][j][i].M3;
	      U1d.E = pG->U[k][j][i].E;

	      W = Cons1D_to_Prim1D (&U1d);
	      Lfac = sqrt (1 + SQR (W.Ux));
	      d = W.d;
	      fac = pow (1 + W.Ux * W.Ux, -0.5);
	      Vx = W.Ux * fac;
	      P = W.P;

	      fprintf (fptr, " %10d %10d %10d %20.14e %20.14e %20.14e", i, 0, 0, x1, 0.05, 0.05);
	      fprintf (fptr, " %13.6e %13.6e %13.6e %13.6e %13.6e", U1d.d, U1d.Mx, 0.0, 0.0, U1d.E);
	      fprintf (fptr, " %13.6e %13.6e %13.6e %13.6e %13.6e %17.6e\n", d, Vx, 0.0, 0.0, P, Lfac);
	    }
	}
    }
}

static void
output_2d (DomainS * pD, FILE * fptr)
{
  GridS *pG = pD->Grid;

  Prim1DS W;
  Cons1DS U1d;

  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  double x1, x2, x3;
  double fac;
  double d, Vx, Vy, P, Lfac;

  fprintf (fptr, "# i=[1], j=[2], x1=[3], x2=[4], d=[5], Vx=[6], Vy=[7], P=[8], D=[9], Mx=[10], My=[11], E=[12], Lorentz fac.=[13]\n");

  for (k = ks; k <= ke; k++)
    {
      for (j = js; j <= je; j++)
	{
	  for (i = is; i <= ie; i++)
	    {
	      cc_pos (pG, i, j, k, &x1, &x2, &x3);

	      U1d.d = pG->U[k][j][i].d;
	      U1d.Mx = pG->U[k][j][i].M1;
	      U1d.My = pG->U[k][j][i].M2;
	      U1d.Mz = pG->U[k][j][i].M3;
	      U1d.E = pG->U[k][j][i].E;

	      W = Cons1D_to_Prim1D (&U1d);
	      Lfac = sqrt (1 + SQR (W.Ux) + SQR (W.Uy));
	      d = W.d;
	      fac = pow (1 + W.Ux * W.Ux + W.Uy * W.Uy, -0.5);
	      Vx = W.Ux * fac;
	      Vy = W.Uy * fac;
	      P = W.P;

	      fprintf (fptr, "%3d %3d  %10.5e  %10.5e %10.5e  %10.5e  %10.5e  %10.5e %10.5e  %10.5e  %10.5e %10.5e %10.5e\n",
		       i, j, x1, x2, d, Vx, Vy, P, U1d.d, U1d.Mx, U1d.My, U1d.E, Lfac);
	    }
	}
    }
}

static void
output_3d (DomainS * pD, FILE * fptr)
{
  GridS *pG = pD->Grid;

  Prim1DS W;
  Cons1DS U1d;

  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  double x1, x2, x3;
  double fac;

  double d, Vx, Vy, Vz, P, Lfac;

  fprintf (fptr, "#%10s %10s %10s %20s %20s %20s", "i[1]", "j[2]", "k[3]", "x[4]", "y[5]", "z[6]");
  fprintf (fptr, "%14s%14s%14s%14s%14s", "Dens[7]", "MomX[8]", "MomY[9]", "MomZ[10]", "Engy[11]");
  fprintf (fptr, "%14s%14s%14s%14s%14s%19s","PrimDens[12]", "Vx[13]", "Vy[14]", "Vz[15]", "Pressure[16]", "Lorentz Fac[17]\n");

  for (k = ks; k <= ke; k++)
    {
      for (j = js; j <= je; j++)
	{
	  for (i = is; i <= ie; i++)
	    {
	      cc_pos (pG, i, j, k, &x1, &x2, &x3);

	      U1d.d = pG->U[k][j][i].d;
	      U1d.Mx = pG->U[k][j][i].M1;
	      U1d.My = pG->U[k][j][i].M2;
	      U1d.Mz = pG->U[k][j][i].M3;
	      U1d.E = pG->U[k][j][i].E;

	      W = Cons1D_to_Prim1D (&U1d);

	      Lfac = sqrt (1 + SQR (W.Ux) + SQR (W.Uy) + SQR (W.Uz));
	      d = W.d;
	      fac = pow (1 + W.Ux * W.Ux + W.Uy * W.Uy + W.Uz * W.Uz, -0.5);
	      Vx = W.Ux * fac;
	      Vy = W.Uy * fac;
	      Vz = W.Uz * fac;
	      P = W.P;

	      //if ((k == 4) && (j==4))
	      if (k==67)
		{
		  fprintf (fptr, " %10d %10d %10d %20.14e %20.14e %20.14e", i, j, k, x1, x2, x3);
		  fprintf (fptr, " %13.6e %13.6e %13.6e %13.6e %13.6e", U1d.d, U1d.Mx, U1d.My, U1d.Mz, U1d.E);
		  fprintf (fptr, " %13.6e %13.6e %13.6e %13.6e %13.6e %17.6e\n", d, Vx, Vy, Vz, P, Lfac);
		}

	    }
	}
    }
}
