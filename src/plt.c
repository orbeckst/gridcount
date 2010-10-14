/* 
   $Id$

   writing density in plt format gOpenMol
   http://www.csc.fi/gopenmol/developers/plt_format.phtml

   Copyright (C) 2004 Oliver Beckstein <oliver@biop.ox.ac.uk>
   This program is made available under the terms of the GNU Public License. 
   See the file LICENSE or http://www.gnu.org/copyleft/gpl.html

   This file contains code contributed by Christoph Freudenberger from 
   his program g_sdf.c [g_sdf.c,v 1.25 2003/08/05 11:00:00 cfreuden]
   Copyright (C) 2003, 2004 Christoph Freudenberger:

   i_write() and f_write()

   density_write_plt() is based on code in g_sdf.c

*/
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "futil.h"
#include "statutil.h"
#include "grid3D.h"

static char *SRCID_plt_c = "$Id$";

static void i_write(FILE *output, int value)
{
  fwrite(&value,sizeof(int),1L,output);
}

static void f_write(FILE *output,float value)
{
  fwrite(&value,sizeof(float),1L,output);
}

/* write plt binary density file for gOpenMol
   Normalise each grid point by dividing by unit.
*/
void density_write_plt (char *fnm, t_tgrid *tgrid, real unit)
{
  FILE *fp;
  int  i,j,k;
  real norm;

  norm = 1./unit;
  fp = ffopen (fnm, "wb");

  /* rank */
  i_write(fp,3);

  /* Type of surface -- use 42 for gridcount */
  i_write(fp,42);

  /* Zdim, Ydim, Xdim */
  for (i=ZZ; i>=XX; i--)
    i_write(fp,tgrid->mx[i]);

  /* [Z,Y,X][min,max] (box corners in Angstroem)
     see grid3D.h for definition of t_grid (note: distances in nm)
  */
  for (i=ZZ; i>=XX; i--) {
    f_write(fp,10. * tgrid->a[i]);
    f_write(fp,10. * tgrid->b[i]);
  }

  for(k=0;k<tgrid->mx[ZZ];k++) {
    for(j=0;j<tgrid->mx[YY];j++) {
      for(i=0;i<tgrid->mx[XX];i++) {
        f_write(fp,norm * tgrid->grid[i][j][k]);
      }
    }
  }
  ffclose(fp);
};


void density_write_ascii (char *fnm, t_tgrid *tgrid, real unit)
{
  FILE *fp;
  int  i,j,k;
  real norm;

  norm = 1./unit;
  fp = ffopen (fnm, "w");

  for(k=0;k<tgrid->mx[ZZ];k++) {
    for(j=0;j<tgrid->mx[YY];j++) {
      for(i=0;i<tgrid->mx[XX];i++) {
        fprintf(fp,"%.6f ",norm * tgrid->grid[i][j][k]);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

