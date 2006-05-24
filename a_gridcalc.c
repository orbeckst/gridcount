/*
 * $Id$
 * analyse grid densities produced by g_ri3Dc
 
   Copyright (C) 2003, 2004 Oliver Beckstein <oliver@biop.ox.ac.uk>
   This program is made available under the terms of the GNU Public License. 
   See the file LICENSE or http://www.gnu.org/copyleft/gpl.html

 */
static char *SRCID_a_ri3Dc_c = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <assert.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "xvgr.h"
#include "gstat.h"
#include "names.h"
#include "filenm.h"
#include "gridcount.h"

/* a_gridcalc
   ----------

   Purpose: 

   manipulte occupation maps (3D grids) produced by g_ri3Dc

   Load a numbar of density grids and combine them with algebraic
   operations, or average them.

   The grids

   - must be of the same shape (dimensions, number of cells)

   - should be file format >= 2 because from then on we store the
     simulation time that was ran to produce the data (Tweight). This
     time is used as weight when averaging densities.
   

*/


void update_tgrid (t_tgrid *);

/* complete tgrid from data hold there */
void update_tgrid (t_tgrid *tg) {
  int i;
  /* top right corner */
  for(i=0;i<DIM;i++)
    tg->b[i] = tg->a[i] + tg->mx[i] * tg->Delta[i];
  return;
}

    
    
int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]a_gridcalc[TT] combines  3D grids produced by [TT]g_ri3Dc[TT] [1]. "
    "[PAR]Known limitations:"
  };

  static char *bugs[] = {
    "requires grid files with format >= 2 (including Tweight field)",
    "only does the average, no other operations implemented"
  };

  static char buf[HEADER_MAX];           /* additional text for graphs */
  static char *header = buf;



  t_pargs pa[] = {
    { "-subtitle", FALSE, etSTR, {&header},
      "Some text to add to the output graphs"},
  };

  t_filenm fnm[] = {
    { efDAT, "-f", "gridxdr",      ffRDMULT },
    { efDAT, "-o", "gridout", ffWRITE },
  };

  int     nfile;          /* number of input files */
  char    **fnms;         /* filenames of input files */
  FILE    *fGrid;         /* 3D grid with occupation numbers */
  FILE    *fOut;          /* reusable fp */
  t_tgrid *tgrids;        /* all information about the input grids */

  int i,j,k, ifile;
  double dV;        /* volume of a cell */
  real sumT;        /* total time for data (ps) */
  real volume = 0;

  char s_tmp[STRLEN];        /* utility buffer, sprintf on the fly */


#define NFILE asize(fnm)
#define NPA   asize(pa)

  strncpy(header,"",HEADER_MAX);

  CopyRight(stderr,argv[0]);
  fprintf (stderr, "\nVersion: %s\n\n", SRCID_a_ri3Dc_c);

  parse_common_args(&argc,argv, 0,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);

  if (bDebugMode()) {
    dfprintf ("%s -- debugging...\n\n", Program());
    dfprintf ("Version of the grid format (software): %d\n",GRID_FF_VERSION);
  };


  /* read in multiple xdr files */
  nfile = opt2fns(&fnms,"-f",NFILE,fnm);
  if(!nfile)
      fatal_error(0,"No input files!");

  /* alloc array of tgrids; grid number nfile is the result */
  snew(tgrids,nfile+1);

  for(ifile=0; (ifile<nfile); ifile++) {
    /* open input file */
    msg("Reading grid 3D file %s\n", fnms[ifile]);
    fGrid    = ffopen (fnms[ifile], "r");
    if (!grid_read(fGrid,&tgrids[ifile],header)) 
      fatal_error(0,"Error reading the 3D grid---no point in continuing!\n");
    fclose(fGrid);
    update_tgrid(&tgrids[ifile]);

    msg("Tweight = %g ps\n",tgrids[ifile].tweight);
    msg("mx = %d %d %d\n",tgrids[ifile].mx[XX],tgrids[ifile].mx[YY],tgrids[ifile].mx[ZZ]);

    sumT += tgrids[ifile].tweight;
    dmsg("sumT[%d] = %g\n",ifile,sumT);
    msg(".. done!\n");
  }

  /* sanity checks would be nice, at least for grid dimensions ... */
  msg("No sanity checks... proceeding with the above %d grids.\n",nfile);
  
  t_tgrid *a;          /* result: averaged grid */
  t_tgrid *g;          /* first grid, used as template */
  a = &tgrids[nfile];
  g = &tgrids[0];

  for(i=XX;i<DIM;i++) {
    a->mx[i] = g->mx[i];
    a->a[i]  = g->a[i];
    a->b[i]  = g->b[i];
  }
  copy_rvec(g->Delta,a->Delta);
  a->tweight = sumT;   /* store total time in final result */
  a->grid = grid3_alloc(a->mx[XX],a->mx[YY],a->mx[ZZ]);

  /* quick results are needed: hard code the averaging */
  /* The weight for each data point from this grid is computed and
     stored IN PLACE of the absolute time value.  By doing the
     normalisation a priori we hope to get around overflow problems.
  */

  for(ifile=0;ifile<nfile;ifile++) {
    g = &tgrids[ifile];
    g->tweight = g->tweight/sumT;    /* ATTENTION: changed the input data */
    for(k=0;k<g->mx[ZZ];k++) 
      for(j=0;j<g->mx[YY];j++) 
	for(i=0;i<g->mx[XX];i++) {
	  a->grid[i][j][k] += g->grid[i][j][k] * g->tweight;
	}
  }
  
  /* done */
  snprintf (header, HEADER_MAX,
	    "Average density: {Ttot}%gps{Ngrid}%d{<Delta>}%g"
	    "{mx,my,mz}%d,%d,%d",
	    sumT,nfile,(a->Delta[XX]+a->Delta[YY]+a->Delta[ZZ])/3,
	    a->mx[XX],a->mx[YY],a->mx[ZZ]
	    );
  header[HEADER_MAX-1] = '\0';   /* better safe than sorry */

  fOut = ffopen(opt2fn("-o",NFILE,fnm),"w");
  if (!grid_write(fOut,a,header))
    msg("Writing the 3D grid failed.\n");

  msg("Wrote averaged density from %g ps of density data.\n",sumT);

  for(i=0;i<=nfile;i++) free_tgrid(&tgrids[i]);

  fprintf(stderr,"\n");
  thanx(stdout);

  return 0;
}

