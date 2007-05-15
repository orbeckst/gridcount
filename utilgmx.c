/* 
   $Id$
   
   everyday utility functions

   Copyright (C) 2004 Oliver Beckstein <oliver@biop.ox.ac.uk>
   This program is made available under the terms of the GNU Public License. 
   See the file LICENSE or http://www.gnu.org/copyleft/gpl.html
*/

static char *SRCID_utilgmx_c = "$Id$";

#include "utilgmx.h"
#include "fatal.h"
#include "futil.h"
#include "maths.h"
#include "smalloc.h"

real ***grid3_alloc(const int nx,const int ny,const int nz) {
  /* 
  This is a trimmed down version of f3tensor() from
  Press et al, "Numerical Recipes in C", Appendix B, p944.
  The following 'copyright' applies to all routines from the Appendix:

  "Non-Copyright Notice: This Appendix and its utility routines are herewith
  placed into the public domain. Anyone may copy them freely for any purpose.
  We of course accept no liability whatsoever for any such use."

  */

  /* allocate a 3D array, size nx*ny*nz which is referenced through pointers.
     By using snew() the memory is already 0-initialised.
     range is g[0..nx-1,0..ny-1,0..nz-1].
     x: rows
     y: columns
     z: depth

     Basically, g[] is a vector, pointing to the rows of an array
     g[][], nx*ny (allocated as one block). These entries again point
     to rows in the 'real' chunk of memory nx*ny*nz, g[][][].  Because
     we have allocate the memory in blocks we can initialise all
     pointers in a fairly tricky fashion and also free them easily.
  */
  int  i,j;
  real ***g = NULL;

  /* allocate a row of 'ptr to ptr to rows' (= index by x values), length nx */
  snew(g,nx);
  if (!g) gmx_fatal(FARGS,"grid3_alloc(): Cannot allocate memory for nx=%d\n",nx);
  
  /* allocate nx 'ptrs to rows' (length ny each) and set ptrs to them */
  snew(g[0],nx*ny);
  if (!g) gmx_fatal(FARGS,"grid3_alloc(): Cannot allocate memory for nx*ny=%d*%d\n",
		      nx,ny);
    
  /* allocate nx*ny 'rows' (length nz each; these memory locations hold the reals) */
  snew(g[0][0],nx*ny*nz);
  if (!g) gmx_fatal(FARGS,"grid3_alloc(): Cannot allocate memory for data "
		      "(nx*ny*nz=%d*%d%d)\n",nx,ny,nz);

  /* initialise all pointers (starting from the already initialised
     initial locations), "leaping" through the array 
  */
  for(j=1; j<ny; j++) g[0][j]=g[0][j-1]+nz;
  for(i=1; i<nx; i++) {
    g[i]=g[i-1]+ny;
    g[i][0]=g[i-1][0]+ny*nz;
    for(j=1; j<ny; j++) g[i][j]=g[i][j-1]+nz;
  }

  /* return ptr to an array of ptr to rows */
  return g;
}


real **grid2_alloc(const int nx,const int ny) {
  /* 
  This is a butchered version of f3tensor() from
  Press et al, "Numerical Recipes in C", Appendix B, p944. 
  The following 'copyright' applies to all routines from the Appendix:

  "Non-Copyright Notice: This Appendix and its utility routines are herewith
  placed into the public domain. Anyone may copy them freely for any purpose.
  We of course accept no liability whatsoever for any such use."

  */

  /* allocate a 2D array, size nx*ny which is referenced through pointers.
     By using snew() the memory is already 0-initialised.
     range is g[0..nx-1,0..ny-1].
     x: rows
     y: columns
  */
  int  i;
  real **g = NULL;

  /* allocate a row of 'ptr to ptr to rows' (= index by x values), length nx */
  snew(g,nx);
  if (!g) gmx_fatal(FARGS,"grid2_alloc(): Cannot allocate memory for nx=%d\n",nx);
  
  /* allocate nx 'ptrs to rows' (length ny each) and set ptrs to them */
  snew(g[0],nx*ny);
  if (!g) gmx_fatal(FARGS,"grid2_alloc(): Cannot allocate memory for nx*ny=%d*%d\n",
		      nx,ny);

  for(i=1;i<nx;i++) g[i]=g[i-1]+ny;
  
  /* return ptr to an array of ptr to rows */
  return g;
}

void grid3_free(real ***g) {
  sfree(g[0][0]);
  sfree(g[0]);
  sfree(g);
}

void grid2_free(real **g) {
  sfree(g[0]);
  sfree(g);
}


void msg (const char *format, ...)
{
  va_list args;
  va_start(args, format);

  vfprintf (stderr, format, args);
  if (bDebugMode()) {
    vfprintf(debug, format, args);
  };

  va_end(args);
  return;
};

void dmsg (const char *format, ...)
{
  va_list args;
  va_start(args, format);

  if (bDebugMode()) {
    vfprintf(stderr, format, args);
    vfprintf(debug,  format, args);
  };

  va_end(args);
  return;
};
  

void dprintf (const char *format, ...)
{
  va_list args;
  va_start(args, format);
  
  if (bDebugMode()) {
    vfprintf(stderr, format, args);
  };
  va_end(args);
  return;
};

void dfprintf (const char *format, ...)
{
  /* write to debug file (see fatal.h) */
  va_list args;
  va_start(args, format);
  
  if (bDebugMode()) {
    vfprintf(debug, format, args);
  };
  va_end(args);
  return;
};


real ldist (rvec x, rvec p, rvec c) {
  /* calculate the (perpendicular) distance of point x from the line
     through c pointing along p 
  */
  real lambda;
  rvec d, u;

  /* u = x - c */
  rvec_sub (x, c, u);
  
  /* d = [(u.p)/p^2]p  - u */
  lambda = iprod (u, p) / iprod (p, p);
  svmul (lambda, p, d);
  rvec_dec (d, u);

  return norm (d);
};

real dt_tpx (char *fn) {
  int step, natoms;
  real t, lambda;
  t_inputrec  ir;  

  /* discard almost all info */
  read_tpx(fn,&step,&t,&lambda, &ir, NULL, &natoms, NULL, NULL, NULL, NULL);

  return ir.delta_t * ir.nstxtcout;
};

int list_add_atomid (const atom_id id, int *nid, atom_id *list) {
  /* append id to list and increase nid if not already in there 
     return the new number of ids, nid, if any changes occurred, 0 otherwise
     
     ATTENTION: nid is modified in the caller!!
   */
  int j;
  for (j=0; j < *nid  &&  id != list[j]; j++);
  if (j >= *nid) {
    list[(*nid)++] = id;
    return *nid;
  };
  return 0;
};
