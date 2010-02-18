/*
  $Id$

  xdr functions to write the 3D grid

  Copyright (C) 2004 Oliver Beckstein <oliver@biop.ox.ac.uk>
  This program is made available under the terms of the GNU Public License. 
  See the file LICENSE or http://www.gnu.org/copyleft/gpl.html

*/

#include "xdr_grid.h"

char *eGridType_names[egtyNR+1] = {"regular","other"};
char *version_ok_names[verokNR+1] = {"INCOMPATIBLE","OLDER","OK"};

bool xdr_grid_v1 (XDR *xdrs, t_XDRgrid *xg);
bool xdr_grid_v2 (XDR *xdrs, t_XDRgrid *xg);


#ifdef DOUBLE
#define xdr_real xdr_double
#else
#define xdr_real xdr_float 
#endif


enum version_ok version_check (int version) {
  /* implement version checks here, return INCOMPATIBLE if necessary */
  //if (version < GRID_FF_VERSION) return OLDER;
  if (version < GRID_FF_VERSION) return INCOMPATIBLE; // until v2 reads v1
  if (version > GRID_FF_VERSION) return INCOMPATIBLE;
  return OK;
}
    
  
/*
  Simple grid format as described in grid3D.h
  (Note: all xdr_* return 1 on success, 0 on error)
*/

bool xdr_grid (XDR *xdrs, t_XDRgrid *xg) {
  /* This is not working properly yet because the version check is
     actually happening while reading the file. For the time being,
     just run it with version = 2 for testing and figure out how to do
     it later. Perhaps 2-pass read or something clever with version
     check (chain in alternative xdr_grid_vX for the rest? */
  
  int version = GRID_FF_VERSION;

  /* For compatibility with older formats have versioned reading routines.
     The older ones fill in the missing bits with defaults.
   */

  switch (version) {
  case 1:
    return xdr_grid_v1(xdrs, xg);
    break;
  case 2:
    return xdr_grid_v2(xdrs, xg);
    break;
  default:
    gmx_fatal(FARGS,"xdr_grid(): unknown version %d of the grid file format",version);
  }
}


bool xdr_grid_v1 (XDR *xdrs, t_XDRgrid *xg) {
  int ngrid;
  
  ngrid = xg->size[0] * xg->size[1] * xg->size[2];
  xg->tweight = 1.0;    /* tweight default for format < 2 */
  return (xdr_int    (xdrs, &xg->version)           &&
	  (version_check(xg->version) > INCOMPATIBLE) &&
	  xdr_string (xdrs, &xg->header, HEADER_MAX) &&
	  xdr_enum   (xdrs,(enum_t *) &xg->type)    &&
	  xdr_int    (xdrs, &xg->dim)               &&
	  xdr_vector (xdrs, (char *) xg->size,   
		      xg->dim, sizeof(int), (xdrproc_t) xdr_int)   &&
	  xdr_vector (xdrs, (char *) xg->delta,  xg->dim, 
		      sizeof(real), (xdrproc_t)xdr_real)   &&
	  xdr_vector (xdrs, (char *) xg->origin, xg->dim, 
		      sizeof(real), (xdrproc_t)xdr_real)   &&
	  xdr_array  (xdrs, (char **) &xg->grid, &ngrid, NGRID_MAX,
		      sizeof(real), (xdrproc_t)xdr_real) );
} 

bool xdr_grid_v2 (XDR *xdrs, t_XDRgrid *xg) {
  int ngrid;
  
  ngrid = xg->size[0] * xg->size[1] * xg->size[2];
  return (xdr_int    (xdrs, &xg->version)           &&
	  (version_check(xg->version) > INCOMPATIBLE) &&
	  xdr_string (xdrs, &xg->header, HEADER_MAX) &&
	  xdr_enum   (xdrs,(enum_t *) &xg->type)    &&
	  xdr_real   (xdrs, &xg->tweight)           &&
	  xdr_int    (xdrs, &xg->dim)               &&
	  xdr_vector (xdrs, (char *) xg->size,   
		      xg->dim, sizeof(int), (xdrproc_t) xdr_int)   &&
	  xdr_vector (xdrs, (char *) xg->delta,  xg->dim, 
		      sizeof(real), (xdrproc_t)xdr_real)   &&
	  xdr_vector (xdrs, (char *) xg->origin, xg->dim, 
		      sizeof(real), (xdrproc_t)xdr_real)   &&
	  xdr_array  (xdrs, (char **) &xg->grid, &ngrid, NGRID_MAX,
		      sizeof(real), (xdrproc_t)xdr_real) );
} 


/* 
   serialise a 3D array of given sizes mx[] 
*/
real *grid3_serialise(real ***g,int *mx) {
  real *v,*tmp;        /* serialised grid */
  int i,j,k;           /* indices of the 3D grid */ 
  
  snew(v,mx[0]*mx[1]*mx[2]);
  if (!v) return NULL;

  tmp = v;
  for(k=0;k<mx[2];k++)
    for(j=0;j<mx[1];j++)
      for(i=0;i<mx[0];i++)
	*(tmp++) = g[i][j][k];
  return v;
}

/* 
   serialise a 2D array of given sizes mx[] 
*/
real *grid2_serialise(real **g,int *mx) {
  real *v,*tmp;        /* serialised grid */
  int i,j;             /* indices of the 2D grid */ 
  
  snew(v,mx[0]*mx[1]);
  if (!v) return NULL;

  tmp = v;
  for(j=0;j<mx[1];j++)
    for(i=0;i<mx[0];i++)
      *(tmp++) = g[i][j];
  return v;
}


real ***grid3_unserialise(real *v,int *mx) {
  real ***g;           /* 3D grid */
  int i,j,k;           /* indices of the 3D grid */ 
  
  g=grid3_alloc(mx[0],mx[1],mx[2]);
  if (!g) return NULL;

  for(k=0;k<mx[2];k++)
    for(j=0;j<mx[1];j++)
      for(i=0;i<mx[0];i++)
	g[i][j][k] = *(v++);
  return g;
}


real **grid2_unserialise(real *v,int *mx) {
  real **g;            /* 2D grid */
  int i,j;             /* indices of the 3D grid */ 
  
  g=grid2_alloc(mx[0],mx[1]);
  if (!g) return NULL;

  for(j=0;j<mx[1];j++)
    for(i=0;i<mx[0];i++)
      g[i][j] = *(v++);
  return g;
}
