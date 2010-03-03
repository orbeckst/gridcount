/*
  $Id$

  xdr functions to write the grid

  Copyright (C) 2004 Oliver Beckstein <oliver@biop.ox.ac.uk>
  This program is made available under the terms of the GNU Public License. 
  See the file LICENSE or http://www.gnu.org/copyleft/gpl.html

*/

#ifndef _xdr_grid_h
#define _xdr_grid_h

static char *SRCID_xdr_grid_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include "gmx_fatal.h"
#include "names.h"
#include "smalloc.h"
#include "typedefs.h"
#include "utilgmx.h"


/* Version number of the file format is an integer, which is incremented 
   when changes are made. If the change breaks something then one has to add 
   specific tests on the format in version_check() and return INCOMPATIBLE 
*/
#define GRID_FF_VERSION 2      /* grid xdr file format version */
#define HEADER_MAX 256         /* descriptive header string */
#define NGRID_MAX  512*512*512 /* max size for a grid (so that I can use
				  xdr auto-allocate in xdr_array()*/

#ifdef DOUBLE
#define xdr_real xdr_double
#else
#define xdr_real xdr_float
#endif

/* these must correspond to eGridType_names[] */
enum eGridType {egtyREGULAR,egtyOTHER,egtyNR};
enum version_ok {INCOMPATIBLE,OLDER,OK,verokNR};

extern char *eGridType_names[egtyNR+1];
extern char *version_ok_names[verokNR+1];

typedef struct {
  int   version;
  char  *header;
  enum eGridType type;
  real  tweight;   /* total simulation time in ps (used as a weight for calculations) */
  int   dim;
  int   *size;
  real  *delta;
  real  *origin;   
  real  *grid;     /* serialised grid (one vector) because I do not
                      know how to put there a type which handles 2D
                      (**real, 3D ***real etc) */
} t_XDRgrid;

extern enum version_ok version_check(int);
extern bool xdr_grid (XDR *, t_XDRgrid *);
extern real *grid3_serialise(real ***,int *);
extern real *grid2_serialise(real  **,int *);
extern real ***grid3_unserialise(real *,int *);
extern real  **grid2_unserialise(real *,int *);

#define GRIDTYPE(e)    ENUM_NAME((e),egtyNR,eGridType_names)
#define VERSION_OK(e)  ENUM_NAME((e),verokNR,version_ok_names)


#endif /* _xdr_grid_h */
