/*
  $Id$

  collect stuff common to g_ri3Dc and a_ri3Dc

  Copyright (C) 2003, 2004 Oliver Beckstein <orbeckst@jhmi.edu>
  This program is made available under the terms of the GNU Public License. 
  See the file LICENSE or http://www.gnu.org/copyleft/gpl.html


  grid: use XDR routines (similar to the gromacs libxdrf.c but less elaborate)
        advantages:
	- binary format (more efficient than huge ascii files)
	- could implement compression similar to xdr3dfcoord()
	- machine independent

  Format: 
  
  VERSION   int              # version info (to catch the inevitable 
                               format change in the future)
  HEADER    str[HEADER_MAX]  # descriptive txt
  GRIDTYPE  enum             # currently we only have regular grids
  TWEIGHT   real             # length of simulation time (ps) for the data
  DIMENSION int              # dim of array
  SIZE      int int ...      # number of cells per dim  \
  DELTA     real real ...    # grid width per dim        } allows reconstruction 
  ORIGIN    real real ...    # lower left corner of grid/
  GRID      real ..................................
                 (0,0,0) (1,0,0) (2,0,0) ... (N1,0,0)       \
                 (0,1,0) (1,1,0) (2,1,0) ... (N1,1,0)       \
		 ...                                        \
                 (0,N2,0) (1,N2,0) (2,N2,0) ... (N1,N2,0)   \
                 (0,0,1) (1,0,1) (2,0,1) ... (N1,0,1)       \
                 (0,1,1) (1,1,1) (2,1,1) ... (N1,1,1)       \
		 ...                                        \
                 (0,N2,0) (1,N2,0) (2,N2,0) ... (N1,N2,0)   \
		 ... 
		 ...
		 .... (N1,N2,N3)

		 (ordered by z-slices)

*/


#ifndef _grid3D_h
#define _grid3D_h

static char *SRCID_grid3D_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "gmx_fatal.h"
#include "xdr_grid.h"
#include "count.h"
#include "utilgmx.h"

typedef struct {
  /* everything about the grid */
  real ***grid;    /* MX x MY x MZ */
  int  mx[3];      /* MX, MY, MZ */
  real a[3], b[3]; /* bottom left and top right corner of the grid in nm */ 
  rvec Delta;      /* spatial resolution in nm */
  real tweight;    /* length of simulation time (ps) for data; used as a */
                   /* weight for combining grids in grid_calc.c */
} t_tgrid;

typedef struct {
  /* all results in one place */
  real      volume;      /* calculated volume */
  real      reff;        /* effective radius */
  t_tgrid *tgrid;        /* full grid structure */
  real   **profile;      /* (z, R*(z)) */
  real   **xyproj;       /* projection of the occupancy in the xy plane */
  real   **xzproj;       /* projection of the occupancy in the xz plane */
} t_result;

extern gmx_bool grid_write (FILE *,t_tgrid *,char *);
extern gmx_bool grid_read  (FILE *,t_tgrid *,char *);
extern gmx_inline double DeltaV (t_tgrid *);
extern gmx_inline double DeltaA (t_tgrid *);
extern gmx_inline double gridV  (t_tgrid *);
extern gmx_inline double gridA  (t_tgrid *);
extern gmx_bool setup_tgrid (t_tgrid *, t_cavity *, rvec);
extern void setup_corners (rvec [], t_cavity *);
extern void free_tgrid (t_tgrid *);

#endif  /* _grid3D_h */
