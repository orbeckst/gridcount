/* $Id$
  
  Functions for the gridcounting programs, including some utility functions
 
  Copyright (C) 2004 Oliver Beckstein <oliver@biop.ox.ac.uk>
  This program is made available under the terms of the GNU Public License. 
  See the file LICENSE or http://www.gnu.org/copyleft/gpl.html

*/

#ifndef _libgridcount_h
#define _libgridcount_h

static char *SRCID_gridcount_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xdr_grid.h" /* read/write of portable xdr file with 3D grid data (native format)*/
#include "grid3D.h"   /* higher level 3D grid routines */
#include "count.h"    /* counting particles in defined volumes */
#include "xf.h"       /* xfarbe data and parameter file */
#include "plt.h"      /* gOpenMol 3D plot format */
#include "utilgmx.h"  /* utility functions (array allocation, debug messages) */

#endif /* _libgridcount_h */
