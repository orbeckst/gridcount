/* 
   $Id$

   writing density in plt format gOpenMol
   http://www.csc.fi/gopenmol/developers/plt_format.phtml

   Copyright (C) 2004 Oliver Beckstein <oliver@biop.ox.ac.uk>
   This program is made available under the terms of the GNU Public License. 
   See the file LICENSE or http://www.gnu.org/copyleft/gpl.html

*/

#ifndef _plt_h
#define _plt_h

static char *SRCID_plt_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "grid3D.h"

extern void density_write_plt   (char *, t_tgrid *, real);
extern void density_write_ascii (char *, t_tgrid *, real);

#endif /* _plt_h */
