/* $Id$

   Copyright (C) 2003, 2004 Oliver Beckstein <oliver@biop.ox.ac.uk>
   This program is made available under the terms of the GNU Public License. 
   See the file LICENSE or http://www.gnu.org/copyleft/gpl.html
   
*/

#ifndef _count_h
#define _count_h

static char *SRCID_count_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include "typedefs.h"
#include "names.h"
#include "vec.h"
#include "princ.h"
#include "rmpbc.h"
#include "readinp.h"
#include "utilgmx.h"

/* types of index: atoms or molecules */
enum ndxtype    {etxATOM, etxMOL, etxNR};  
enum eDensUnit {eduUNITY, eduSPC, eduMOLAR, eduANG, eduVOXELPROB, eduNR};
/* units to mesure a number density n (given in nm^-3) in:
   n' = n/DensUnit

   UNITY:
    1 nm^-3 = 1 x 1 nm^-3

   SPC: 
    rho0 = 0.9669;      # bulk SPC water at 300K (Gromacs 3, ffgmx)
                        # new: from max of P(N), labbook II, p159:
                        #            0.966894±0.000462 g cm^-3
                        #      from N/<V>, II, p158
                        #            0.966458±0.00458  g cm^-3
    # volume of one water molecule at T=300K, P=1bar (Labbook II p159)
    # in nm^3 
    v_water   = 0.030938 nm^3
    1/v_water = 32.3227  nm^-3

   MOLAR:
    1 nm^-3 = 1/N_Avogadro * (10^-8 dm)^-3

   ANG: in Angstroem^-3
    1 nm^-3 = 1 (10 A)^-3 

   VOXELPROB:
    probability for this voxel to be occupied (a value between 0 and 1):
    Simply multiply the density by the cell volume.

    n' = n / (1/(deltaX*deltaY*deltaZ))

    (Note that this is NOT a probability yet but just the average number
    per voxel. For a probability, normalise by the total sum.)
 */
extern char *endxtype_names[etxNR+1];
extern char *eDensUnit_names[eduNR+1];
extern real DensUnit[eduNR];


#define ENDXTYPE(e)      ENUM_NAME(e,etxNR,endxtype_names)
#define EDENSUNITTYPE(e) ENUM_NAME(e,eduNR,eDensUnit_names)


typedef struct {
  rvec    axis;       /* vector parallel to pore axis */
  rvec    cpoint;     /* any point on the axis */
  real    radius;     /* max radius of the cylindrical pore */
  real    z1, z2;     /* confined between z1 and z2 */
  real    vol;        /* volume V=r^2 pi (z2-z1) */
} t_cavity;


typedef struct {      /* times in the simulation */
  int   tot_frames;   /* number of frames read */
  real  t_tot;        /* simulation time (ps) */
  real  delta_t;      /* time step (ps) */
  real  tau;          /* threshold for tracking a particle: 
			 t_cav > tau */
} t_simtime;


extern gmx_bool autoset_cavity (t_cavity *,matrix,int,t_pargs []);
extern gmx_bool bInCavity (const rvec x, const t_cavity *cavity);
extern real npbc2com (t_topology *, atom_id *, matrix, rvec *, 
		      atom_id, rvec);
extern int mols_from_index (atom_id *index, int gnx, t_block *mols, 
		     atom_id *molndx, int max_mol);
extern void print_ldist (const rvec x, const t_cavity *cavity, const atom_id idx, const real mass);
extern void fwrite_index (FILE *, atom_id *, int, t_topology *, char *);

#endif   /* _count_h */
