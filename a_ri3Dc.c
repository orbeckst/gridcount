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

/* a_ri3Dc
   -------

   Purpose: 

   analyse occupation maps (3D grids) produced by g_ri3Dc

   Typically application: count water molecules, and determine the
   true solvent accessible volume of a pore.

     V = Sum_{ijk: N(i,j,k)>N_min} deltaV

   ijk       cell (0=<i<M for XX,YY,ZZ)
   N(ijk)    occupancy of cell
   N_min     OCCUPIED_MIN
   delta_V   volume of a cell = delta_x * delta_y * delta_z

   R*: effective radius: Radius of gyration of the density in the slice at z

   Output:

   profile:  (z, R*(z))
   projection (xfarbe format) on
   xy:       title
             MX MY
	     y=0  x
	     y=1  x
	     ...
   xz:       title
             MX MZ
	     z=0  x
	     z=1  x
	     ...
*/

/* constants for RAD_BA_STEP moving average on the first RAD_BA_NR
   bins of the radial distribution functions 
   Doesnt look very promising; use from the hidden options
*/
#define RAD_BA_NR   0
#define RAD_BA_STEP 1

#define OCCUPIED_MIN 0       /* occupancy must be greater than this to
                                  be counted as occupied for volume
                                  calculations */
#define RWATER 0.14          /* canonical radius of a water molecule */

void update_tgrid (t_tgrid *);
void tgrid2cavity (t_tgrid *, t_cavity *);
static gmx_inline bool bInCircle(real,real,real,real,real);
static gmx_inline bool bInRing(real,real,real,real,real,real);

/* complete tgrid from data hold there */
void update_tgrid (t_tgrid *tg) {
  int i;
  /* top right corner */
  for(i=0;i<DIM;i++)
    tg->b[i] = tg->a[i] + tg->mx[i] * tg->Delta[i];
  return;
}

void tgrid2cavity (t_tgrid *tg, t_cavity *c) {
  rvec u,v;

  new_rvec(c->axis, 0,0,1);           /* fixed along z */
  rvec_add(tg->a,tg->b, u);
  svmul(0.5,u, c->cpoint);            /* point on axis */
  rvec_sub(tg->b,tg->a, u);
  svmul(0.5,u, v);
  c->radius = MIN(fabs(v[XX]),fabs(v[YY]));
  c->z1     = MIN(tg->a[ZZ],tg->b[ZZ]);
  c->z2     = MAX(tg->a[ZZ],tg->b[ZZ]);
  c->vol    = PI * sqr(c->radius) * (c->z2 - c->z1);
  return;
}


static gmx_inline bool bInCircle(real x,real y, real cx, real cy, real r) {
  return (x-cx)*(x-cx) + (y-cy)*(y-cy) <= r*r;
}

/* unused */
static gmx_inline bool bInRing(real x,real y, real cx, real cy, real r, real dr) {
  real u = (x-cx)*(x-cx) + (y-cy)*(y-cy); 
  return u >= r*r  &&  u < (r+dr)*(r+dr);
}

/* snap a change on the nearest grid point: 
     x_new = x_old + discrete_adjust 
   (this sign convention is more convenient in what follows)

   (Hmph; I could have written this mostly with integers ... all this
   casting/rounding is quite horrible); in this stupid hack I fudge the
   rounding
*/
real discrete_adjust (real x_old, real x_new, real delta) {
  return -(int)((x_old-x_new)/delta + (x_old > x_new ? 0.5 : -0.5)) * delta;
}

/* Change the grid according to user's radius, z1, z2 */
bool readjust_tgrid(t_tgrid *tg,t_cavity *g,int npa,t_pargs pa[]) {
  bool     bReadjusted = FALSE;
  real     dd, R_new;
  real     DeltaR = 0.5*(tg->Delta[XX]+tg->Delta[YY]); /* lacking elegance... */
  t_tgrid  old;
  int      i,j,k;
  int      base[3];    /* offsets of the indices in old.grid compared
                          to the smaller new one */

  R_new = g->radius;   /* g->radius will be overwritten later on */

  /* save the needed bits of the old grid */
  old.grid = tg->grid;            
  copy_ivec(tg->mx,    old.mx);
  copy_rvec(tg->Delta, old.Delta);
  copy_rvec(tg->a,     old.a);
  copy_rvec(tg->b,     old.b);

  /* tgrid values are the base values, the geometry is updated later */
  if(opt2parg_bSet("-z1",npa,pa)) {
    tg->a[ZZ] += discrete_adjust(tg->a[ZZ],g->z1,tg->Delta[ZZ]);
    msg("Readjusted: z1 = %g nm (supplied: %g)\n", tg->a[ZZ], g->z1);
    bReadjusted = TRUE;
  }
  if(opt2parg_bSet("-z2",npa,pa)) {
    tg->b[ZZ] += discrete_adjust(tg->b[ZZ],g->z2,tg->Delta[ZZ]);
    msg("Readjusted: z2 = %g nm (supplied: %g)\n", tg->b[ZZ], g->z2);
    bReadjusted = TRUE;
  }
  /* sanity check .. post mortem ? */
  if (tg->a[ZZ] >= tg->b[ZZ]) 
    fatal_error(0, "FAILURE: Apparently z2 =< z1, or either of them was too large "
		"for the given grid.");

  tgrid2cavity(tg,g);  /* overwrite radius from grid, finally set z1 
			  and z2 and cpoint */

  if(opt2parg_bSet("-R",npa,pa)) {
    g->radius += (dd = discrete_adjust(g->radius,R_new,DeltaR));
    msg("Readjusted: R = %g nm (supplied: %g)\n", g->radius, R_new);
    bReadjusted = TRUE;
  }

  if (!bReadjusted) return FALSE;
  /* no harm done, go back ... */

  /* (1) change grid from the current values (self consistently ;-) )
         -- this includes a new malloc for the new grid. A bit
         wasteful of memory but I cannot be bothered to use Numerical
         Recipes clever grid-reinterpretation routines with offsets
     (2) copy pieces of the old grid into the new one
  */
  if (!setup_tgrid(tg,g,old.Delta))
    fatal_error(0,"FAILURE while readjusting the grid.\n");
  
  /* Base = ((old.a - new->a)/Delta) -- beware rounding errors? */
  for(i=XX;i<DIM;i++) {
    base[i] = (int) ((tg->a[i] - old.a[i])/tg->Delta[i]);
    assert(base[i] >= 0);
    assert(base[i] + tg->mx[i] <= old.mx[i]);
  }
  for(i=0;i<tg->mx[XX];i++)
      for(j=0;j<tg->mx[YY];j++)
	  for(k=0;k<tg->mx[ZZ];k++)
	    tg->grid[i][j][k] = old.grid[i+base[XX]][j+base[YY]][k+base[ZZ]]; 

  grid3_free(old.grid);
  return TRUE;
}
/* 
   flag filename option status (from filenm.c) */
/*
#define UN_SET(fn) (fn.flag = (fn.flag & ~ffSET))
#define DO_SET(fn) (fn.flag = (fn.flag |  ffSET))
*/

/* set all options given in opts (and set to default filenames) */
/* requires set_filenms which is not in a library 

--- leave it out for the moment

void force_fnopts (int nopts, char *opts[], int nfile, t_filenm fnm[]) {
  int i,o;
  for (o=0; o<nopts; o++) {
    for(i=0; (i<nfile); i++) {
      if (strcmp(opts[o],fnm[i].opt)==0) {
	DO_SET(fnm[i]);
	break;
      }
    }
  }
  set_filenms(nfile,fnm);
  return;
}
    
*/
    
    
int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]a_rid3Dc[TT] analyses a 3D grid produced by [TT]g_ri3Dc[TT] [1]."
    "The 2D projections are written in a format suitable for the fast 2D density "
    "plotter xfarbe [2]. "
    "It also writes a parameter file `" XFARBE "' which can be used automatically "
    "by setting the environment variable `[TT]XAPPLRESDIR=.[TT]', and then running "
    "xfarbe, e.g. `[TT]xfarbe rzp.dat[TT]'\n"
    "All geometry options (radius, z1, z2) that are not set and appear as "
    "'0' are set from the grid dimensions once it is read in.\n"
    "The [TT]-dump[TT] option converts a binary grid file into ascii txt. "
    "The [TT]-plt[TT] option produces a binary density file for gOpenMol, "
    "which you might have to rename to xxx.plt. VMD [3] since version 1.8.2 "
    "can also read plt files and render them. "
    "The average density in the test cylinder is printed on std out, using the set units."
    "[PAR]"
    "Output\n"
    "------\n"
    "xy, z-averaged:        1/Lz Int_z n(x,y,z)\n"
    "xz, y-averaged:        1/Ly Int_y n(x,y,z)\n"
    "yz, x-averaged:        1/Lx Int_x n(x,y,z)\n"
    "radially averaged:  1/2pi r Int_phi n(r,phi,z)\n"
    "pore profile:       (z, R*(z), R*(z)/R)\n"
    "radial distribution function (rdf): \n"
    "               P(r) = 1/Lz 1/2pi r Int dphi dz n(r,phi,z)\n"
    "axial distribution function (zdf):\n"
    "               P(z) = 1/Lx 1/Ly Int dx dy n(x,y,z)\n"
    "P(r) and P(z) are actually averaged densities (normalisation is straightforward "
    "in order to turn them into 'real' probability distributions).\n"
    "The local density axial distribution function -lzdf averages over all  "
    "occupied grid cells per z-slice and divides by an effective area which is "
    "determined from the 'radius of gyration' of the density:\n" 
    "               R*(z)^2 = 2 Int_0^2pi dphi Int_0^R dr r^3 n(r,phi,z)/"
    "Int drho dphi r n(r,phi,z)\n"
    "This radius is a good approximation to the "
    "pore profile (eg compared to HOLE). It describes the solvent accessible surface (SAS) "
    "of the pore because the density is based on the centres of the molecules. "
    "The density itself is in the lzdf.xvg file."
    "[PAR]For diagnostic purposes one can also plot the radial distributions of "
    "the unoccupied cells (holes in the grid) in order to find suitable grid "
    "spacings.\n"
    "[PAR]",
    "---------\n"
    "[1] Oliver Beckstein, http://sansom.biop.ox.ac.uk/oliver/software/#gridcount\n"
    "[2] A. Preusser, http://www.fhi-berlin.mpg.de/~grz/pub/xfarbe.html\n",
    "[3] Visual Molecular Dynamics http://www.ks.uiuc.edu/Research/vmd/\n"
    "[PAR]Known limitations:"
  };

  static char *bugs[] = {
    "The radial bin width DeltaR is fixed to (Delta[XX]+Delta[YY])/2). "
    "In any case one should never have different bin widths in X and Y.",
    "There are still a few hidden options of questionable usefulness. "
    "Resampling (=changing Delta) is not implemented yet.",
    "gOpenMol plt binary file comes out with wrong suffix",
    "note: -minocc also influences -lzdf"
  };

  static t_cavity geometry = {   /* describes the cylinder */
    {0, 0, 1},            /* axis -- cannot be changed */
    {0, 0, 0},            /* cpoint -- set from grid */
    0,                    /* radius */
    0, 0,                 /* z1 < z2 */
    0                     /* volume - calculate later */
  };
  static real min_occ = OCCUPIED_MIN;    /* count cell as occupied 
					    if p > min_occ */
  static bool bMirror = TRUE;            /* pretty picture of P(r,z) */
  static bool bDoHoles = FALSE;          /* use instead of giving fns
                                            explixitly */
  static char *DensUnitStr[] = 
         { NULL, "unity", "SPC", "molar", "Angstrom", NULL };
  static rvec Delta = {0.02,0.02,0.02};  /* resolution in nm for
                                            coarse graining */
  static int  rad_ba_nr =   RAD_BA_NR;   /* block average of radial
                                            bins */
  static int  rad_ba_step = RAD_BA_STEP; 
  static char buf[HEADER_MAX];           /* additional text for graphs */
  static char *header = buf;

  /* maximum xfarbe level (in nm^-3) (this is motivated by the
   observation that densities for water in excess of 1.5 n_bulk rarely
   occur) */
  static real maxDensity = 48.48405;
  static real rsolvent=RWATER; /* radius of a solvent molecule or ion */


  t_pargs pa[] = {
    { "-R",      FALSE, etREAL, {&(geometry.radius)},
      "Maximum radius that should be contained in the grid"},
    { "-z1",     FALSE, etREAL, {&(geometry.z1)},
      "Center grid between z1, and ..."},
    { "-z2",     FALSE, etREAL, {&(geometry.z2)},
      "z2 (these boundaries are kept fixed!)"},
    { "-delta",  FALSE, etRVEC, {&Delta},
      "HIDDENSpatial resolution in X, Y, and Z for resampling (in nm)"},
    { "-minocc",   FALSE, etREAL, {&min_occ},
      "HIDDENThe occupancy of a  cell must be larger than this number so that it is "
      "counted as occupied when calculating the volume, effective radius "
      "and local density axial distribution -lzdf. This is given in the chosen "
      "units (see -unit)."},
    { "-rsolvent", FALSE, etREAL, {&rsolvent},
      "HIDDENradius of a solvent moelecule or ion; used to incorporate "
      "the excluded volume in n(z) and is shown in the 3rd column in lzdf"},
    { "-subtitle", FALSE, etSTR, {&header},
      "Some text to add to the output graphs"},
    { "-mirror", FALSE, etBOOL, {&bMirror},
      "mirror the radial projection P(r,z) to create the impression of a "
      "full view of the pore"},
    { "-unit", FALSE, etENUM, {DensUnitStr},
      "divide number density (in nm^-3) by 1, the density of SPC "
      "water at 300K and 1 bar, in mol/l or in "
      "Angstrom^-3. Allowed values"},
    { "-xfarbe-maxlevel", FALSE, etREAL, {&maxDensity},
      "xfarbe will plot 15 equally space level up to this density (the unit must "
      "be the same as for the -unit option!) Default is 1.5 SPC bulk." },
    { "-holes",  FALSE, etBOOL,  {&bDoHoles},
      "HIDDENcalculate all hole distribution functions, using "
      "default filenames" },
    { "-radnr",   FALSE, etINT, {&rad_ba_nr},
      "HIDDENblock-average the first rad_ba_nr radial bins..." },        
    { "-radstep", FALSE, etINT, {&rad_ba_step},
      "HIDDEN...in blocks of rad_ba_step"},
  };

  t_filenm fnm[] = {
    { efDAT, "-grid",    "gridxdr", ffREAD },
    { efXVG, "-profile", "profile", ffWRITE },
    { efDAT, "-xyp",  "xyp",  ffWRITE },
    { efDAT, "-xzp",  "xzp",  ffOPTWR },
    { efDAT, "-yzp",  "yzp",  ffOPTWR },
    { efDAT, "-rzp",  "rzp",  ffWRITE },
    { efXVG, "-rdf",  "rdf",  ffWRITE },
    { efXVG, "-zdf",  "zdf",  ffWRITE },
    { efXVG, "-lzdf", "lzdf", ffWRITE },
    { efDAT, "-hxyp", "hxyp", ffOPTWR },
    { efDAT, "-hrzp", "hrzp", ffOPTWR },
    { efXVG, "-hrdf", "hrdf", ffOPTWR },
    { efDAT, "-dump", "gridasc",  ffOPTWR },
    { efDAT, "-plt",  "plt",  ffOPTWR },       /* gOpenMol density, see plt.c */
  };

  FILE       *fGrid;         /* 3D grid with occupation numbers */
  FILE       *fGrid2;        /* resampled */
  FILE       *fOut;          /* reusable fp */
  FILE       *fDump;         /* dump the grid as ascii -- use sparringly */
  t_tgrid  tgrid;            /* all information about the grid */
  real       **profile=NULL;   /* (z, R*(z)) */
  real       **xyp=NULL;       /* Sum_z P(x,y,z) */
  real       **xzp=NULL;       /* Sum_y P(x,y,z) */
  real       **yzp=NULL;       /* Sum_x P(x,y,z) */
  real       **rzp=NULL;       /* Sum_phi P(r,phi,z) */
  real       **hxyp=NULL;   /* distribution of unoccupied cells */
  real       **hrzp=NULL;   /* distribution of unoccupied cells */
  double     sumP,sumH;     /* Sum_xyz P(x,y,z); sumH is over
                               unoccupied cells */
  real       **rdf=NULL;          /* n(r) */
  real       **zdf=NULL;          /* n(z) */  		
  real       **lzdf=NULL;         /* n(z), but from occupied cells only */
  real       **hrdf=NULL;

  enum eDensUnit nunit = eduUNITY;       /* how to measure the density */

  int i,j,k;
  int nocc;         /* number of occupied (> min_occ) cells */
  real Rgyr2;       /* 'radius of gyration' squared for occupied cells */
  real Agyr;        /* area Agyr = pi * Rgyr^2 */  
  /* used to approximate the excluded volume */
  double dV;        /* volume of a cell */
  real height;      /* z2 - z1 */
  real volume = 0;
  real reff;        /* effective radius */

  real DeltaR;      /* radial resolution */
  real ci, cj;      /* centered grid coordinates */
  int irad,iradNR;  /* index for the radius */
  real *rweight;    /* weights for radial bins */
  real average = 0; /* average density in count cylinder */
  int ivol = 0;     /* number of grid cells in the cylinder ("volume") */

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

  assert(fmod((real)rad_ba_nr, (real)rad_ba_step) < GMX_REAL_EPS);

  if (bDoHoles) {
    char *holeopts[] = {"-hrdf", "-hrzp", "-hxyp"};
#define NHOPTS asize(holeopts)
    msg("Sorry, this option does not work at the moment.\n");
    /*
    force_fnopts(NHOPTS,holeopts,NFILE,fnm);
    */
  }

  /* select the unit (DensUnit[nunit]) for the density plots */
  switch (DensUnitStr[0][0]) {
  case 'S': /*SPC water at 300K, 1 bar */
    nunit = eduSPC;
    break;
  case 'm': /* in mol/l */
    nunit = eduMOLAR;
    break;
  case 'A': /* in Angstrom^-3 */
    nunit = eduANG;
    break;
  case 'u': /* unity */
  default:
    nunit = eduUNITY;
  } /*end switch */

  min_occ *= DensUnit[nunit]; /* compare filling of cell in chosen unit */

  if(!opt2parg_bSet("-xfarbe-maxlevel",NPA,pa)) {
    maxDensity /= DensUnit[nunit];
    dmsg("Converted the default xfarbe-maxlevel: %f %s\n",
	 maxDensity, EDENSUNITTYPE(nunit));
  }

  /* open input file */
  msg("Reading grid 3D file...\n");
  fGrid    = ffopen (opt2fn("-grid", NFILE, fnm), "r");
  if (!grid_read(fGrid,&tgrid,header)) 
    fatal_error(0,"Error reading the 3D grid---no point in continuing!\n");
  update_tgrid(&tgrid);
  fclose(fGrid);
  msg(".. done!\n");

  /* 
     (1) for each user-changed parameter radius, z1, z2 change the grid
     (2) set 'geometry' from the grid box
  */
  readjust_tgrid(&tgrid,&geometry,NPA,pa);
  tgrid2cavity(&tgrid,&geometry);

  DeltaR = 0.5*(tgrid.Delta[XX]+tgrid.Delta[YY]);
  if (tgrid.Delta[XX] != tgrid.Delta[YY]) 
    dmsg("Warning: Delta[XX]=%.4f != Delta[YY]=%.4f, setting DeltaR=%.4f\n",
	 tgrid.Delta[XX], tgrid.Delta[YY], DeltaR);

  /* setup proj */
  if (!(xyp=grid2_alloc(tgrid.mx[XX],tgrid.mx[YY])) ||
      !(xzp=grid2_alloc(tgrid.mx[XX],tgrid.mx[ZZ])) ||
      !(yzp=grid2_alloc(tgrid.mx[YY],tgrid.mx[ZZ])))
    fatal_error (-1,"FAILED: allocating memory for the projection maps\n");

  /* setup axial distribution (zdf) and profile  
     We are wasting memory because the x-column is identical for all of them
     but its conceptually simpler.
   */
  if (! (zdf=grid2_alloc(tgrid.mx[ZZ],2)) ||
      ! (lzdf=grid2_alloc(tgrid.mx[ZZ],2)) ||
      ! (profile=grid2_alloc(tgrid.mx[ZZ],2)))
    fatal_error (-1,"FAILED: allocating memory for the axial "
		 "distribution functions or the profile\n");

  /* setup radial distributions (rdf and rzprojection)   */
  iradNR = (int)floor(geometry.radius/DeltaR);
  snew(rweight,iradNR);
  if (!rweight || !(rdf=grid2_alloc(iradNR,2)))
    fatal_error (-1,"FAILED: allocating memory for the radial "
		 "distribution function\n");
	
  if (!(rzp=grid2_alloc(iradNR,tgrid.mx[ZZ])))
    fatal_error (-1,"FAILED: allocating memory for the rz projection map\n");

  /* hole distributions */
  if (!(hrdf=grid2_alloc(iradNR,2)) || 
      !(hrzp=grid2_alloc(iradNR,tgrid.mx[ZZ])) ||
      !(hxyp=grid2_alloc(tgrid.mx[XX],tgrid.mx[YY])))
    fatal_error(-1,"FAILED: allocating memory for the hole distributions\n");

  /* for the volume calculation only look at cells within the circle
     of the original radius so that a comparison between this radius
     and the effective radius is meaningful */

  /* the grid output is a number density. Unit nm^-3   */
  /* sum over the whole grid (for the integral, multiply by dV) */
  sumP = 0; sumH = 0;
  for(k=0;k<tgrid.mx[ZZ];k++) 
    for(j=0;j<tgrid.mx[YY];j++) 
      for(i=0;i<tgrid.mx[XX];i++) {
	sumP += (double) tgrid.grid[i][j][k];
	if (tgrid.grid[i][j][k] <= min_occ)  sumH++;
      }
  dV = DeltaV(&tgrid);
  sumP *= dV;     /* dont forget the dV ... */
  sumH *= dV;
  dmsg("Integrals over the whole grid:\n"
       "  sumP = <<N>>     = Sum_ijk dV*n[ijk] = %.2f\n"
       "  sumH = <<holes>> = Sum_ijk dV*h[ijk] = %.2f\n",
       sumP,sumH);
    
  /* NB:     Integral_L dx f(x) =       Sum_i DELTAx f_i
         1/L Integral_L dx f(x) = 1/n_x Sum_i        f_i
  */
  for(k=0;k<tgrid.mx[ZZ];k++) {
    nocc=0;      /* number of occupied cells in this xy layer */
    Rgyr2=0;      /* radius of gyration for occupied cells */
    for(j=0;j<tgrid.mx[YY];j++) {
      for(i=0;i<tgrid.mx[XX];i++) {
	/* projections on planes */
	xyp[i][j] += tgrid.grid[i][j][k] / 
	             ((real)tgrid.mx[ZZ] * DensUnit[nunit]);
	xzp[i][k] += tgrid.grid[i][j][k] / 
	             ((real)tgrid.mx[YY] * DensUnit[nunit]);
	yzp[j][k] += tgrid.grid[i][j][k] / 
	             ((real)tgrid.mx[XX] * DensUnit[nunit]);
	if (tgrid.grid[i][j][k] <= min_occ) 
	  hxyp[i][j] += 1.0/((real)tgrid.mx[ZZ] * DensUnit[eduUNITY]);

	/* translate to coordinates centered on the axis and
	   pointing to the center of the cell 
           These are also needed for radial binning:
        */
	ci = ((real)i+0.5) - (real)tgrid.mx[XX]/2.0;
	cj = ((real)j+0.5) - (real)tgrid.mx[YY]/2.0;

	if (tgrid.grid[i][j][k] > min_occ && 
	    bInCircle((real)i, (real)j, 
		      (real)tgrid.mx[XX]/2.0, (real)tgrid.mx[YY]/2.0, 
		      geometry.radius/DeltaR))  {
	  /* multiply by dV later */
	  volume++;        /* total number of occupied cells */
          nocc++;          /* per k-slice; not used, just for debugging or future use */

          /* local density axial distribution function
             the area is estimated as the 'radius of gyration' of the occupied 
             cells 

             Rgyr = sqrt( Sum_i,occ n(i,z)*r(i,z)^2 / Sum_i,occ n(i,z)  )
             where r(i,z) is the distance from the center of the pore
             Agyr = pi * Rgyr^2
           
             units of DeltaR for radius
          */
        
          /* Rgyr2 could equally be calculated from ALL cells; empty
             grid cells don't have any weight in the sum)
          */
          lzdf[k][1] += tgrid.grid[i][j][k]; 
          Rgyr2 += tgrid.grid[i][j][k] * (ci*ci + cj*cj);
	}


	/* ------------------------------------------------------------
           radial distributions 
           ------------------------------------------------------------ 
        */

	/* radius in units DeltaR (really:
	   sqrt((i*DX)^2+(j*DY)^2)/DeltaR but with the fixed setting
	   of DeltaR this cancels out) 
        */
	irad = (int)floor(sqrt(ci*ci + cj*cj));
	
	if (irad < iradNR) {
          /* grand average in the cylinder 
             (divide by cylinder volume later (which is NOT volume---'volume' 
             only contains the occupied volume). ivol counts the number of
             grid cells in the cylinder and is the discrete volume.
          */
          average += tgrid.grid[i][j][k];
          ivol++;
	  /* count grid cells in this radial bin. Because we are
             looping over z as well we need to average over all z
             slices 
          */
	  rweight[irad] += 1.0/(real)tgrid.mx[ZZ];     
	  if (tgrid.grid[i][j][k] > min_occ) {
	    /* multiply by norm later */
	    rdf[irad][1] += tgrid.grid[i][j][k];
	    rzp[irad][k] += tgrid.grid[i][j][k];
	  } else {
	    /* distribution of unoccupied cells */
	    hrdf[irad][1]++;    
	    hrzp[irad][k]++; 
	  }
	}
	/* axial distribution function zdf (normalise later) Note:
           this integrates square disks out and collapses them onto
           the z-axis, NOT circles. (should be fixed) The 'radius' R
           determines the maximum grid dimensions ie length of the
           'square disk'.

           Could be fixed by moving it into the 'profile' loop but I
           leave it for consistency with previous calculations for the
           time being. See lzdf for alternative.

        */
	zdf[k][1] += tgrid.grid[i][j][k];

      }
    }
    /* zdf --  z at the center of each cell --> +0.5 
       & normalise by the area (square disks -> Lx * Ly)
    */
    zdf[k][0]  = tgrid.a[ZZ] + (k+0.5)*tgrid.Delta[ZZ];
    zdf[k][1] /= (real)tgrid.mx[XX]*(real)tgrid.mx[YY] * DensUnit[nunit];
    lzdf[k][0] = profile[k][0] = zdf[k][0];     /* z coordinate */

    /* normalisation; factor 2 from comparison with result for const
       density; see Labbook V, p11 */  
    Rgyr2 /= 0.5*lzdf[k][1];  

    /* pore profile (z, R*(z))  */
    profile[k][1] = sqrt(Rgyr2)*DeltaR;         /* Rgyr, nm  */

    Agyr = PI * Rgyr2;                          /* effective area */

    /* Axial density profile (radially averaged) using local density
       and an approximation to the pore radius as the radius of
       gyration of the density, Rgyr
       
       n(z) = Sum n(x,y,z) * dV / (A*dz)
       lzdf[k] * dV: avg number of particles in slice.  
       Divide by total volume of the slice with the effective area
       A=pi*Rgyr2

       dV/(Delta[ZZ]*sqr(DeltaR) == 1 but we do it explicitly because
       it might break when code related to DeltaR is changed

       lzdf[k][1] /= Agyr * DensUnit[nunit];
     */
    lzdf[k][1] = lzdf[k][1]*dV/(tgrid.Delta[ZZ]*Agyr*DeltaR*DeltaR * DensUnit[nunit]);

#ifdef DEBUG
      printf("[k=%d] dV=%g dz=%g dr=%g nocc=%d Rgyr=%g Rgyr/R=%g Agyr=%g nm^2 Agyr=%g\n",
             k,
             dV,tgrid.Delta[ZZ],DeltaR,nocc,
             profile[k][1], lzdf[k][2]/geometry.radius, Agyr*DeltaR*DeltaR, Agyr);
#endif
  } /* end loop over k-slices (ie z) */

  /* radial distributions 
     (1) normalisation (important!!) 
     (2) block average the innermost bins: 0+1 and 2+3 to smooth 
         out artifactual  wiggles
      --- (2) does not really help but I leave it in for the time being
  */

  /* (1)
     norm for rdf:  dz * dx*dy/ring_area(r) 
     norm for rzp:       dx*dy/ring_area(r) 
     unit for ring area == dx * dy !
   */
  for(irad=0; irad<iradNR; irad++) {
    rdf[irad][0] = hrdf[irad][0] = irad * DeltaR;
    assert(rweight[irad] > 0);
    rdf[irad][1]     /= rweight[irad]*(real)tgrid.mx[ZZ] * DensUnit[nunit];
    hrdf[irad][1]    /= rweight[irad]*(real)tgrid.mx[ZZ] * DensUnit[eduUNITY];
    for(k=0; k<tgrid.mx[ZZ]; k++) {
      rzp[irad][k]   /= rweight[irad] * DensUnit[nunit];
      hrzp[irad][k]  /= rweight[irad] * DensUnit[eduUNITY];
    }
  }
  /* (2) do the block average 
   ..and Im sure one can do it a lot more elegantly than this m(a|e)ss 
   of for loops
  */ 
  for(irad=0;irad < rad_ba_nr;irad += rad_ba_step) {
    for(i=1;i<rad_ba_step;i++) {
      /* accumulate in first element of the average */
      rdf[irad][1]  +=  rdf[irad+i][1];
      hrdf[irad][1] += hrdf[irad+i][1];
      for(k=0; k<tgrid.mx[ZZ]; k++) {
	 rzp[irad][k]  +=  rzp[irad+1][k];   
         hrzp[irad][k] += hrzp[irad+1][k];      
      }
    }
    /* average */
    rdf[irad][1]  /= (real) rad_ba_step;
    hrdf[irad][1] /= (real) rad_ba_step;
    for(k=0; k<tgrid.mx[ZZ]; k++) {
      rzp[irad][k]  /= (real) rad_ba_step;
      hrzp[irad][k] /= (real) rad_ba_step;
    }
    /* fill the remaining bins after the first one with the average */
    for(i=1;i<rad_ba_step;i++) {
      rdf[irad+1][1]  = rdf[irad][1];
      hrdf[irad+1][1] = hrdf[irad][1];
      for(k=0; k<tgrid.mx[ZZ]; k++) {
        rzp[irad+1][k]  =  rzp[irad][k];   
	hrzp[irad+1][k] = hrzp[irad][k];      
      }
    }
  }

  volume *= dV;                      /* volume of occupied cells */
  height = tgrid.b[ZZ] - tgrid.a[ZZ];
  reff = sqrt(volume/(PI*height));
  average /= ivol * DensUnit[nunit]; /* discrete cylinder volume = ivol * dV */

  msg("\n# ---> Result <-----\n"
      "# Densities have the unit %s\n"
      "  V=%g nm³  R*=%g nm sumP=%g sumH=%g dV=%f\n",
      EDENSUNITTYPE(nunit),volume,reff,sumP,sumH,dV);
  printf("\nDATA_RESULT:  V %f   R* %f   sumP %f sumH %f  dV %f # %s\n",
	 volume,reff,sumP,sumH,dV,header);
  printf("DATA_RESULT:  average %f  iVolume %f\n", average, (real)ivol*dV);

  /* 
     write output 
  */

  /* pore radius profile, based on how many particles fill the probed cylinder 
     Format:
     z/nm   R* /nm   R* /R 
                     (should approach 1 for bulk)
   */ 
  fOut = xmgropen (opt2fn("-profile",NFILE,fnm),
		       "Effective radius R*(z) pore profile",
		       header, "z [nm]", "R* [nm]");
  for(k=0;k<tgrid.mx[ZZ];k++) {
    fprintf (fOut,"%.6f   %.6f  %.6f\n",
	     profile[k][0],profile[k][1],profile[k][1]/geometry.radius);
  };
  fclose(fOut);


  /* projections */
  /* see 
     http://www.fhi-berlin.mpg.de/gnz/pub/xfarbe/xfarbe-2.6/html/#Format%20of%20Data%20Files 
  */
  
  fOut = xf_open(opt2fn("-xyp",NFILE,fnm),header,tgrid.mx[XX],tgrid.mx[YY]);
  for(j=0;j<tgrid.mx[YY];j++) {
    for(i=0;i<tgrid.mx[XX];i++) 
      fprintf(fOut,"%f ",xyp[i][j]);
    fprintf(fOut,"\n");
  }
  xf_append_levels(fOut,XFARBE_NCOLS,maxDensity);
  xf_append_axes_annotation(fOut,tgrid.a[XX],tgrid.b[XX],
			    tgrid.a[YY],tgrid.b[YY]);
  fclose(fOut);

  if (opt2bSet("-hxyp",NFILE,fnm)) {
    fOut = xf_open(opt2fn("-hxyp",NFILE,fnm),header,tgrid.mx[XX],tgrid.mx[YY]);
    for(j=0;j<tgrid.mx[YY];j++) {
      for(i=0;i<tgrid.mx[XX];i++) 
	fprintf(fOut,"%f ",hxyp[i][j]);
      fprintf(fOut,"\n");
    }
    xf_append_levels(fOut,XFARBE_NCOLS,1);
    xf_append_axes_annotation(fOut,tgrid.a[XX],tgrid.b[XX],
			      tgrid.a[YY],tgrid.b[YY]);
    fclose(fOut);
  }

  if (opt2bSet("-xzp",NFILE,fnm)) {      
    fOut = xf_open(opt2fn("-xzp",NFILE,fnm),header,tgrid.mx[XX],tgrid.mx[ZZ]);
    for(k=0;k<tgrid.mx[ZZ];k++) {
      for(i=0;i<tgrid.mx[XX];i++)
	fprintf(fOut,"%f ",xzp[i][k]);
      fprintf(fOut,"\n");
    }
    xf_append_levels(fOut,XFARBE_NCOLS,maxDensity);
    xf_append_axes_annotation(fOut,tgrid.a[XX],tgrid.b[XX],
			    tgrid.a[ZZ],tgrid.b[ZZ]);
    fclose(fOut);
  }

  if (opt2bSet("-yzp",NFILE,fnm)) {      
    fOut = xf_open(opt2fn("-yzp",NFILE,fnm),header,tgrid.mx[XX],tgrid.mx[ZZ]);
    for(k=0;k<tgrid.mx[ZZ];k++) {
      for(j=0;j<tgrid.mx[YY];j++)
	fprintf(fOut,"%f ",yzp[j][k]);
      fprintf(fOut,"\n");
    }
    xf_append_levels(fOut,XFARBE_NCOLS,maxDensity);
    xf_append_axes_annotation(fOut,tgrid.a[YY],tgrid.b[YY],
			    tgrid.a[ZZ],tgrid.b[ZZ]);
    fclose(fOut);
  }


  fOut = xf_open(opt2fn("-rzp",NFILE,fnm),header,
                   (bMirror ? 2*iradNR - 1 : iradNR), tgrid.mx[ZZ]);
  for(k=0;k<tgrid.mx[ZZ];k++) {
    for(irad=(bMirror ? -iradNR + 1 : 0); irad<iradNR; irad++) 
      fprintf(fOut,"%f ",rzp[abs(irad)][k]);
    fprintf(fOut,"\n");
  }
  xf_append_levels(fOut,XFARBE_NCOLS,maxDensity);
  xf_append_axes_annotation(fOut, (bMirror ? 
	      -geometry.radius + DeltaR  :  geometry.radius),
	      geometry.radius,    tgrid.a[ZZ],tgrid.b[ZZ]);
  fclose(fOut);

  if (opt2bSet("-hrzp",NFILE,fnm)) {
    fOut = xf_open(opt2fn("-hrzp",NFILE,fnm),header,
                     (bMirror ? 2*iradNR - 1 : iradNR) ,tgrid.mx[ZZ]);
    for(k=0;k<tgrid.mx[ZZ];k++) {
      for(irad=(bMirror ? -iradNR + 1 : 0); irad<iradNR; irad++) 
	fprintf(fOut,"%f ",hrzp[abs(irad)][k]);
      fprintf(fOut,"\n");
    }
    xf_append_levels(fOut,XFARBE_NCOLS,1);
    xf_append_axes_annotation(fOut, (bMirror ? 
	      -geometry.radius + DeltaR  :  geometry.radius),
	      geometry.radius,    tgrid.a[ZZ],tgrid.b[ZZ]);
    fclose(fOut);
  }
  

	  
  /* 
     distribution functions 
  */

  /* (1) axial distribution P(z) */
  snprintf(s_tmp,STRLEN,"n(z) [%s], PMF G(z)/kT",EDENSUNITTYPE(nunit));
  fOut = xmgropen (opt2fn("-zdf",NFILE,fnm),
		       "Axial distribution function n(z)",
		       header, "z [nm]", s_tmp);
  for(k=0;k<tgrid.mx[ZZ];k++) {
    fprintf (fOut,"%.6f   %.6f   %.6f\n",
	     zdf[k][0],zdf[k][1],-log(zdf[k][1]));
  };
  fclose(fOut);

  snprintf(s_tmp,STRLEN,"n\\slocal\\N(z) [%s]",EDENSUNITTYPE(nunit));
  fOut = xmgropen (opt2fn("-lzdf",NFILE,fnm),
		       "Local density axial distribution function n\\slocal\\N(z)",
		       header, "z [nm]", s_tmp);
  for(k=0;k<tgrid.mx[ZZ];k++) {
    fprintf (fOut,"%.6f   %.6f  %.6f  %.6f\n",
	     lzdf[k][0],lzdf[k][1],
             lzdf[k][1]/(1+sqr(rsolvent/profile[k][1])),
	     -log(lzdf[k][1]) ); 
             /* z  n(z) n'(z) G(z)/kT
              where n'(z) = N(z)/pi*(Rgyr+r)^2 with r: radius of solvent */
  };
  fclose(fOut);


  /* (2) radial distribution function */
  if (bDebugMode()) {
    /* print radial weigths */
    fOut = fopen ("rweights.dat", "w");
    for(irad=0; irad<iradNR; irad++) 
      fprintf(fOut,"%f %d %f\n", 
	      rdf[irad][0], (int)rweight[irad], 
	      2*PI*DeltaR*DeltaR*(irad+0.5)/DeltaA(&tgrid));
    fclose(fOut);
  }

  snprintf(s_tmp,STRLEN,"n(r) [%s]",EDENSUNITTYPE(nunit));
  fOut = xmgropen (opt2fn("-rdf",NFILE,fnm),
		       "Radial distribution function P(r)",
		       header, "r [nm]", s_tmp);
  for(irad=0;irad<iradNR;irad++) {
    fprintf (fOut,"%.6f   %.6f\n",
	     rdf[irad][0],rdf[irad][1]);
  };
  fclose(fOut);

  if (opt2bSet("-hrdf",NFILE,fnm)) {
    snprintf(s_tmp,STRLEN,"h(r) [%s]",EDENSUNITTYPE(eduUNITY));
    fOut = xmgropen (opt2fn("-hrdf",NFILE,fnm),
		     "Radial distribution function of unoccupied cells",
		     header, "r [nm]", "h(r)");
    for(irad=0;irad<iradNR;irad++) {
      fprintf (fOut,"%.6f   %.6f\n",
	       hrdf[irad][0],hrdf[irad][1]);
    };
    fclose(fOut);
  }

  /* dump grid as ascii */
  if (opt2bSet("-dump",NFILE,fnm)) {
    msg("Dumping the 3D grid into ascii file %s..."
        "[Unit is %s]",
	opt2fn("-dump",NFILE,fnm), EDENSUNITTYPE(nunit));
    density_write_ascii(opt2fn("-dump",NFILE,fnm),&tgrid,DensUnit[nunit]);
    msg(".. done!\n");
  }

  /* write gOpenMol density plt file (binary) */
  if (opt2bSet("-plt",NFILE,fnm)) {
    msg("Dumping the 3D grid into gOpenMol density file %s... (rename to xxx.plt)\n"   
        "[Unit is %s]",
	opt2fn("-plt",NFILE,fnm), EDENSUNITTYPE(nunit));
    density_write_plt(opt2fn("-plt",NFILE,fnm),&tgrid,DensUnit[nunit]);
    msg("... done!\n");
  }


  /* write parameter file for xfarbe 
     see https://indigo1.biop.ox.ac.uk/xfarbe/html/#xfarbe.file
  */
  xf_write_XFarbe(XFARBE_NCOLS);


  /* clean up a bit */

  free_tgrid(&tgrid);
  sfree(rdf);
  sfree(zdf);
  sfree(lzdf);
  sfree(hrdf);
  grid2_free(xzp);
  grid2_free(xyp);
  grid2_free(hrzp);
  grid2_free(hxyp);
  sfree(profile);

  fprintf(stderr,"\n");
  thanx(stdout);

  return 0;
}

