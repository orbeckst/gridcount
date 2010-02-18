#!/usr/bin/perl -w
# $Id$

# usage
#       $0 path/to/prefix.tpr species [delta||0.1]
#


# for a windows directory, produce density grid for species S
# exclude the restrained particle (last index in P.ndx)
#
# 1. find TPR and deduce prefix P
# 2. XTC, NDX
# 3. create new NDX without sampled particle (2nd entry in NDX)
# 4. run new version of g_ri3Dc

# input: Delta, species
# output: gridxdr, meta data?

use File::Basename;

$density_d = "Density.d";
%nspec = ( NA => 14, CL => 15, SOL => 13 );
%mswitch = ( NA => "-nom", CL => "-nom", SOL => "-m" );

$TPR     = $ARGV[0];
$SPECIES = ($ARGV[1] || "UNK");
$DELTA   = ($ARGV[2] || 0.1);

die "TPR file '$TPR' missing" unless -f $TPR;
die "Unknown/missing species '$SPECIES'" unless ($nspec{$SPECIES});

$dir = dirname $TPR;
$TPR = basename ${TPR};

# prefix from trp
($pref = $TPR) =~ s(\.tpr)();

$XTC="${pref}.xtc";
$NDX="${pref}.ndx";

#echo $pref;

chdir $dir || die "failed to enter ${dir}";

die "NDX $NDX not found" unless -f $NDX;
die "XTC $XTC not found" unless -f $XTC;

mkdir ${density_d};

$tmp_ndx = "${density_d}/tmp.ndx";
$s_ndx   = "${density_d}/${SPECIES}.ndx";
$xdr     = "${density_d}/${SPECIES}xdr.dat";

# setting up the index file
system("echo -e \"keep $nspec{$SPECIES}\nq\" | make_ndx -f ${TPR} -o ${tmp_ndx}");
system("echo -e '0 & ! 2'\"\nkeep 3\nq\" | make_ndx -f ${TPR} -n ${tmp_ndx} ${NDX} -o ${s_ndx}");

# running the grid counter
$opts = $mswitch{$SPECIES};
system("g_ri3Dc  -s ${TPR} -f ${XTC} -n ${s_ndx} ${opts} -delta ${DELTA} -grid ${xdr}"); 

unlink $tmp_ndx;

exit 0;

