# GNU Makefile to compile the grid counter outside the Gromacs source tree
#
#   Copyright (C) 2003-2010 Oliver Beckstein <orbeckst@gmailcom>
#   This program is made available under the terms of the GNU Public License. 
#   See the file COPYING or http://www.gnu.org/copyleft/gpl.html
#
# Edit BIN_DIR and the GMX_* variables to reflect the locations in
# your setup. Run 'make help' for, um, help.
#
# Makefile switches:
# Switch on on commandline by setting them, eg 
#        make DEBUG=1  g_ri3Dc 
#
#  DEBUG        compile with debugging, -g and no optimisations

# Set GMX_TOP_DIR to your Gromacs installation and check the
# INCLUDE and LIB dir.
GMX_TOP_DIR     := $(HOME)/Library/Gromacs/version/4.0.2

# EXEC depends on your machine/OS
GMX_EXEC_PREFIX := $(GMX_TOP_DIR)/`config.guess`
GMX_LIB_DIR     := $(GMX_EXEC_PREFIX)/lib#
GMX_INCLUDE_DIR := $(GMX_TOP_DIR)/include/gromacs#

# install binaries into
BIN_DIR := $(GMX_EXEC_PREFIX)/bin

# This is only necessary for the creation of etags and can be safely ignored
# (for compilation it is not important).
GMX_SOURCE_DIR  := 
#
##########################################################################


##########################################################################
# Probably no need to touch anything below this line. All important paths
# can be set above, and you can use CFLAGS and LDFLAGS in the environment.
#
# New numbering scheme: NAME-GMXBASE-MAJOR.MINOR
NAME     := gridcount#
GMXBASE  := gmx4.0
MAJOR    := 1#
MINOR    := 0#

CPPFLAGS += -I$(GMX_INCLUDE_DIR)

ifdef DEBUG
CFLAGS   +=  -DDEBUG -g -Wall  -Wno-unused 
else
CFLAGS   +=  -g -O2 -fomit-frame-pointer -finline-functions \
             -funroll-loops -Wall -Wno-unused 
endif

LDFLAGS  +=  -lm -L$(GMX_LIB_DIR) -lmd -lgmx 

CC      := gcc
LD      := $(CC)
INSTALL := install

LIBGC     := libgridcount.a
LIBGC_NAMES :=  utilgmx xf plt count grid3D xdr_grid
LIBGC_SRC := $(addsuffix .c, $(LIBGC_NAMES))
LIBGC_H   := $(addsuffix .h, $(LIBGC_NAMES)) gridcount.h
LIBGC_OBJ := $(addsuffix .o, $(LIBGC_NAMES))


# g_ri3Dc
G_RI3DC     := g_ri3Dc 
G_RI3DC_SRC := g_ri3Dc.c 
G_RI3DC_H   := 
G_RI3DC_OBJ := g_ri3Dc.o 

# a_ri3Dc
A_RI3DC     := a_ri3Dc 
A_RI3DC_SRC := a_ri3Dc.c 
A_RI3DC_H   := 
A_RI3DC_OBJ := a_ri3Dc.o

# a_gridcalc
A_GRIDCALC     := a_gridcalc 
A_GRIDCALC_SRC := a_gridcalc.c 
A_GRIDCALC_H   := 
A_GRIDCALC_OBJ := a_gridcalc.o

ALL_PROG := $(G_RI3DC) $(A_RI3DC) $(A_GRIDCALC)

DIST_GRIDCOUNT_H   := $(sort  $(A_RI3DC_H) $(A_GRIDCALC_H) \
	$(G_RI3DC_H) $(LIBGC_H))
DIST_GRIDCOUNT_SRC := $(sort  $(A_RI3DC_SRC) $(A_GRIDCALC_SRC) \
	$(G_RI3DC_SRC) $(LIBGC_SRC))
DIST_GRIDCOUNT_MAKE     := Makefile config.guess
DIST_GRIDCOUNT_EXAMPLES := examples/Makefile.grid examples/slicer.pl \
	examples/window_density.pl
DIST_GRIDCOUNT_DOCS     := README INSTALL FAQ CHANGELOG
DIST_GRIDCOUNT_DIRS     := contrib

DIST_NAME      := $(NAME)-$(GMXBASE)-$(MAJOR).$(MINOR)
DIST_DIR       := $(DIST_NAME)
DIST_GRIDCOUNT := $(DIST_NAME).tar.gz 


define usage
\nIn order to compile programs edit the Makefile for paths to the Gromacs\
\nlibraries (release $(GMXBASE)). Then do \
\n   make clean; make PROGRAM\
\nwhere PROGRAM can be one of \`$(ALL_PROG)'.\
\nPerhaps you have to edit variables at the top of the Makefile \
\nto make it work. Programs are statically linked so you can take the \
\nbinaries wherever you like (hopefully...).\
\n\
\nIn order to install all programs, type\
\n   make BIN_DIR=<your_target_dir> install\
\nor set BIN_DIR at the top of the Makefile.\
\n\
\nInstall targets: \
\n   all          compile \`$(ALL_PROG)' (default)\
\n   install      install all compiled programs in BIN_DIR\
\n   clean        clean object files etc\
\n   distclean    remove every generated file\
\n\
\nMakefile switches:\
\nSwitch on on commandline by setting them, eg \`make DEBUG=1 g_ri3Dc' \
\n  DEBUG        compile with debugging, -g and no optimisations\
\n\
\nImportant Makefile parameters:\
\nGMX_TOP_DIR...............$(GMX_TOP_DIR)\
\nGMX_EXEC_PREFIX...........$(GMX_EXEC_PREFIX)\
\nGMX_LIB_DIR...............$(GMX_LIB_DIR)\
\nGMX_INCLUDE_DIR...........$(GMX_INCLUDE_DIR)\
\n\
\nBIN_DIR...................$(BIN_DIR)\
\n\
\nCC........................$(CC)\
\nLD........................$(LD)\
\nCPPFLAGS..................$(CPPFLAGS)\
\nCFLAGS....................$(CFLAGS)\
\nLDFLAGS...................$(LDFLAGS)
endef
# 'emacs font-lock

.PHONY: all help 
all:	$(ALL_PROG)

help:
	@echo -e "$(usage)"


$(LIBGC): $(LIBGC)($(LIBGC_OBJ)) $(LIBGC_H)
	ranlib $@

$(G_RI3DC): $(G_RI3DC_OBJ) $(LIBGC)
	$(LD) -o $@ $^ $(LDFLAGS)
$(G_RI3DC_OBJ): $(G_RI3DC_SRC) $(G_RI3DC_H)


$(A_RI3DC): $(A_RI3DC_OBJ) $(LIBGC) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(A_RI3DC_OBJ): $(A_RI3DC_SRC) $(A_RI3DC_H)


$(A_GRIDCALC): $(A_GRIDCALC_OBJ) $(LIBGC) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(A_GRIDCALC_OBJ): $(A_GRIDCALC_SRC) $(A_GRIDCALC_H)


TAGS: $(ALL_SOURCES)
	etags *.c *.h 
	find $(GMX_SOURCE_DIR) -name '*.[ch]' | xargs etags --append 


.PHONY: clean distclean tar-clean install dist rsync

install:  $(G_RI3DC) $(A_RI3DC) $(A_GRIDCALC)
	for p in $^; do \
	    if [ -e $$p ]; then  \
	       echo ">>> Installing file \`$$p' ..."; \
	       $(INSTALL) -v -m 755 $$p $(BIN_DIR); \
	    fi; \
	done;


clean:
	-rm -f core *.o *.a *~ 

distclean: clean tar-clean
	-rm -f $(ALL_PROG) TAGS

.cvsignore: 
	echo $(ALL_PROG) > $@
	echo "*.log *~ core *.a TAGS tmp" >> $@
	echo "*.xtc *.trr *.tpr *.edr *.ndx" >> $@
	echo "*Test*" >> $@

dist: $(DIST_GRIDCOUNT)

$(DIST_GRIDCOUNT): FORCE
	-rm -rf $(DIST_DIR)
	mkdir $(DIST_DIR) $(DIST_DIR)/examples
	cp $(DIST_GRIDCOUNT_H) $(DIST_GRIDCOUNT_SRC) $(DIST_DIR)
	cp $(DIST_GRIDCOUNT_DOCS) $(DIST_GRIDCOUNT_MAKE) $(DIST_DIR)
	cp $(DIST_GRIDCOUNT_EXAMPLES) $(DIST_DIR)/examples
	cp -r $(DIST_GRIDCOUNT_DIRS) $(DIST_DIR)
	tar -zcvf $@ $(DIST_DIR)
	-rm -r $(DIST_DIR)

tar-clean:
	-rm -rf $(DIST_DIR) $(DIST_GRIDCOUNT)

rsync: $(DIST_GRIDCOUNT)
	rsync -avP $^ clathrin:/sansom/public_html/sbcb/oliver/download/Gromacs

.phony: FORCE
FORCE:
