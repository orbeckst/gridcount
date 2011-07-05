# GNU Makefile to compile the grid counter
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

#====================================================================
# 
# THE FOLLOWING VARIABLES MUST BE CHECKED/SET BY THE USER
# -------------------------------------------------------
#
# See the file INSTALL for instructions.
#
# Include directory for Gromacs header files:
GMX_INCLUDE_DIR = /usr/local/gromacs/include/gromacs

# Set the directories where Gromacs libraries are to be found:
GMX_LIB_DIR     = /usr/local/gromacs/lib

# Install binaries into:
BIN_DIR = /usr/local/gromacs/bin
#
#
#====================================================================


# GMX_SOURCE_DIR is only necessary for the creation of etags and 
# can be safely ignored (for compilation it is not important).
GMX_SOURCE_DIR  := 
#
##########################################################################


##########################################################################
# Probably no need to touch anything below this line. All important paths
# can be set above, and you can use CFLAGS and LDFLAGS in the environment.
#
# New numbering scheme: NAME-GMXBASE-MAJOR.MINOR
NAME     := gridcount#
GMXBASE  := gmx4.5#
MAJOR    := 1#
MINOR    := 4#

export PROJECTDIR := $(realpath .)
export INCLUDEDIR := $(PROJECTDIR)/include
export SRCDIR     := $(PROJECTDIR)/src
export SCRIPTSDIR := $(PROJECTDIR)/scripts
export BIN_DIR

ARCH := $(shell $(SCRIPTSDIR)/config.guess)

CPPFLAGS += -I$(INCLUDEDIR) -I$(GMX_INCLUDE_DIR)
export CPPFLAGS

ifdef DEBUG
CFLAGS   +=  -DDEBUG -g -Wall  -Wno-unused 
else
CFLAGS   +=  -g -O2 -fomit-frame-pointer -finline-functions \
             -funroll-loops -Wall -Wno-unused 
endif
export CFLAGS

LDFLAGS  +=  -lm -L$(GMX_LIB_DIR) -lmd -lgmx 
export LDFLAGS

export CC      := gcc
export LD      := $(CC)
export INSTALL := install


DIST_GRIDCOUNT_H   := $(basename $(INCLUDEDIR))
DIST_GRIDCOUNT_SRC := $(basename $(SRCDIR))
DIST_GRIDCOUNT_MAKE     := Makefile $(basename $(SCRIPTSDIR))
DIST_GRIDCOUNT_EXAMPLES := examples/Makefile.grid examples/slicer.pl \
	examples/window_density.pl
DIST_GRIDCOUNT_DOCS     := README INSTALL FAQ CHANGELOG

DIST_NAME      := $(NAME)-$(GMXBASE)-$(MAJOR).$(MINOR)
DIST_DIR       := $(DIST_NAME)
DIST_GRIDCOUNT := $(DIST_NAME).tar.gz 


define usage
\nIn order to compile programs edit the Makefile for paths to the Gromacs\
\nlibraries (release $(GMXBASE)). Then do \
\n   make clean; make all\
\n\
\nPerhaps you have to edit variables at the top of the Makefile \
\nto make it work. Programs are statically linked so you can take the \
\nbinaries wherever you like (hopefully...).\
\n\
\nIn order to install all programs, type\
\n   make BIN_DIR=<your_target_dir> install\
\nor set BIN_DIR at the top of the Makefile.\
\n\
\nInstall targets: \
\n   all          compile everything (default)\
\n   install      install all compiled programs in BIN_DIR\
\n   clean        clean object files etc\
\n   dist-clean   remove every generated file\
\n\
\nMakefile switches:\
\nSwitch on on commandline by setting them, eg \`make DEBUG=1 g_ri3Dc' \
\n  DEBUG        compile with debugging, -g and no optimisations\
\n\
\nImportant Makefile parameters:\
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

.PHONY: all install help clean dist-clean check
all:
	cd $(SRCDIR) && $(MAKE) $@

install:
	cd $(SRCDIR) && $(MAKE) $@

clean:
	cd $(SRCDIR) && $(MAKE) $@

dist-clean: tar-clean
	cd $(SRCDIR) && $(MAKE) $@
	rm -f TAGS

help:
	@echo -e "$(usage)"

define include_check
if test -e $(GMX_INCLUDE_DIR)/xtcio.h; then \
    echo "OK   Gromacs include directory found: $(GMX_INCLUDE_DIR)"; \
else \
    echo "BAD  Gromacs include directory missing. Set GMX_INCLUDE_DIR in Makefile!"; \
fi
endef

_LIBS = $(wildcard $(GMX_LIB_DIR)/libgmx.*)
define lib_check
if test -n "$(_LIBS)"; then \
    echo "OK   Gromacs lib directory found: $(GMX_LIB_DIR)"; \
else \
    echo "BAD  Gromacs lib directory missing. Set GMX_LIB_DIR in Makefile!"; \
fi
endef

check:
	@echo "============================================================"
	@echo "Checking if GMX_INCLUDE_DIR and GMX_LIB_DIR are set properly"
	@echo "============================================================"
	@$(call include_check)
	@$(call lib_check)
	@echo "============================================================"
	@echo "If you got any 'BAD' entries please read INSTALL."

TAGS: FORCE
	find $(PROJECTDIR) -name '*.[ch]' | xargs etags
	find $(GMX_SOURCE_DIR) -name '*.[ch]' | xargs etags --append 

.cvsignore: 
	echo $(ALL_PROG) > $@
	echo "*.log *~ core *.a TAGS tmp" >> $@
	echo "*.xtc *.trr *.tpr *.edr *.ndx" >> $@
	echo "*Test*" >> $@

dist: $(DIST_GRIDCOUNT)

$(DIST_GRIDCOUNT): FORCE
	-rm -rf $(DIST_DIR)
	mkdir $(DIST_DIR) $(DIST_DIR)/examples
	cp -r $(DIST_GRIDCOUNT_H) $(DIST_GRIDCOUNT_SRC) $(DIST_DIR)
	cp -r $(DIST_GRIDCOUNT_DOCS) $(DIST_GRIDCOUNT_MAKE) $(DIST_DIR)
	cp $(DIST_GRIDCOUNT_EXAMPLES) $(DIST_DIR)/examples
	tar --exclude='*~' --exclude='*.o' --exclude='*.a' \
            --exclude=a_gridcalc --exclude=a_ri3Dc --exclude=g_ri3Dc \
            -zcvf $@ $(DIST_DIR)
	-rm -r $(DIST_DIR)

tar-clean:
	-rm -rf $(DIST_DIR) $(DIST_GRIDCOUNT)

rsync: $(DIST_GRIDCOUNT)
	rsync -avP $^ clathrin:/sansom/public_html/html/sbcb/oliver/download/Gromacs

.phony: FORCE
FORCE:
