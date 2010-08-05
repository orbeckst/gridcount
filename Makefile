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

# Set GMX_TOP_DIR to your Gromacs installation and check the
# INCLUDE and LIB dir.
GMX_TOP_DIR     := $(HOME)/Library/Gromacs/version/4.0.2

# Set the directories where Gromacs executables and libraries are 
# to be found. GMX_EXEC_PREFIX will probably depend on your OS;
# use "$(ARCH)" if you have versions for different architectures installed:
GMX_EXEC_PREFIX = $(GMX_TOP_DIR)/$(ARCH)
GMX_LIB_DIR     = $(GMX_EXEC_PREFIX)/lib
# Include directory for Gromacs header files:
GMX_INCLUDE_DIR = $(GMX_TOP_DIR)/include/gromacs

# install binaries into
BIN_DIR = $(GMX_EXEC_PREFIX)/bin

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
GMXBASE  := gmx4.0
MAJOR    := 1#
MINOR    := 2#

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

.PHONY: all install help clean dist-clean
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
	rsync -avP $^ clathrin:/sansom/public_html/sbcb/oliver/download/Gromacs

.phony: FORCE
FORCE:
