# Makefile in gridcount/src
# Do not run this directly. Instead, go to the top directory
# and use the top level Makefile.


LIBGC     := libgridcount.a
LIBGC_NAMES :=  utilgmx xf plt count grid3D xdr_grid
LIBGC_SRC := $(addsuffix .c, $(LIBGC_NAMES))
LIBGC_H   := $(addprefix $(INCLUDEDIR)/, $(addsuffix .h, $(LIBGC_NAMES)) gridcount.h)
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

.phony: all
all: $(ALL_PROG)

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


.PHONY: clean dist-clean install

install:  $(G_RI3DC) $(A_RI3DC) $(A_GRIDCALC)
	test -d $(BIN_DIR) || mkdir -p $(BIN_DIR)
	for p in $^; do \
	    if [ -e $$p ]; then  \
	       echo ">>> Installing file \`$$p' ..."; \
	       $(INSTALL) -v -m 755 $$p $(BIN_DIR); \
	    fi; \
	done;


clean:
	-rm -f core *.o *.a *~ 

dist-clean: clean
	-rm -f $(ALL_PROG)

