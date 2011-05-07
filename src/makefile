SRCDIR=$(HOME)/src
IDIR = $(SRCDIR)/include
COMMON_DIR = $(SRCDIR)/common

#CC = mpicc
CC = gcc ## -pg
CPPFLAGS = -I $(IDIR) ## -I /usr/lib/openmpi/include
CFLAGS = -O2 -Wall
LDFLAGS = -L /usr/lib/R/lib -lm -lRmath

OBJS = psi.o emam.o mcmcam.o vb.o vbam.o vbnam.o \
	$(COMMON_DIR)/dan.o $(COMMON_DIR)/rdan.o $(COMMON_DIR)/io.o $(SRCDIR)/iogeno/iogeno.o $(COMMON_DIR)/ran.o ## vb.o digamma.o

EXECS = psi

$(EXECS):			$(OBJS)
psi.o:				psi.h
vb.o:				vb.h psi.h
vbnam.o:			vb.h psi.h
vbam.o:				vb.h psi.h
emam.o:				em.h psi.h
mcmcam.o:			mcmc.h psi.h
$(SRCDIR)/iogeno/iogeno.o:	$(IDIR)/iogeno.h
$(COMMON_DIR)/dan.o:		$(IDIR)/dan.h
$(COMMON_DIR)/rdan.o:		$(IDIR)/rdan.h

$(COMMON_DIR)/io.o:		$(COMMON_DIR)/io.c $(COMMON_DIR)/io-templates.c $(IDIR)/io.h
				$(CC) $(CFLAGS) $(CPPFLAGS) -I $(COMMON_DIR) -c -o $@ $< ## different: io-templates.c is #included

clean:			
				$(RM) $(EXECS) $(OBJS) *.o
