FFTWPATH        = /opt/fftw

# ----------

VPATH           = source:test

LIBRARY         = fftw3_wrapper
LIBRARY_SOURCES = fftw3.f90 fftw3_wrapper.f90
LIBRARY_OBJECTS = $(LIBRARY_SOURCES:.f90=.o)
PROGRAM         = fftw3_wrapper_test
PROGRAM_SOURCES = fftw3_wrapper_test.f90
PROGRAM_OBJECTS = $(PROGRAM_SOURCES:.f90=.o)

AR              = ar
ARFLAGS         = -rcDs
CPP             = cpp -C -nostdinc -P
CPPFLAGS        = -DHAVE_LIBFFTW3F -I.
FC              = gfortran

# IEEE 754 Compliance
# https://gcc.gnu.org/wiki/FloatingPointMath

FCFLAGS         = -I. -I$(FFTWPATH)/include -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -fpic -pipe -std=f2008 -Wall -Wextra -Wpedantic
LDFLAGS         = -L$(FFTWPATH)/lib -fpic
LIBS            = -lfftw3 -lfftw3f

ifeq ($(strip $(DEBUG)), 1)
FCFLAGS        := $(FCFLAGS) -fbacktrace -fcheck=all -g -Og
else
FCFLAGS        := $(FCFLAGS) -flto -march=haswell -mtune=native -O3
LDFLAGS        := $(LDFLAGS) -flto -march=haswell -mtune=native -O3
endif

LIBRARY_LDFLAGS = $(LDFLAGS) -shared
PROGRAM_LDFLAGS = $(LDFLAGS)

CP              = cp -f
CPDIR           = cp -f -r
MKDIR           = mkdir -p
RM              = rm -f
RMDIR           = rm -f -r

# ----------

.SUFFIXES:

.PHONY: all

all: lib$(LIBRARY).a lib$(LIBRARY).so
	@ls -l *.a *.mod *.so 2>/dev/null

.PHONY: check

check: all $(PROGRAM).exe
	@PROGRAM="$(PROGRAM)" FFTWPATH="$(FFTWPATH)" "./$(PROGRAM).sh"

.PHONY: clean

clean:
	@$(RM) *.a *.exe *.mod *.o *.so *.txt

# ----------

lib$(LIBRARY).a: $(LIBRARY_OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

lib$(LIBRARY).so: $(LIBRARY_OBJECTS)
	$(FC) $(LIBRARY_LDFLAGS) -o $@ $^ $(LIBS)

$(PROGRAM).exe: $(PROGRAM_OBJECTS) lib$(LIBRARY).a
	$(FC) $(PROGRAM_LDFLAGS) -o $@ $^ $(LIBS)

source/fftw3_wrapper.f90: fftw3_wrapper_template.F90 fftw3_wrapper_forward_implementation.F90 fftw3_wrapper_inverse_implementation.F90
	$(CPP) $(CPPFLAGS) $< >$@

fftw3_wrapper.o: fftw3_wrapper.f90 fftw3.o
	$(FC) $(FCFLAGS) -c -o $@ $<

test/fftw3_wrapper_test.f90: fftw3_wrapper_test_template.F90 fftw3_wrapper_test_implementation.F90
	$(CPP) $(CPPFLAGS) $< >$@

fftw3_wrapper_test.o: fftw3_wrapper_test.f90 fftw3_wrapper.o
	$(FC) $(FCFLAGS) -c -o $@ $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
