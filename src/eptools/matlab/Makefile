# -------------------------------------------------------------------
# ApBsInT: Approximate Bayesian Inference Toolkit
# -------------------------------------------------------------------
# Makefile for Matlab interface (MEX functions)
# Author: Matthias Seeger
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Customization:
# System-dependent variables are set in the file make.inc.def, or
# alternatively in the file make.inc.$$$ if "make" is called with
# the option "where=$$$". See LHOTSE installation guide for
# details.
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Control variables (pass on command line)
#
# - where: Where does the compilation take place?
#   Use this variable if you maintain the LHOTSE system in different
#   environments in parallel. For a legal value XXX of 'where', a
#   file make.inc.XXX must be present.
#   The default is 'def'. make.inc.def must be present.
# - opt: Level of optimisation:
#		debug:  No optimisation, generate debug info
#		some:   Some optimisation
#		more:   More optimisation. Inlining enforced
#		all:	Highest generic level of compiler. Inlining
#                       enforced
#   The default is 'debug'. NOTE: In level 'debug', the compiler does not
#   do inlining, and the matrix/vector methods will run very slow. Do not
#   use this level for final experiments!
# - mex: Used only for targets resulting in Matlab MEX libraries being
#   built. By def., the value is 'no', and a built will fail when using
#   "mex". This is because LHOTSE enforces static linkage. Values:
#   		yes:	  Use this with targets res. in MEX libraries.
#			  Will try to use dynamic libraries throughout
#               yes_nomm: Same as 'yes', but USE_MATLAB_MM is not defined,
#                         so STL allocators are not replaced. This is
#                         simpler and works whenever memory does not have
#                         to be kept persistent between MEX calls
# - fort: By def., the value is 'yes'. If 'no', none of the external
#   Fortran code is compiled. The corr. LHOTSE wrappers will not work,
#   in general an 'NotImplemException' is thrown. The macro HAVE_FORTRAN
#   is defined iff 'fort' is 'yes'
#
# Generic targets:
#
# - clean:     Remove all objects, but not the libraries
# - clean_all: Remove all objects and libraries
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Default values for some of the system-dependent variables
# -------------------------------------------------------------------

CC=		gcc
CXX=		g++
FC=		g77
#FORCE_STATIC=	-static
MEXCMD=		mex
MEXSUFFIX=	mexa64

EXLIBS_BLAS=	-llapack -lf77blas -lcblas -latlas -lg2c
EXLDOPTS_BLAS=

EXLIBS_GSL=	-lgsl
EXLDOPTS_GSL=
EXINCS_GSL=
DEFINES_GSL=	-DHAVE_LIBGSL

EXINCS_MATLAB=
EXLDOPTS_MATLAB=	-cxx

# -------------------------------------------------------------------
# Include make.inc.$(where) file with system-dependent variables
# -------------------------------------------------------------------

blas=	no
mex=	yes_nomm
where=	def
prof=	no
opt=	debug
fort=	no
gsl=	no

include make.inc.$(where)

# -------------------------------------------------------------------
# Most relevant variables:
# GCCOPTS:	  CPP compiler options
# DEFINES:	  Macro def. passed to compiler
# INCS:		  Include path directives ("-I...")
# LIBS:		  Libraries for linking (excludes the LHOTSE main
# 		  library, which is not used for all targets)
# LDFLAGS:	  Linker flags, library path directives
# These are assembled based on the control variables
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Settings depending on control variables.
# Most of these are additive, they are appended to the corr.
# variable
# -------------------------------------------------------------------

GCCOPTS_optdebug=	-g
GCCOPTS_optsome=	-O
GCCOPTS_optmore=	-O2 -funroll-loops
GCCOPTS_optall=		-O3 -funroll-loops

FFLAGS_optdebug=	-c -g -fno-f2c
FFLAGS_optsome=		-c -fno-f2c -O3
FFLAGS_optmore=		-c -funroll-all-loops -fno-f2c -O3
FFLAGS_optall=		-c -funroll-all-loops -fno-f2c -O3

DEFINES_optdebug=	-DHAVE_DEBUG
DEFINES_optsome=
DEFINES_optmore=	-DHAVE_INLINE
DEFINES_optall=		-DHAVE_INLINE

DEFINES_blasno=		-DHAVE_NO_BLAS
GCCOPTS_profyes=	-pg
LDFLAGS_profyes=	-pg

DEFINES_mexyes=		-DMATLAB_MEX -DUSE_MATLAB_MM -DHAVE_BLAS64BIT
DEFINES_mexyes_nomm=	-DMATLAB_MEX -DHAVE_BLAS64BIT
GCCOPTS_mexyes=		-fPIC

# -------------------------------------------------------------------
# Main LHOTSE directories
# -------------------------------------------------------------------

LHOTSELIBDIR=	$(ROOTDIR)/lhotse
LHOTSEPROJDIR=	$(ROOTDIR)/src
LHOTSEEXECDIR=	$(ROOTDIR)

# -------------------------------------------------------------------
# Most high-level definitions
# -------------------------------------------------------------------

RANLIB=	ranlib

CPPFLAGS=	$(DEFINES) $(INCS) $(GCCOPTS)

CFLAGS=		$(DEFINES) $(INCS) $(GCCOPTS)

INCS=		-I$(ROOTDIR) $(CXXINCS)

LDFLAGS=	-L$(LHOTSELIBDIR) -L$(LHOTSEPROJDIR) $(CXXLDOPTS)

LIBS=		-lm

# -------------------------------------------------------------------
# Module objects
# -------------------------------------------------------------------

GLOBALDIR=	$(LHOTSELIBDIR)

SPECFUNDIR=	$(LHOTSELIBDIR)/specfun

# Bare minimum for compiling basic code (f.ex. basic MEX files). In
# particular, the matrix/vector classes are avoided.

ESSMINIMUMOBJS=	$(GLOBALDIR)/global.o \
		$(GLOBALDIR)/StandardException.o \
		$(GLOBALDIR)/FileUtils.o \
		$(GLOBALDIR)/IntVal.o \
		$(GLOBALDIR)/Interval.o \
		$(GLOBALDIR)/Range.o

ESSMINIMUMOBJS_gslyes=	$(SPECFUNDIR)/Specfun.o

# -------------------------------------------------------------------
# Project objects
# A project consists of one or more project modules and a main_XXX.cc
# code file. Objects of project modules are not written into the main
# library.
# -------------------------------------------------------------------

EPTOOLSDIR=	$(LHOTSEPROJDIR)/eptools
_EPTOOLSOBJS=	potentials/DefaultPotManager \
		potentials/EPPotentialFactory \
		potentials/EPPotentialNamedFactory \
		potentials/PotManagerFactory \
		potentials/SpecfunServices \
		FactorizedEPRepresentation \
		FactorizedEPDriver
EPTOOLSOBJS=	$(_EPTOOLSOBJS:%=$(EPTOOLSDIR)/%.o)

EPTOOLSMEXDIR=	$(EPTOOLSDIR)/matlab/mex
_EPTOOLSMEXOBJS=	mex_helper
EPTOOLSMEXOBJS=	$(_EPTOOLSMEXOBJS:%=$(EPTOOLSMEXDIR)/%.o)

EPTOOLSWRAPDIR=	$(EPTOOLSDIR)/wrap
_EPTOOLSWRAPOBJS=	eptools_helper \
			eptools_helper_basic
EPTOOLSWRAPOBJS=	$(_EPTOOLSWRAPOBJS:%=$(EPTOOLSWRAPDIR)/%.o)

# MEX library files are written here:
EPTOOLSLIBDIR=	$(EPTOOLSDIR)/matlab/bin

# -------------------------------------------------------------------
# The make process:
#
# 1. The main target <main> defines TARGET as the name of the
#     internal target, which typically is <main>_int, then runs:
# 2. make_opt$(opt)   (compiler optimization)
# 3. make_blas$(blas) (use BLAS, LAPACK?)
# 4. make_prof$(prof) (profiling by gprof?)
# 5. make_fort$(fort) (use external Fortran code?)
# 6. make_gsl$(gsl)
# 9. make_mex$(mex)   (MEX library?)
#
# NOTE: When extending this, make sure that make_mex is always
# called last, otherwise the mex=stat hack does not work.
# -------------------------------------------------------------------

all:	eptools_all eptools_native_all

# -------------------------------------------------------------------
# Main targets
# 'mainlib' creates the main LHOTSE library, which all projects
# require. The other main targets correspond to projects.
# -------------------------------------------------------------------

# ApBsInT MEX functions (eptools_all compiles all of them, except below)
# Compile with 'mex=yes_nomm blas=no', to avoid linking in large libraries
# and any advanced code.
# 'yes_nomm' avoids the use of lhotse/matif/MatlabAlloc, which is a bit
# shaky anyway.
# NOTE: This excludes MEX functions which require BLAS/LAPACK.
# See eptools_choluprk1 below.

eptools_getpotid:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_getpotname:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_epupdate_single:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_epupdate_parallel:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_fact_sequpdates:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_fact_compmarginals:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_fact_compmaxpi:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_potmanager_isvalid:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_all:	eptools_getpotid eptools_getpotname eptools_epupdate_single \
	eptools_epupdate_parallel eptools_fact_sequpdates \
	eptools_fact_compmarginals eptools_fact_compmaxpi \
	eptools_potmanager_isvalid

# ApBsInT MEX functions which require BLAS/LAPACK.
# Compile with 'mex=yes_nomm blas=no'.
# These are stand-alone C files which do not use LHOTSE, and where $(MEXCMD)
# sorts out everything.
# They are not included in eptools_all, but collected in eptools_native_all.

eptools_choluprk1:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_choldnrk1:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

eptools_native_all:	eptools_choluprk1 eptools_choldnrk1

# -------------------------------------------------------------------
# 'opt'-specific   make commands
# 'prof'-specific  make commands
# 'fort'-specific  make commands
# 'gsl'-specific   make commands
# 'mex'-specific   make commands
# -------------------------------------------------------------------

make_optdebug:
	@$(MAKE) make_blas$(blas) GCCOPTS="$(GCCOPTS_optdebug)" DEFINES="$(DEFINES_optdebug)" FFLAGS="$(FFLAGS_optdebug)"

make_optsome:
	@$(MAKE) make_blas$(blas) GCCOPTS="$(GCCOPTS_optsome)" DEFINES="$(DEFINES_optsome)" FFLAGS="$(FFLAGS_optsome)"

make_optmore:
	@$(MAKE) make_blas$(blas) GCCOPTS="$(GCCOPTS_optmore)" DEFINES="$(DEFINES_optmore)" FFLAGS="$(FFLAGS_optmore)"

make_optall:
	@$(MAKE) make_blas$(blas) GCCOPTS="$(GCCOPTS_optall)" DEFINES="$(DEFINES_optall)" FFLAGS="$(FFLAGS_optall)"

make_blasno:
	@$(MAKE) make_prof$(prof) DEFINES="$(DEFINES) $(DEFINES_blasno)"

#make_blasyes:
#	@$(MAKE) make_prof$(prof) LDFLAGS="$(EXLDOPTS_BLAS) $(LDFLAGS)" LIBS="$(EXLIBS_BLAS) $(LIBS)"

make_profno:
	@$(MAKE) make_fort$(fort)

make_profyes:
	@$(MAKE) make_fort$(fort) GCCOPTS="$(GCCOPTS) $(GCCOPTS_profyes)" LDFLAGS="$(LDFLAGS) $(LDFLAGS_profyes)"

make_fortno:
	@$(MAKE) make_gsl$(gsl)

#make_fortyes:
#	@$(MAKE) make_gsl$(gsl) DEFINES="$(DEFINES) -DHAVE_FORTRAN" OBJECTS="$(OBJECTS) $(FORTRANOBJS)" ESSENTIALOBJS="$(ESSENTIALOBJS) $(FORTRANOBJS)"

make_gslno:
	@$(MAKE) make_mex$(mex)

make_gslyes:
	@$(MAKE) make_mex$(mex) LIBS="$(EXLIBS_GSL) $(LIBS)" LDFLAGS="$(LDFLAGS) $(EXLDOPTS_GSL)" INCS="$(INCS) $(EXINCS_GSL)" DEFINES="$(DEFINES) $(DEFINES_GSL)" ESSMINIMUMOBJS="$(ESSMINIMUMOBJS) $(ESSMINIMUMOBJS_gslyes)"

#make_fftno:
#	@$(MAKE) make_fftw$(fftw)

#make_fftyes:
#	@$(MAKE) make_fftw$(fftw) LIBS="$(EXLIBS_FFT) $(LIBS)" LDFLAGS="$(LDFLAGS) $(EXLDOPTS_FFT)" INCS="$(INCS) $(EXINCS_FFT)" DEFINES="$(DEFINES) $(DEFINES_FFT)" OBJECTS="$(OBJECTS) $(EXTNFFTOBJS)" ESSENTIALOBJS="$(ESSENTIALOBJS) $(EXTNFFTOBJS)"

#make_fftwno:
#	@$(MAKE) make_gpu$(gpu)

#make_fftwyes:
#	@$(MAKE) make_gpu$(gpu) LIBS="$(EXLIBS_FFTW) $(LIBS)" LDFLAGS="$(LDFLAGS) $(EXLDOPTS_FFTW)" INCS="$(INCS) $(EXINCS_FFTW)" DEFINES="$(DEFINES) $(DEFINES_FFTW)"

#make_gpuno:
#	@$(MAKE) make_mex$(mex)

make_mexno:
	@$(MAKE) $(TARGET)
#	@$(MAKE) LDFLAGS="$(FORCE_STATIC) $(LDFLAGS)" $(TARGET)

make_mexyes:
	@$(MAKE) DEFINES="$(DEFINES) $(DEFINES_mexyes)" GCCOPTS="$(GCCOPTS) $(GCCOPTS_mexyes)" INCS="$(INCS) $(EXINCS_MATLAB)" LDFLAGS="$(LDFLAGS) $(EXLDOPTS_MATLAB)" FFLAGS="$(FFLAGS) $(EXINCS_MATLAB)" $(TARGET)

make_mexyes_nomm:
	@$(MAKE) DEFINES="$(DEFINES) $(DEFINES_mexyes_nomm)" GCCOPTS="$(GCCOPTS) $(GCCOPTS_mexyes)" INCS="$(INCS) $(EXINCS_MATLAB)" LDFLAGS="$(LDFLAGS) $(EXLDOPTS_MATLAB)" FFLAGS="$(FFLAGS) $(EXINCS_MATLAB)" $(TARGET)

# -------------------------------------------------------------------
# Internal main targets (MEX libraries)
# -------------------------------------------------------------------

eptools_getpotid_int: $(ESSMINIMUMOBJS) $(EPTOOLSOBJS) $(EPTOOLSWRAPOBJS) $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptwrap_getpotid.o $(EPTOOLSMEXDIR)/eptools_getpotid.o
	$(MEXCMD) -o $(EPTOOLSLIBDIR)/eptools_getpotid.$(MEXSUFFIX) $^ $(LDFLAGS) $(LIBS)

eptools_getpotname_int: $(ESSMINIMUMOBJS) $(EPTOOLSOBJS) $(EPTOOLSWRAPOBJS) $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptwrap_getpotname.o $(EPTOOLSMEXDIR)/eptools_getpotname.o
	$(MEXCMD) -o $(EPTOOLSLIBDIR)/eptools_getpotname.$(MEXSUFFIX) $^ $(LDFLAGS) $(LIBS)

eptools_epupdate_single_int: $(ESSMINIMUMOBJS) $(EPTOOLSOBJS) $(EPTOOLSWRAPOBJS) $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptwrap_epupdate_single.o $(EPTOOLSMEXDIR)/eptools_epupdate_single.o
	$(MEXCMD) -o $(EPTOOLSLIBDIR)/eptools_epupdate_single.$(MEXSUFFIX) $^ $(LDFLAGS) $(LIBS)

eptools_epupdate_parallel_int: $(ESSMINIMUMOBJS) $(EPTOOLSOBJS) $(EPTOOLSWRAPOBJS) $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptwrap_epupdate_parallel.o $(EPTOOLSMEXDIR)/eptools_epupdate_parallel.o
	$(MEXCMD) -o $(EPTOOLSLIBDIR)/eptools_epupdate_parallel.$(MEXSUFFIX) $^ $(LDFLAGS) $(LIBS)

eptools_fact_sequpdates_int: $(ESSMINIMUMOBJS) $(EPTOOLSOBJS) $(EPTOOLSWRAPOBJS) $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptwrap_fact_sequpdates.o $(EPTOOLSMEXDIR)/eptools_fact_sequpdates.o
	$(MEXCMD) -o $(EPTOOLSLIBDIR)/eptools_fact_sequpdates.$(MEXSUFFIX) $^ $(LDFLAGS) $(LIBS)

eptools_fact_compmarginals_int: $(ESSMINIMUMOBJS) $(EPTOOLSOBJS) $(EPTOOLSWRAPOBJS) $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptwrap_fact_compmarginals.o $(EPTOOLSMEXDIR)/eptools_fact_compmarginals.o
	$(MEXCMD) -o $(EPTOOLSLIBDIR)/eptools_fact_compmarginals.$(MEXSUFFIX) $^ $(LDFLAGS) $(LIBS)

eptools_fact_compmaxpi_int: $(ESSMINIMUMOBJS) $(EPTOOLSOBJS) $(EPTOOLSWRAPOBJS) $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptwrap_fact_compmaxpi.o $(EPTOOLSMEXDIR)/eptools_fact_compmaxpi.o
	$(MEXCMD) -o $(EPTOOLSLIBDIR)/eptools_fact_compmaxpi.$(MEXSUFFIX) $^ $(LDFLAGS) $(LIBS)

eptools_potmanager_isvalid_int: $(ESSMINIMUMOBJS) $(EPTOOLSOBJS) $(EPTOOLSWRAPOBJS) $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptwrap_potmanager_isvalid.o $(EPTOOLSMEXDIR)/eptools_potmanager_isvalid.o
	$(MEXCMD) -o $(EPTOOLSLIBDIR)/eptools_potmanager_isvalid.$(MEXSUFFIX) $^ $(LDFLAGS) $(LIBS)

eptools_choluprk1_int: $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptools_helper_basic.o $(EPTOOLSWRAPDIR)/eptwrap_choluprk1.o $(EPTOOLSMEXDIR)/eptools_choluprk1.cc
	$(MEXCMD) -v -largeArrayDims -o $(EPTOOLSLIBDIR)/eptools_choluprk1.$(MEXSUFFIX) $^ $(DEFINES) $(INCS) $(LDFLAGS) -lmwblas -lm

eptools_choldnrk1_int: $(EPTOOLSMEXOBJS) $(EPTOOLSWRAPDIR)/eptools_helper_basic.o $(EPTOOLSWRAPDIR)/eptwrap_choldnrk1.o $(EPTOOLSMEXDIR)/eptools_choldnrk1.cc
	$(MEXCMD) -v -largeArrayDims -o $(EPTOOLSLIBDIR)/eptools_choldnrk1.$(MEXSUFFIX) $^ $(DEFINES) $(INCS) $(LDFLAGS) -lmwblas -lm

# -------------------------------------------------------------------
# Clean targets
# -------------------------------------------------------------------

clean:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

clean_all:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

clean_int:	clean_global clean_specfun clean_src

clean_src:	clean_PROJECTMAIN clean_eptools

clean_all_int: clean_int clean_MEXLIBS

CLEAN_FILES=	-f *.o *.flc *.rpo core errs *~ \#*

clean_MEXLIBS:
	@cd $(EPTOOLSLIBDIR); \
	rm $(CLEAN_FILES); \
	rm *.$(MEXSUFFIX); \
	cd $(ROOTDIR)

clean_global:
	@cd $(GLOBALDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_specfun:
	@cd $(SPECFUNDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_PROJECTMAIN:
	@cd $(LHOTSEPROJDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_eptools:
	@cd $(EPTOOLSDIR); \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/potentials; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/potentials/quad; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/wrap; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/scripts; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/doc; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/doc/design; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/doc/python; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/+ept; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/+ept/@Mat; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/+ept/@MatContainer; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/+ept/@MatDef; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/+ept/@MatDiag; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/+ept/@MatEye; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/+ept/@MatSub; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/mex; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/test; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/test/basic; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/matlab/test/binclass; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/python; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/python/cython; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/python/apbsint; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/python/test; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/python/test/basic; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/python/test/binclass; \
	rm $(CLEAN_FILES); \
	cd $(EPTOOLSDIR)/python/test/potentials; \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)
