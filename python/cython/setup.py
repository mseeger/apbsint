#! /usr/bin/env python

# Build ApBsInT extension modules (C++ code). Use '--workaround' option
# in order to build workaround code (this needs private part not contained
# in the public repo).

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import sys
import numpy

work_around = False
if '--workaround' in sys.argv:
    work_around = True
    sys.argv.remove('--workaround')

# Basic information passed to compiler/linker
# - df_include_dirs: Include path(s). The repo root must be in there
# - df_library_dirs: Library path(s) for your system
df_include_dirs = [numpy.get_include(), '/home/seeger/apbsint']
nwa_include_dirs = df_include_dirs[:]
df_define_macros = [('HAVE_NO_BLAS', None), ('HAVE_FORTRAN', None)]
nwa_define_macros = df_define_macros[:]
df_libraries = ['m']
nwa_libraries = df_libraries[:]
df_library_dirs = ['/usr/lib64']
nwa_library_dirs = df_include_dirs[:]

if work_around:
    df_define_macros.append(('HAVE_LIBGSL', None))
    df_libraries.append('gsl')

# eptools_ext: Main API to C++ functions
eptools_ext_sources = [
    'eptools_ext.pyx',
    'base/lhotse/global.cc',
    'base/lhotse/StandardException.cc',
    'base/lhotse/FileUtils.cc',
    'base/lhotse/IntVal.cc',
    'base/lhotse/Interval.cc',
    'base/lhotse/Range.cc',
    'base/lhotse/optimize/OneDimSolver.cc',
    'base/src/eptools/potentials/DefaultPotManager.cc',
    'base/src/eptools/potentials/EPPotentialFactory.cc',
    'base/src/eptools/potentials/EPPotentialNamedFactory.cc',
    'base/src/eptools/potentials/PotManagerFactory.cc',
    'base/src/eptools/potentials/SpecfunServices.cc',
    'base/src/eptools/FactorizedEPRepresentation.cc',
    'base/src/eptools/FactorizedEPDriver.cc',
    'base/src/eptools/wrap/eptools_helper_basic.cc',
    'base/src/eptools/wrap/eptools_helper.cc',
    'base/src/eptools/wrap/eptwrap_choldnrk1.cc',
    'base/src/eptools/wrap/eptwrap_choluprk1.cc',
    'base/src/eptools/wrap/eptwrap_epupdate_parallel.cc',
    'base/src/eptools/wrap/eptwrap_epupdate_single.cc',
    'base/src/eptools/wrap/eptwrap_fact_compmarginals.cc',
    'base/src/eptools/wrap/eptwrap_fact_compmaxpi.cc',
    'base/src/eptools/wrap/eptwrap_fact_sequpdates.cc',
    'base/src/eptools/wrap/eptwrap_getpotid.cc',
    'base/src/eptools/wrap/eptwrap_getpotname.cc',
    'base/src/eptools/wrap/eptwrap_potmanager_isvalid.cc'
]
if work_around:
    eptools_ext_sources.append(
        'base/lhotse/specfun/Specfun.cc'
    )

# apbtest_ext: API for test code (excluding the workaround)
apbtest_ext_sources = [
    'apbtest_ext.pyx',
    'base/lhotse/global.cc',
    'base/lhotse/StandardException.cc',
    'base/lhotse/FileUtils.cc',
    'base/lhotse/IntVal.cc',
    'base/lhotse/Interval.cc',
    'base/lhotse/Range.cc',
    'base/lhotse/optimize/OneDimSolver.cc',
    'base/src/eptools/potentials/SpecfunServices.cc'
]

# apbtest_workaround_ext: API for test code (workaround part)
apbtest_workaround_ext_sources = [
    'apbtest_workaround_ext.pyx',
    'base/lhotse/global.cc',
    'base/lhotse/StandardException.cc',
    'base/lhotse/FileUtils.cc',
    'base/lhotse/IntVal.cc',
    'base/lhotse/Interval.cc',
    'base/lhotse/Range.cc',
    'base/lhotse/specfun/Specfun.cc'
]

# ptannotate_ext: Potential annotation classes. Right now, they are all in
#     the workaround
ptannotate_ext_sources = [
    'ptannotate_ext.pyx',
    'base/lhotse/global.cc',
    'base/lhotse/StandardException.cc',
    'base/lhotse/FileUtils.cc',
    'base/lhotse/IntVal.cc',
    'base/lhotse/Interval.cc',
    'base/lhotse/Range.cc',
    'base/src/eptools/potentials/quad/AdaptiveQuadPackServices.cc',
    'base/src/eptools/potentials/quad/AdaptiveQuadPackDebugServices.cc'
]

df_ext_modules = [
    Extension(
        'eptools_ext',
        sources = eptools_ext_sources,
        include_dirs = df_include_dirs,
        define_macros = df_define_macros,
        libraries = df_libraries,
        library_dirs = df_library_dirs,
        language = 'c++'
    ),
    # NOTE: apbtest_ext has to be build in default mode, even if the workaround
    # is used. This is because we use apbtest_ext vs. apbtest_workaround_ext
    # for comparison tests.
    Extension(
        'apbtest_ext',
        sources = apbtest_ext_sources,
        include_dirs = nwa_include_dirs,
        define_macros = nwa_define_macros,
        libraries = nwa_libraries,
        library_dirs = nwa_library_dirs,
        language = 'c++'
    )
]
if work_around:
    df_ext_modules.append(
        Extension(
            'apbtest_workaround_ext',
            sources = apbtest_workaround_ext_sources,
            include_dirs = df_include_dirs,
            define_macros = df_define_macros,
            libraries = df_libraries,
            library_dirs = df_library_dirs,
            language = 'c++'
        )
    )
    df_ext_modules.append(
        Extension(
            'ptannotate_ext',
            sources = ptannotate_ext_sources,
            include_dirs = df_include_dirs,
            define_macros = df_define_macros,
            libraries = df_libraries,
            library_dirs = df_library_dirs,
            language = 'c++'
        )
    )

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = df_ext_modules
)
