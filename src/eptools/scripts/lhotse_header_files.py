#! /usr/bin/env python

# Processes the output of compiling MEX files with -MM option. The preprocessor
# writes all files included into the .o files.
# We collect all, remove duplicates, and output the result

of_lst = ['/home/seeger/lhotse/lhotse/global.o',
          '/home/seeger/lhotse/lhotse/StandardException.o',
          '/home/seeger/lhotse/lhotse/FileUtils.o',
          '/home/seeger/lhotse/lhotse/MachDep.o',
          '/home/seeger/lhotse/lhotse/IntVal.o',
          '/home/seeger/lhotse/lhotse/Interval.o',
          '/home/seeger/lhotse/lhotse/Range.o',
          '/home/seeger/lhotse/lhotse/specfun/Specfun.o',
          '/home/seeger/lhotse/src/eptools/potentials/DefaultPotManager.o',
          '/home/seeger/lhotse/src/eptools/potentials/EPPotentialFactory.o',
          '/home/seeger/lhotse/src/eptools/potentials/EPPotentialNamedFactory.o',
          '/home/seeger/lhotse/src/eptools/potentials/PotManagerFactory.o',
          '/home/seeger/lhotse/src/eptools/FactorizedEPRepresentation.o',
          '/home/seeger/lhotse/src/eptools/FactorizedEPDriver.o',
          '/home/seeger/lhotse/src/eptools/wrap/eptools_helper.o',
          '/home/seeger/lhotse/src/eptools/wrap/eptools_helper_basic.o',
          '/home/seeger/lhotse/src/eptools/matlab/mex/mex_helper.o',
          '/home/seeger/lhotse/src/eptools/wrap/eptwrap_getpotid.o',
          '/home/seeger/lhotse/src/eptools/wrap/eptwrap_choluprk1.o',
          '/home/seeger/lhotse/src/eptools/wrap/eptwrap_choldnrk1.o',
          '/home/seeger/lhotse/src/eptools/matlab/mex/eptools_getpotid.o']

all_files = set()
for fname in of_lst:
    fid = open(fname,'r')
    lst = fid.read().replace('\n','').replace('\\','').split()
    all_files.update(lst[1:])
#print all_files
lst = list(all_files)
lst.sort()
for x in lst:
    print x
