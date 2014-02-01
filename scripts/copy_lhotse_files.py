#! /usr/bin/env python

# Copies required LHOTSE source files to a new destination, retaining the
# original directory structure. There are two types of files:
# - All files in lh_files: These files are copied here
# - All (visible) files in projdir and below (apart from certain
#   exceptions): These are not copied, but a command is output which does
#   the copying

import os
import shutil
#import subprocess

src_root = '/home/seeger/lhotse'
trg_root = '/home/seeger/apbsint'
projdir = 'src/eptools'
rootadd_dir = os.path.join(projdir,'doc/rootfiles')

tmp_fname = 'apbsint.tgz'
# What is excluded here:
# - Temporary files
#   NOTE: Due to a bug, --exclude-backups does not exclude all *~ files
# - doc/rootfiles (see 'rootadd_dir')
# - Matlab build: Sym link and workaround
# - Python build: Temp. files, setup.py symlink and workaround
# - python/tmp: Contains private stuff
# - C++ files belonging to workaround
tar_ccmd = ("tar cfz " + tmp_fname + " --exclude-vcs --exclude-backups " +
            "--exclude '*~' " +
            "--exclude '*.pyc' --exclude '*.o' --exclude doc/rootfiles " +
            "--exclude doc/design/migrate2github.txt " +
            "--exclude python/cython/build --exclude 'matlab/bin/*' " +
            "--exclude matlab/make.inc.def " +
            "--exclude matlab/make.inc.workaround " +
            "--exclude python/cython/eptools_ext.so " +
            "--exclude python/cython/eptools_ext.cpp " +
            "--exclude python/cython/setup.py " +
            "--exclude python/cython/setup.py.workaround " +
            "--exclude python/tmp " +
            "--exclude scripts/workaround_symlinks.tcl " +
            "--exclude potentials/SpecfunServices_workaround.h " +
            "--exclude doc/howto_workaround.txt " +
            "--exclude potentials/quad *")

tar_xcmd = 'tar xzf ' + tmp_fname

# LHOTSE files (not part of project)
lh_dirs = ['lhotse', 'src', 'lhotse/matif']
lh_files = ['lhotse/FileUtils.cc',
            'lhotse/IntVal.cc',
            'lhotse/Interval.cc',
            'lhotse/Range.cc',
            'lhotse/StandardException.cc',
            'lhotse/global.cc',
            'lhotse/AccumulFunc.h',
            'lhotse/ArrayHandle.h',
            'lhotse/ArrayPtrHandle.h',
            'lhotse/AssertMethod.h',
            'lhotse/DebugVars.h',
            'lhotse/DefaultLogs.h',
            'lhotse/FileUtils.h',
            'lhotse/FuncObjects.h',
            'lhotse/Handle.h',
            'lhotse/IntVal.h',
            'lhotse/Interval.h',
            'lhotse/LogFile.h',
            'lhotse/NullaryFunc.h',
            'lhotse/NumberFormats.h',
            'lhotse/Range.h',
            'lhotse/StandardException.h',
            'lhotse/exceptions.h',
            'lhotse/global.h',
            'lhotse/global_mem.h',
            'lhotse/matif/mex_for_cpp.h',
            'src/main.h']

# Copy LHOTSE files
print 'Copy LHOTSE files'
if not os.path.exists(trg_root):
    os.mkdir(trg_root)
for dn in lh_dirs:
    os.mkdir(os.path.join(trg_root,dn))
for fn in lh_files:
    shutil.copy2(os.path.join(src_root,fn),os.path.join(trg_root,fn))

# The files in 'rootadd_dir' are copied directly to the target root
# NOTE: We assume there are no subdirectories
print 'Copy files from', rootadd_dir
for root, dirs, files in os.walk(os.path.join(src_root,rootadd_dir)):
    if '.svn' in dirs:
        dirs.remove('.svn')
    files = [x for x in files if x[-1] != '~']
    for fn in files:
        shutil.copy2(os.path.join(root,fn),os.path.join(trg_root,fn))

# Copy project files
# Somehow, I cannot call 'tar' from Python properly: giving up.
os.mkdir(os.path.join(trg_root,projdir))
print '\nOK: Now do the following:'
pstr = ('  cd ' + os.path.join(src_root,projdir) + '; ' + tar_ccmd + '; ' +
        'cd ' + os.path.join(trg_root,projdir) +'; ' +
        'mv ' + os.path.join(src_root,projdir,tmp_fname) + ' .; ' +
        tar_xcmd + '; rm ' + tmp_fname)
print pstr

#os.chdir(os.path.join(src_root,projdir))
#subprocess.call(tar_cmd.split())

#for root, dirs, files in os.walk(projdir):
#    if '.svn' in dirs:
#        dirs.remove('.svn')
    # Remove undesired files
#    files = [x for x in files if os.path.splitext(x)[1] not in bad_end]
#    files = [x for x in files if x[-1] != '~']
#    for fn in files:
#        ffn = os.path.join(root,fn)
#        shutil.copy2(os.path.join(src_root,ffn),os.path.join(trg_root,ffn))
    # Create directories
#    for dn in dirs:
#        os.mkdir(os.path.join(trg_root,root,dn))
