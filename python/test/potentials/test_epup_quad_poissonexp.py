#! /usr/bin/env python

# EPTOOLS Python Interface
# Test of quadrature implementation of EP updates. Laplace transformation
# and adaptive quadrature.
#
# Potential: Poisson with exponential rate function.
# NOTE: This is not a real test against ground truth, but rather checking
# a consistency constraint between results for different y.

import numpy as np
import apbsint as abt

# Helper functions

def reldiff(a,b):
    return np.abs(a-b)/np.maximum(np.maximum(np.abs(a),np.abs(b)),1e-8)

# Main code

# Parameters (potential, quadrature)
y_lst = [0., 1., 2., 4., 7., 15.]
q_maxiv = 50
q_epsabs = 1e-6
q_epsrel = 1e-6
# Cavity moments
cmu = np.arange(-5.,5.,0.1)
#cmu = np.arange(-20.,20.,0.1)
crho_val = 1.
#cmu = np.arange(-200.,200.,1.)
#crho_val = 1e+3
#crho_val = 1e-6
crho = crho_val*np.ones_like(cmu)

m = cmu.shape[0]

# Quadrature services
qserv2 = abt.ptannotate_ext.AdaptiveQuadPackServices(q_maxiv,q_epsabs,q_epsrel,
                                                     q_verbose)
qs2_ptr = np.uint64(qserv2.getptr())
# Setup potential managers
num_pm = len(y_lst)
potman = []
for y in y_lst:
    pm_elem = abt.ElemPotManager('PoissonExpRate',m,y,qserv2)
    potman.append(abt.PotManager(pm_elem))
    potman[-1].check_internal()
# Run EP updates for different PMs
# We always compare y=0 against one of the others
rstat = np.empty(m,dtype=np.int32)
alpha = np.empty(m)
nu = np.empty(m)
logz = np.empty(m)
cmu2 = np.empty(m)
rstat2 = np.empty(m,dtype=np.int32)
alpha2 = np.empty(m)
nu2 = np.empty(m)
logz2 = np.empty(m)
print 'Poisson potential, exponential link: Consistency check (different y)'
print 'Quadrature: maxiv=%d, epsabs=%f, epsrel=%f' % (q_maxiv,q_epsabs,
                                                      q_epsrel)
for i in range(1,num_pm):
    y = y_lst[i]
    pm = potman[i]
    abt.eptools_ext.epupdate_parallel(pm.potids,pm.numpot,pm.parvec,
                                      pm.parshrd,pm.annobj,cmu,crho,rstat,
                                      alpha,nu,logz)
    mm = np.nonzero(rstat)[0].shape[0]
    if mm<m:
        print 'i=%d: %d updates failed' % (i,m-mm)
    cmu2[:] = cmu
    cmu2 += y*crho
    pm = potman[0]
    abt.eptools_ext.epupdate_parallel(pm.potids,pm.numpot,pm.parvec,
                                      pm.parshrd,pm.annobj,cmu2,crho,rstat2,
                                      alpha2,nu2,logz2)
    mm = np.nonzero(rstat2)[0].shape[0]
    if mm<m:
        print 'i=%d[y=0]: %d updates failed' % (i,m-mm)
    alpha2 += y
    cnst = 0.
    for j in np.arange(2.,y+1.):
        cnst -= np.log(j)
    logz2 += (y*(cmu + 0.5*y*crho) + cnst)
    # Comparison:
    # For each of logz, alpha, nu, we always show the 3 cases with largest
    # relative difference. We compute rd(exact,quad_normal),
    # rd(exact,quad_extreme), then take the max of both
    mat1 = np.vstack((logz, alpha, nu))
    mat2 = np.vstack((logz2, alpha2, nu2))
    names = ('logz', 'alpha', 'nu')
    print '\ny=%d vs. y=0' % np.int32(y)
    for k in range(3):
        v1 = mat1[k]; v2 = mat2[k];
        rdf = reldiff(v1,v2)
        ind = np.argsort(rdf)
        print '%s:' % names[k]
        for j in ind[-1:-4:-1]:
            print ('  rdf=%.4e (v1=%f,v2=%f):' +
                   ' j=%d,cmu=%f,crho=%f') % (rdf[j],v1[j],v2[j],j,cmu[j],
                                              crho[j])
