#! /usr/bin/env python

# EPTOOLS Python Interface
# Test of quadrature implementation of EP updates. Laplace transformation
# and adaptive quadrature.
#
# Potential: Probit. We compare the exact implementation with adaptive
# quadrature and Newton root finding.

import numpy as np
import apbsint as abt

# Helper functions

def reldiff(a,b):
    return np.abs(a-b)/np.maximum(np.maximum(np.abs(a),np.abs(b)),1e-8)

# Main code

# Parameters (potential, quadrature)
# Switch on 'use_parallel' for larger range of cavity moments, with
# verbosity switched off.
use_parallel = True
q_verbose = 0
#use_parallel = False
#q_verbose = 1
prb_y = 1.
prb_soff = 0.
q_maxiv = 50
q_epsabs = 1e-6
q_epsrel = 1e-6
# Cavity moments
#cmu = np.arange(-5.,5.,0.1)
cmu = np.arange(-20.,20.,0.1)
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

pars = np.empty(2)
pars[0] = prb_y; pars[1] = prb_soff
if not use_parallel:
    # Main loop, using epupdate_single
    for i in xrange(m):
        print '\nIteration %d:' % i
        # Exact computation
        (rstat1, alpha1, nu1, logz1) = \
            abt.eptools_ext.epupdate_single('Probit',pars,0,cmu[i],crho[i])
        # Adaptive quadrature (normal)
        (rstat2, alpha2, nu2, logz2) = \
            abt.eptools_ext.epupdate_single('DebugQuadProbit',pars,qs2_ptr,
                                            cmu[i],crho[i])
        print 'Exact:   al=%.4e, nu=%.4e, logz=%.4e (rstat=%d)' % \
            (alpha1,nu1,logz1,rstat1)
        print 'AdaQuad: al=%.4e, nu=%.4e, logz=%.4e (rstat=%d)' % \
            (alpha2,nu2,logz2,rstat2)
else:
    # Setup potential managers
    num_pm = 2
    pm_elem1 = abt.ElemPotManager('Probit',m,tuple(pars))
    potman1 = abt.PotManager(pm_elem1)
    potman1.check_internal()
    pm_elem2 = abt.ElemPotManager('DebugQuadProbit',m,tuple(pars),qserv2)
    potman2 = abt.PotManager(pm_elem2)
    potman2.check_internal()
    # Run EP updates for different PMs
    rstat = np.empty((num_pm,m),dtype=np.int32)
    alpha = np.empty((num_pm,m))
    nu = np.empty((num_pm,m))
    logz = np.empty((num_pm,m))
    pm_lst = (potman1, potman2)
    for i in range(num_pm):
        pm = pm_lst[i]
        abt.eptools_ext.epupdate_parallel(pm.potids,pm.numpot,pm.parvec,
                                          pm.parshrd,pm.annobj,cmu,crho,
                                          rstat[i],alpha[i],nu[i],logz[i])
        mm = np.nonzero(rstat[i])[0].shape[0]
        if mm<m:
            print 'i=%d: %d updates failed' % (i,m-mm)
    # Comparison:
    # For each of logz, alpha, nu, we always show the 3 cases with largest
    # relative difference. We compute rd(exact,quad_normal),
    # rd(exact,quad_extreme), then take the max of both
    mat1 = np.vstack((logz[0], alpha[0], nu[0]))
    mat2 = np.vstack((logz[1], alpha[1], nu[1]))
    names = ('logz', 'alpha', 'nu')
    print 'Probit potential: Exact vs. adaptive quadrature'
    print 'y = %f, soff=%f' % (prb_y, prb_soff)
    print 'Quadrature: maxiv=%d, epsabs=%f, epsrel=%f' % (q_maxiv,q_epsabs,
                                                          q_epsrel)
    for k in range(3):
        v1 = mat1[k]; v2 = mat2[k];
        rdf = reldiff(v1,v2)
        ind = np.argsort(rdf)
        print '%s:' % names[k]
        for j in ind[-1:-4:-1]:
            print ('  rdf=%.4e (v1=%f,v2=%f):' +
                   ' j=%d,cmu=%f,crho=%f') % (rdf[j],v1[j],v2[j],j,cmu[j],
                                              crho[j])
