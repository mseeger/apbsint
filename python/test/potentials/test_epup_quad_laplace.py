#! /usr/bin/env python

# EPTOOLS Python Interface
# Test of quadrature implementation of EP updates. Laplace transformation
# and adaptive quadrature.
#
# Potential: Laplace. We compare the exact implementation with two
# different adaptive quadrature variants.

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
lap_y = 0.
lap_tau = 1.
q_maxiv = 50
q_epsabs = 1e-6
q_epsrel = 1e-6
# Cavity moments
#cmu = np.arange(-2.,2.,0.1)
cmu = np.arange(-20.,20.,0.1)
crho_val = 1.
#cmu = np.arange(-200.,200.,1.)
#crho_val = 1e+3
#crho_val = 1e-6
crho = crho_val*np.ones_like(cmu)

m = cmu.shape[0]

# Quadrature services
# NOTE: We use q_maxiv/2 for 'AdaptiveQuadPackDebugServices', because this
# runs two adaptive quad calls
qserv2 = abt.ptannotate_ext.AdaptiveQuadPackServices(q_maxiv,q_epsabs,q_epsrel,
                                                     q_verbose)
qs2_ptr = np.uint64(qserv2.getptr())
qserv3 = abt.ptannotate_ext.AdaptiveQuadPackDebugServices(q_maxiv/2,q_epsabs,
                                                          q_epsrel,q_verbose)
qs3_ptr = np.uint64(qserv3.getptr())

pars = np.empty(2)
pars[0] = lap_y; pars[1] = lap_tau
if not use_parallel:
    # Main loop, using epupdate_single
    for i in xrange(m):
        print '\nIteration %d:' % i
        # Exact computation
        (rstat1, alpha1, nu1, logz1) = \
            abt.eptools_ext.epupdate_single('Laplace',pars,0,cmu[i],crho[i])
        # Adaptive quadrature (normal)
        (rstat2, alpha2, nu2, logz2) = \
            abt.eptools_ext.epupdate_single('DebugQuadLaplace',pars,qs2_ptr,
                                            cmu[i],crho[i])
        # Adaptive quadrature (extreme)
        (rstat3, alpha3, nu3, logz3) = \
            abt.eptools_ext.epupdate_single('DebugQuadLaplace',pars,qs3_ptr,
                                            cmu[i],crho[i])
        print 'Exact:          al=%.4e, nu=%.4e, logz=%.4e (rstat=%d)' % \
            (alpha1,nu1,logz1,rstat1)
        print 'Quad (normal):  al=%.4e, nu=%.4e, logz=%.4e (rstat=%d)' % \
            (alpha2,nu2,logz2,rstat2)
        print 'Quad (extreme): al=%.4e, nu=%.4e, logz=%.4e (rstat=%d)' % \
            (alpha2,nu2,logz2,rstat2)
else:
    # Setup potential managers
    pm_elem1 = abt.ElemPotManager('Laplace',m,tuple(pars))
    potman1 = abt.PotManager(pm_elem1)
    potman1.check_internal()
    pm_elem2 = abt.ElemPotManager('DebugQuadLaplace',m,tuple(pars),qserv2)
    potman2 = abt.PotManager(pm_elem2)
    potman2.check_internal()
    pm_elem3 = abt.ElemPotManager('DebugQuadLaplace',m,tuple(pars),qserv3)
    potman3 = abt.PotManager(pm_elem3)
    potman3.check_internal()
    # Run EP updates for 3 different PMs
    rstat = np.empty((3,m),dtype=np.int32)
    alpha = np.empty((3,m))
    nu = np.empty((3,m))
    logz = np.empty((3,m))
    pm_lst = (potman1, potman2, potman3)
    for i in range(3):
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
    mat3 = np.vstack((logz[2], alpha[2], nu[2]))
    names = ('logz', 'alpha', 'nu')
    print 'Laplace potential: Exact vs. quadrature (two versions)'
    print 'y = %f, tau=%f' % (lap_y, lap_tau)
    print 'Quadrature: maxiv=%d, epsabs=%f, epsrel=%f' % (q_maxiv,q_epsabs,
                                                          q_epsrel)
    for k in range(3):
        v1 = mat1[k]; v2 = mat2[k]; v3 = mat3[k];
        rdf1 = reldiff(v1,v2)
        rdf2 = reldiff(v1,v3)
        rdf = np.maximum(rdf1,rdf2)
        ind = np.argsort(rdf)
        print '%s:' % names[k]
        for j in ind[-1:-4:-1]:
            print ('  rdf12=%.4e, rdf13=%.4e (v1=%f,v2=%f,v3=%f):' +
                   ' j=%d,cmu=%f,crho=%f') % (rdf1[j],rdf2[j],v1[j],v2[j],
                                              v3[j],j,cmu[j],crho[j])
