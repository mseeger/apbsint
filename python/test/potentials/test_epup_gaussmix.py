#! /usr/bin/env python

# EPTOOLS Python Interface
# Unit tests: EP updates with Gaussian mixture potential
#
# TODO:
# - Create a single comparison function, which can be used to debug
#   the updates under realistic conditions

import numpy as np
import apbsint as abt

# Helper functions

# 'x' must be a matrix. The logsumexp function is computed for each column,
# and a row vector is returned
def logsumexp(x):
    mx = np.max(x,0)
    return np.log(np.sum(np.exp(x-mx),0)) + mx

def reldiff(a,b):
    return np.abs(a-b)/np.maximum(np.maximum(np.abs(a),np.abs(b)),1e-8)

# Main code

# Specify mixture component parameters and cavity moments
gm_p = np.array([0.25, 0.5, 0.25])
#gm_v = np.array([0.01, 1., 100.])
gm_v = np.array([0.0001, 1., 10000.])
#cmu = np.arange(-20.,20.,0.1)
#crho_val = 1.
cmu = np.arange(-200.,200.,1.)
crho_val = 1e+6
#crho_val = 1e-6
crho = crho_val*np.ones_like(cmu)

numl = gm_p.shape[0]
m = cmu.shape[0]

# Setup potential manager
pvec = np.zeros(2*numl)
pvec[0] = numl
tvec = np.log(gm_p)
pvec[1:numl] = tvec[:numl-1]-tvec[numl-1]
pvec[numl:] = gm_v
pm_elem = abt.ElemPotManager('GaussMixture',m,tuple(pvec))
potman = abt.PotManager(pm_elem)
potman.check_internal()

# Compute log Z and new moments via 'epupdate_parallel'
rstat = np.empty(m,dtype=np.int32)
alpha = np.empty(m)
nu = np.empty(m)
logz = np.empty(m)
abt.eptools_ext.epupdate_parallel(potman.potids,potman.numpot,potman.parvec,
                                  potman.parshrd,cmu,crho,rstat,alpha,nu,logz)
indok = np.nonzero(rstat)[0]
if indok.shape[0] < m:
    print 'epupdate_parallel: %d updates failed' % (m-indok.shape[0])
tvec = 1. - nu[indok]*crho[indok]
indok2 = np.nonzero(tvec >= 1e-12)[0]
m2 = indok2.shape[0]
if m2 < m:
    print 'epupdate_parallel: %d updates give invalid results' % (m-m2)
    indok = indok[indok2].copy()
    hmu = cmu[indok].copy()
    hrho = crho[indok].copy()
    hmu += alpha[indok]*hrho
    hrho *= tvec[indok2]
    logz = logz[indok].copy()
else:
    hmu = cmu + alpha*crho
    hrho = crho*tvec

# Compute log Z and new moments by direct code
# We do everything with numl-by-m2 matrices and NumPy broadcasting.
if m2 < m:
    cmu2 = cmu[indok]
    crho2 = crho[indok]
else:
    cmu2 = cmu
    crho2 = crho
# Z_l = p_l N(0 | cmu, crho + v_l), Z = sum_l Z_l, r_l = Z_l/Z
#import pdb
#pdb.set_trace()
tmat = np.reshape(crho2,(1,m2)) + np.reshape(gm_v,(numl,1))
logzl = -0.5*(np.log(tmat) + np.reshape(cmu2**2,(1,m2))/tmat +
              np.log(2.*np.pi)) + np.reshape(np.log(gm_p),(numl,1))
logz2 = logsumexp(logzl)
rprob = np.exp(logzl-logz2)
# rho_l = crho/(1 + crho/v_l), mu_l = cmu/(1 + crho/v_l)
tmat = np.reshape(crho2,(1,m2))/np.reshape(gm_v,(numl,1)) + 1.
rhol = np.reshape(crho2,(1,m2))/tmat
mul = np.reshape(cmu2,(1,m2))/tmat
hmu2 = np.sum(rprob*mul,0)
mul -= hmu2
hrho2 = np.sum(rprob*(rhol + mul**2),0)
logz2 = logz2.ravel()
hmu2 = hmu2.ravel()
hrho2 = hrho2.ravel()

# Comparison
# For each of logz, mu_hat, rho_hat, we always show the 3 cases with largest
# relative difference
mat1 = np.vstack((logz, hmu, hrho))
mat2 = np.vstack((logz2, hmu2, hrho2))
names = ('logz', 'mu_hat', 'rho_hat')
print 'Gaussian mixture potential:'
print 'p_l =', gm_p
print 'v_l =', gm_v
for k in range(3):
    v1 = mat1[k]; v2 = mat2[k]
    rdf = reldiff(v1,v2)
    ind = np.argsort(rdf)
    print '%s:' % names[k]
    for j in ind[-1:-4:-1]:
        print ('  rdf=%.4e (v1=%f,v2=%f):' +
               ' j=%d,cmu=%f,crho=%f') % (rdf[j],v1[j],v2[j],j,cmu2[j],
                                          crho2[j])
