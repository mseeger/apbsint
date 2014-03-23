#! /usr/bin/env python

# EPTOOLS Python Interface
# Example: Binary classification, coupled or factorized posterior
# mode.
# Can be configured with Gaussian or Laplace prior.
#
# Very similar to Matlab
#   testcode eptools/matlab/test/binclass/eptest_binclass.m
# Same data as glm-ie_v1.5/doc/classify_inference.m:
# Adult (a9a): Dimension 123, first 30000 for testing, last 2561
# for training.

import numpy as np
import scipy.sparse as ssp
import scipy.io  # DEBUG: Load Matlab files
import numbers
import time  # Profiling

import apbsint as abt

# Helper functions

# Main code

# Load dataset: CSV format. Inputs: Each row contains indices of nonzero
# features (value: 1).
# NOTE: Building a dense nested list is faster (by factor 10!) than
# building a ssp.lil_matrix row by row.
# NOTE: We remove 3 attributes (45, 116, 122), these are not touched by
# any patterns in the training set.
num_feat = n = 120  # After removing 3
tmat = []
fid = open('adult_a9a_inputs_comp.csv','r')
for line in fid:
    ind = [int(x) for x in line.split(',')]
    v = np.zeros(n+3,dtype=np.float64)
    v[ind] = 1.
    # Remove attributes 45, 116, 122
    tmat.append(list(np.hstack((v[:45], v[46:116], v[117:122]))))
fid.close()
num_cases = len(tmat)
print 'Dataset: Read %d cases.' % num_cases
inp_all = ssp.csr_matrix(tmat)
del tmat
num_test = 30000
num_train = num_cases-num_test
fid = open('adult_a9a_targets.csv','r')
targ_all = np.array([float(x) for x in fid.readline().split(',')],
                    dtype=np.float64)
fid.close()
if targ_all.size != num_cases:
    raise IndexError('Internal error: Wrong file size')

# Setup
#imode = 'CoupParallel'   # Coupled, parallel updating
#imode = 'CoupSequential' # Coupled, sequential updating
imode = 'Factorized'     # Factorized, sequential updating
is_fact = (imode == 'Factorized')
do_laplace = True        # Laplace prior
tau_lapl = 2./5.
#do_laplace = False       # Gaussian prior
#ssq_gauss = 25./4.
do_seldamp = True
#do_seldamp = False
seldamp_numk = 5
# DEBUG: Compare intermediate results of Matlab and Python code directly.
# We load Matlab variables from files written after init. (%d = 0) and after
# each sweep (%d = 1,2,...). Anything non-deterministic also loaded from there.
# NOTE: Implemented for 'Factorized' mode only.
#debug_matcomp_fname = '/home/seeger/exp/ept/eptb13_vars_it%d.mat'
try:
    debug_matcomp_fname
    do_debug_matcomp = True
except NameError:
    do_debug_matcomp = False

# Coupling factors
if not is_fact:
    bfct_test = abt.MatSparse(inp_all[:num_test,:].copy())
    bfct_train = abt.MatContainer([abt.MatEye(n),
                                   abt.MatSparse(inp_all[num_test:,:].copy())])
else:
    bfct_test = abt.MatFactorizedInf(inp_all[:num_test,:].copy())
    mx_tmp = ssp.vstack([ssp.eye(n,format='csr'), inp_all[num_test:,:]],
                        format='csr')
    bfct_train = abt.MatFactorizedInf(mx_tmp)
m = bfct_train.shape(0)
# Potential managers
if not do_laplace:
    pm_elem1 = abt.ElemPotManager('Gaussian',n,(0., ssq_gauss))
else:
    pm_elem1 = abt.ElemPotManager('Laplace',n,(0., tau_lapl))
pm_elem2 = abt.ElemPotManager('Probit',num_train,
                              (targ_all[num_test:].copy(), 0.))
pman_train = abt.PotManager((pm_elem1, pm_elem2))
pman_test = abt.PotManager(abt.ElemPotManager('Probit',num_test,
                                              (targ_all[:num_test].copy(), 0.)))
print 'Prior:      %s\nLikelihood: %s' % (pman_train.elem[0].name,
                                          pman_train.elem[1].name)
# Initialization of representation parameters
if not is_fact:
    model_train = abt.ModelCoupled(bfct_train,pman_train)
    model_test = abt.ModelCoupled(bfct_test,pman_test)
    repres = abt.RepresentationCoupled(bfct_train,
                                       keep_margs=(imode == 'CoupParallel'))
    if imode == 'CoupParallel':
        inf_driv = abt.EPCoupParallelInfDriver(model_train,repres)
    else:
        inf_driv = abt.EPCoupSequentialInfDriver(model_train,repres)
    if not do_laplace:
        inf_driv.init('ADF')
    else:
        # Same initialization as in glm-ie
        tvec = np.empty(m)
        tvec[:] = 1.
        repres.setpi(tvec)
        tvec[:n] = 0.
        tvec[n:] = 0.5*targ_all[num_test:]
        repres.setbeta(tvec)
        repres.refresh()
else:
    model_train = abt.ModelFactorized(bfct_train,pman_train)
    model_test = abt.ModelFactorized(bfct_test,pman_test)
    repres = abt.RepresentationFactorized(bfct_train)
    inf_driv = abt.EPFactorizedInfDriver(model_train,repres)
    if not do_laplace:
        # 'cav_var' argument does not matter, because coupling matrix for
        # Gaussian potentials is I
        inf_driv.init('ADF')
    else:
        # Below, we skip the prior potentials in the 1st sweep. We set their
        # pi's to 1, beta's to 0 (as in coupled mode and glm-ie)
        tvec = np.zeros(repres.size_pars())
        repres.setbeta(tvec)
        tvec[:n] = 1.
        repres.setpi(tvec)
        repres.refresh()
    if do_seldamp:
        if not do_laplace:
            # Exclude Gaussian potentials from selective damping mechanism
            # (they are not updated on anyway)
            sd_subind = pman_train.filterpots(set(['Gaussian']))
            repres.seldamp_reset(seldamp_numk,sd_subind,True)
        else:
            repres.seldamp_reset(seldamp_numk)
    # DEBUG
    if do_debug_matcomp:
        deb_mc = scipy.io.loadmat(debug_matcomp_fname % 0)
        deb_ep_pi = deb_mc['ep_pi'].ravel()
        deb_ep_beta = deb_mc['ep_beta'].ravel()
        deb_marg_pi = deb_mc['marg_pi'].ravel()
        deb_marg_beta = deb_mc['marg_beta'].ravel()
        print ('DEBUG[%d]: df(ep_pi)=%.4e, df(ep_beta)=%.4e, ' +
               'df(m_ep)=%.4e, df(m_beta)=%.4e') % \
               (0, abt.helpers.maxreldiff(repres.ep_pi,deb_ep_pi),
                abt.helpers.maxreldiff(repres.ep_beta,deb_ep_beta),
                abt.helpers.maxreldiff(repres.marg_pi,deb_marg_pi),
                abt.helpers.maxreldiff(repres.marg_beta,deb_marg_beta))
# Initial test set prediction
opts = abt.helpers.Struct()
opts.imode = imode
opts.ptype = 3
(h_q, rho_q, logz, h_p, rho_p) = inf_driv.predict(model_test,opts)
acc = 100.*float((np.sign(h_q)==targ_all[:num_test]).sum())/num_test
loglh = logz.sum()/num_test
print ('\nInitial predictions:\nAccuracy: %4.2f%%\n'
       'Test set log likelihood: %.6f') % (acc, loglh)

# EP inference
print '\nRunning EP inference (mode: %s)...' % imode
if imode == 'CoupParallel':
    opts = abt.helpers.Struct()
    opts.imode = imode
    opts.maxit = 20       # Max. number sweeps
    #opts.maxit = 40       # Max. number sweeps
    opts.deltaeps = 1e-4  # Convergence threshold
    opts.damp = 0.        # No damping
    opts.caveps = 1e-5
    opts.verbose = 1
    opts.res_det = True
    opts.bc_testmodel = model_test
    t_start = time.time()
    (res, res_det) = inf_driv.inference(opts)
    t_stop = time.time()
    print 'Time(inference): %.6fs' % (t_stop-t_start)
    if res.rstat==1:
        print '\nDone MAXIT iterations.'
    else:
        print '\nConverged to DELTAEPS accuracy.'
    print ('Number sweeps:   %d\n' +
           'Final delta:     %f\n' +
           'Number of skips: %d') % (res.nit,res.delta,res.nskip)
elif imode == 'CoupSequential':
    opts = abt.helpers.Struct()
    opts.imode = imode
    opts.maxit = 15       # Max. number sweeps
    opts.deltaeps = 1e-4  # Convergence threshold
    opts.damp = 0.        # No damping
    opts.caveps = 1e-5
    opts.skipeps = 1e-7
    opts.refresh = True   # Recompute repres. after each sweep
    opts.verbose = 1
    opts.res_det = True
    #opts.bc_testmodel = model_test
    t_start = time.time()
    (res, res_det) = inf_driv.inference(opts)
    t_stop = time.time()
    print 'Time(inference): %.6fs' % (t_stop-t_start)
    if res.rstat==1:
        print '\nDone MAXIT iterations.'
    else:
        print '\nConverged to DELTAEPS accuracy.'
    print ('Number sweeps:   %d\n' +
           'Final delta:     %f') % (res.nit,res.delta)
    print 'Skip histogram: ', res.nskip
elif imode == 'Factorized':
    opts = abt.helpers.Struct()
    opts.imode = imode
    opts.maxit = 100      # Max. number sweeps
    opts.deltaeps = 1e-4  # Convergence threshold
    opts.damp = 0.        # No damping
    opts.piminthres = 1e-7
    opts.refresh = True   # Recompute repres. after each sweep
    opts.verbose = 1
    opts.res_det = True
    if not do_laplace:
        # Gaussian potentials have coupling factor I, so EP updates do
        # not change anything. Skip them
        opts.skip_gauss = True
    else:
        # Skip updates for Laplace prior potentials in 1st sweep
        opts.upd_1stsweep = set(['Probit'])
    if do_debug_matcomp:
        opts.deb_matcomp_fname = debug_matcomp_fname
    opts.bc_testmodel = model_test
    t_start = time.time()
    (res, res_det) = inf_driv.inference(opts)
    t_stop = time.time()
    print 'Time(inference): %.6fs' % (t_stop-t_start)
    if res.rstat==1:
        print '\nDone MAXIT iterations.'
    else:
        print '\nConverged to DELTAEPS accuracy.'
    print ('Number sweeps:   %d\n' +
           'Final delta:     %f') % (res.nit,res.delta)
    print 'Skip histogram:  ', res.nskip
    if do_seldamp:
        print 'Select. damps:   %d' % res.nsdamp
else:
    raise ValueError("IMODE: Unknown value '" + imode + "'")

# Predictions
opts = abt.helpers.Struct()
opts.imode = imode
#opts.ptype = 0
opts.ptype = 3
#h_q = inf_driv.predict(model_test,opts)
(h_q, rho_q, logz, h_p, rho_p) = inf_driv.predict(model_test,opts)
acc = 100.*float((np.sign(h_q)==targ_all[:num_test]).sum())/num_test
loglh = logz.sum()/num_test
print ('\nFinal predictions:\nAccuracy: %4.2f%%\n'
       'Test set log likelihood: %.6f') % (acc, loglh)
