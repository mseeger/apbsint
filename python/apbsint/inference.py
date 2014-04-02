"""
inference
=========

EP inference drivers (in the future, other inference methods could be
implemented here as well).

"""

import numpy as np
import scipy.linalg as sla
import scipy.sparse as ssp
import scipy.io  # DEBUG: Load Matlab files
import numbers
import time  # For profiling

import apbsint.helpers as helpers
import apbsint.coup_fact as cf
import apbsint.utilities as ut
import apbsint.eptools_ext as epx

__all__ = ['InfDriver', 'CoupledInfDriver', 'EPCoupParallelInfDriver',
           'EPCoupSequentialInfDriver', 'EPFactorizedInfDriver']

# Inference driver classes

class InfDriver:
    """
    InfDriver
    =========

    Base class for approximate inference drivers.
    A driver is configured by a Model and a Representation. The former
    contains model descriptions (potential manager, coupling factor), the
    latter contains the representation of the posterior approximation.

    The 'inference' method operates on the representation, fitting the
    approximate to the true posterior. The 'predict' method computes
    predictive moments on a test Model.

    If the model contains bivariate precision potentials, the Representation
    must be associated with the PotManager (see Representation constructor).
    """
    def __init__(self,model,rep):
        if not isinstance(model,ut.Model):
            raise TypeError('MODEL must be instance of apbsint.Model')
        if not isinstance(rep,ut.Representation):
            raise TypeError('REP must be instance of apbsint.Representation')
        if not model.bfact is rep.bfact:
            raise ValueError('MODEL.BFACT, REP.BFACT must be same object')
        pman = model.potman
        if pman.num_bvprec > 0:
            if not (helpers.check_vecsize(rep.ep_taua,pman.num_bvprec) and
                    rep.tauind is pman.tauind):
                raise ValueError('Bivariate precision potentials: REP must be associated to MODEL.POTMAN')
        elif rep.ep_taua is not None:
            raise ValueError('Bivariate precision potentials? REP yes, MODEL.POTMAN no')
        self.model = model
        self.rep = rep

    def inference(self,opts):
        raise NotImplementedError('INFERENCE must be implemented')

    def predict(self,pmodel,opts):
        raise NotImplementedError('PREDICT must be implemented')

    def _predict_epcomp(self,pmodel,h_q,rho_q,a_q=None,c_q=None):
        """
        Helper for 'predict'. Given Gaussian moments in 'h_q', 'rho_q', runs
        local EP computations and returns 'logz', 'h_p', 'rho_p'. See
        docstring of 'CoupledInfDriver.predict'.

        If 'pmodel' contains bivariate precision potentials, the Gamma
        parameters (a, c) for the corresponding marginals have to be passed
        in 'a_q', 'c_q'. In this case, we return
           (logz, h_p, rho_p, a_p, c_p),
        where a_p, c_p are Gamma parameters of the predictive marginals.
        Note that mean = a/c, variance = a/c^2 for Gamma(a,c).
        """
        pbfact = pmodel.bfact
        (pm, n) = pbfact.shape()
        ppotman = pmodel.potman
        ppotman.check_internal()
        rstat = np.empty(pm,dtype=np.int32)
        alpha = np.empty(pm)
        nu = np.empty(pm)
        logz = np.empty(pm)
        has_bvp = (ppotman.num_bvprec>0)
        if not has_bvp:
            epx.epupdate_parallel(ppotman.potids,ppotman.numpot,ppotman.parvec,
                                  ppotman.parshrd,ppotman.annobj,h_q,rho_q,
                                  rstat,alpha,nu,logz)
        else:
            pmb = ppotman.num_bvprec
            if not (helpers.check_vecsize(a_q,pmb) or
                    helpers.check_vecsize(c_q,pmb)):
                raise TypeError('Bivariate precision potentials: A_Q, C_Q must be vectors of size {0}'.format(pmb))
            a_p = np.empty(pmb)
            c_p = np.empty(pmb)
            epx.epupdate_parallel_bvprec(ppotman.potids,ppotman.numpot,
                                         ppotman.parvec,ppotman.parshrd,
                                         ppotman.annobj,h_q,rho_q,a_q,c_q,
                                         rstat,alpha,nu,a_p,c_p,logz)
        indok = np.nonzero(rstat)[0]
        tvec = 1. - nu[indok]*rho_q[indok]
        indok2 = np.nonzero(tvec >= 1e-9)[0]
        if indok2.shape[0] < pm:
            indok = indok[indok2].copy()
            # ATTENTION: This is **very** slow:
            indnok = [x for x in range(pm) if x not in set(indok)]
            logz[indnok] = 0.
            h_p = h_q.copy()
            rho_p = rho_q.copy()
            h_p[indok] += alpha[indok]*rho_q[indok]
            rho_p[indok] *= tvec[indok2]
        else:
            h_p = h_q + alpha*rho_q
            rho_p = rho_q*tvec
        if not has_bvp:
            return (logz, h_p, rho_p)
        else:
            return (logz, h_p, rho_p, a_p, c_p)

    def _infer_check_commonargs(self,opts):
        """
        Checks common arguments of 'inference' implementations in
        subclasses and assigns default values.
        'opts.imode' must be 'CoupParallel', 'CoupSequential' or
        'Factorized'.
        """
        if not (opts.imode == 'CoupParallel' or
                opts.imode == 'CoupSequential' or
                opts.imode == 'Factorized'):
            raise ValueError('OPTS.IMODE has wrong value')
        if not (isinstance(opts.maxit,numbers.Integral) and opts.maxit>=1):
            raise TypeError('OPTS.MAXIT wrong')
        if not (isinstance(opts.deltaeps,numbers.Real) and opts.deltaeps>0.):
            raise TypeError('OPTS.DELTAEPS wrong')
        try:
            if not (isinstance(opts.damp,numbers.Real) and opts.damp>=0. and
                    opts.damp<1.):
                raise TypeError('OPTS.DAMP wrong')
        except AttributeError:
            opts.damp = 0.
        try:
            if not isinstance(opts.res_det,bool):
                raise TypeError('OPTS.RES_DET wrong')
        except AttributeError:
            opts.res_det = False
        try:
            if not isinstance(opts.verbose,numbers.Integral):
                raise TypeError('OPTS.VERBOSE wrong')
        except AttributeError:
            opts.verbose = 0
        if opts.imode != 'Factorized':
            try:
                if not isinstance(opts.bc_testmodel,ut.ModelCoupled):
                    raise TypeError('OPTS.BC_TESTMODEL must be apbsint.ModelCoupled')
            except AttributeError:
                pass
            try:
                if not (isinstance(opts.caveps,numbers.Real) and opts.caveps>0.):
                    raise TypeError('OPTS.CAVEPS wrong')
            except AttributeError:
                opts.caveps = 1e-5
        else:
            try:
                if not isinstance(opts.bc_testmodel,ut.ModelFactorized):
                    raise TypeError('OPTS.BC_TESTMODEL must be apbsint.ModelFactorized')
            except AttributeError:
                pass
        if opts.imode != 'CoupParallel':
            try:
                if not isinstance(opts.refresh,bool):
                    raise TypeError('OPTS.REFRESH wrong')
            except AttributeError:
                opts.refresh = True

    def _binclass_print_teststats(self,pmodel,targets,imode):
        """
        Helper for 'inference'. Only for binary classification right now.
        'targets' must be target vector (values -1, +1).
        Calls 'predict', computes test set statistics (accuracy, avg. log
        likelihood) and prints information.
        """
        popts = helpers.Struct()
        popts.imode = imode
        popts.ptype = 3
        (h_q, rho_q, logz, h_p, rho_p) = self.predict(pmodel,popts)
        nte = targets.shape[0]
        acc = 100.*float((np.sign(h_q)==targets).sum())/nte
        loglh = logz.sum()/nte
        print ('Test set predictions: Accuracy: %4.2f%%, '
               'log likelihood: %.6f') % (acc, loglh)

    def _binclass_assemble_targets(self,opts):
        # Assemble test target vector
        pman = opts.bc_testmodel.potman
        targets = np.empty(pman.size)
        off = 0
        for el in pman.elem:
            sz = el.size
            targets[off:off+sz] = el.pars[0]
            off += sz
        return targets

class CoupledInfDriver(InfDriver):
    """
    CoupledInfDriver
    ================

    Base class of coupled mode inference drivers. 'predict' is implemented
    here, as well as 'init'.

    """
    def __init__(self,model,rep):
        if not isinstance(model,ut.ModelCoupled):
            raise TypeError('MODEL must be instance of apbsint.ModelCoupled')
        if not isinstance(rep,ut.RepresentationCoupled):
            raise TypeError('REP must be instance of apbsint.RepresentationCoupled')
        InfDriver.__init__(self,model,rep)

    def predict(self,pmodel,opts):
        """
        Prediction on test model 'pmodel' (type apbsint.ModelCoupled). 'opts'
        is a struct with attributes:
        - imode: Inference mode ('CoupParallel', 'CoupSequential')
        - ptype: What predictive moments are returned?
          0: Gaussian means h_q
          1: Gaussian moments (h_q, rho_q)
          2: Predictive moments (logz, h_p, rho_p)
          3: Everything (h_q, rho_q, logz, h_p, rho_p)
          Here, the predictive marginal at a test point is
            p(s) = Z^-1 t(s) q(s),
          where t(s) is the potential, q(s) the Gaussian marginal. 'logz'
          returns the log Z values

        NOTE: We assume that the posterior covariance A^-1 is in 'rep.post_cov'
        if 'opts.imode'=='CoupParallel'. In the other modes, A^-1 is
        recomputed.

        If training model 'model' has bivar. precision potentials, 'pmodel'
        may as well (but does not have to). If so, the dimensionality of tau
        must be the same. In this case, further 'ptype' values are:
         4: (h_q, a_q, c_q)
         5: (h_q, rho_q, a_q, c_q)
         6: (logz, h_p, rho_p, a_p, c_p)
         7: (h_q, rho_q, a_q, c_q, logz, h_p, rho_p, a_p, c_p)
        Here, (a_q, c_q) and (a_p, c_p) are Gamma parameters of q(tau) and
        predictive p_hat(tau) respectively.
        """
        model = self.model
        rep = self.rep
        if not isinstance(pmodel,ut.ModelCoupled):
            raise TypeError('PMODEL must be instance of apbsint.ModelCoupled')
        pbfact = pmodel.bfact
        (pm, n) = pbfact.shape()
        if n != model.bfact.shape(1):
            raise TypeError('PMODEL, MODEL: Different number of variables')
        use_cov = (opts.imode == 'CoupParallel')
        if not (isinstance(opts.ptype,numbers.Integral) and opts.ptype>=0 and
                opts.ptype<=3):
            raise ValueError('opts.ptype wrong')
        nbvp = pmodel.potman.num_bvprec
        if nbvp>0:
            if (model.potman.num_bvprec==0 or
                model.potman.num_tau != pmodel.potman.num_tau):
                raise TypeError('Bivariate precision potentials: PMODEL, MODEL not consistent')
        # HIER!
        # Compute Gaussian moments
        h_q = np.empty(pm)
        rho_q = np.empty(pm) if opts.ptype>0 else None
        rep.predict(pbfact,h_q,rho_q,use_cov)
        if opts.ptype==0:
            return h_q
        elif opts.ptype==1:
            return (h_q, rho_q)
        if opts.ptype==2:
            res = ()
        else:
            res = (h_q, rho_q)
        return res + self._predict_epcomp(pmodel,h_q,rho_q)

    def init(self,mode,refresh=True):
        """
        Initialize EP parameters according to mode 'mode'. The representation
        is refreshed afterwards iff 'refresh'==True. Modes:
        - 'ADF': Parameters for all non-Gaussian potentials set to zero. For
          Gaussian potentials, parameters are set to represent them (they do
          not change afterwards).
        """
        if mode.upper() == 'ADF':
            bfact = self.model.bfact
            potman = self.model.potman
            rep = self.rep
            m, n = bfact.shape()
            potman.check_internal()
            ep_pi = np.zeros(rep.size_pars())
            ep_beta = np.zeros(rep.size_pars())
            off = 0
            for el in potman.elem:
                numk = el.size
                if el.name == 'Gaussian':
                    ep_pi[off:off+numk] = 1./el.pars[1]
                    ep_beta[off:off+numk] = el.pars[0]/el.pars[1]
                off += numk
            rep.setpi(ep_pi)
            rep.setbeta(ep_beta)
            if refresh:
                rep.refresh()
        else:
            raise ValueError("Unknown mode '" + mode + "'")

class EPCoupParallelInfDriver(CoupledInfDriver):
    """
    EPCoupParallelInfDriver
    =======================

    Implements parallel updating expectation propagation in 'inference'.

    """
    def __init__(self,model,rep):
        CoupledInfDriver.__init__(self,model,rep)

    def inference(self,opts):
        """
        Update representation by running parallel updating EP. One sweep
        consists of (a) parallel EP updates on all non-Gaussian potentials
        (model.potman.updind), then (b) a recomputation (refresh) of the
        representation. The latter stores the posterior covariance in
        rep.post_cov, which is recycled by 'predict'. 'opts' attributes:
        - maxit: Maximum number of sweeps
        - deltaeps: Threshold for convergence (statistic based on relative
          change of Gaussian means and stddevs.)
        - damp: Damping constant (def.: 0 -> no damping)
        - caveps: Update on k is skipped if
            cavvar_k / margvar_k > 1/caveps,
          or if there is some numerical failure
        - res_det: Return detailed results in 'res_det' (below)? Def.: False
        - verbose: Verbosity level (0: no messages, 1: some messages). Def.: 0
        - bc_testmodel: Optional. Only for binary classification right now.
          Test set model (type apbsint.ModelCoupled). Test set accuracy and
          avg. log likelihood are computed and printed after each sweep.
        Returns 'res' or '(res, res_det)' (latter if 'opts.res_det'==True).
        'res' attributes:
        - rstat: Return status (0: Converged to 'deltaeps'; 1: Done
          'maxit' sweeps)
        - nit: Number of sweeps done
        - delta: Value convergence statistic after last sweep
        - nskip: Total number of skipped updates across all sweeps
        'res_det' attributes (optional):
        - delta: Value after each sweep
        - nskip: Value after each sweep
        """
        #t_start0=time.time()
        if not self.rep.keep_margs:
            raise ValueError('REP.KEEP_MARGS must be True')
        opts.imode = 'CoupParallel'
        self._infer_check_commonargs(opts)
        # Initialization
        res = helpers.Struct()
        res.rstat = 1
        res.nskip = 0
        if opts.res_det:
            res_det = helpers.Struct()
            res_det.delta = []
            res_det.nskip = []
        bfact = self.model.bfact
        potman = self.model.potman
        rep = self.rep
        m, n = bfact.shape()
        potman.check_internal()
        try:
            targets = self._binclass_assemble_targets(opts)
            do_teststats = True
        except AttributeError:
            do_teststats = False
        #t_stop=time.time()
        #print 'Time(inference::init): %.8fs' % (t_stop-t_start0)
        # Loop over sweeps
        mm = potman.updind.shape[0]
        cmu = np.empty(mm)
        crho = np.empty(mm)
        alpha0 = np.empty(mm)
        nu0 = np.empty(mm)
        rstat = np.empty(mm,dtype=np.int32)
        sz = rep.ep_pi.shape[0]
        new_pi = np.empty(sz)
        new_beta = np.empty(sz)
        old_margs = np.empty(2*mm)
        new_margs = np.empty(2*mm)
        for res.nit in range(1,opts.maxit+1):
            #t_start1=time.time()
            # Local EP updates
            # We update only on potentials in 'potman.updind' (excludes
            # Gaussians)
            indok = potman.updind
            # Compute cavity marginals
            #t_start=time.time()
            cmu[:] = rep.marg_means[indok]
            crho[:] = rep.marg_vars[indok]
            tvec = 1. - rep.ep_pi[indok]*crho
            indok2 = np.nonzero(tvec >= opts.caveps)[0]
            if indok2.shape[0] == mm:
                cmu -= crho*rep.ep_beta[indok]
                cmu /= tvec
                crho /= tvec
            else:
                indok = indok[indok2].copy()
                cmu[indok2] -= crho[indok2]*rep.ep_beta[indok]
                tvec = 1./tvec[indok2]
                cmu[indok2] *= tvec
                crho[indok2] *= tvec
            #t_stop=time.time()
            #print 'Time(inference:comp_cav): %.8fs' % (t_stop-t_start)
            #t_start=time.time()
            epx.epupdate_parallel(potman.potids,potman.numpot,potman.parvec,
                                  potman.parshrd,potman.annobj,cmu,crho,rstat,
                                  alpha0,nu0,None,potman.updind)
            #t_stop=time.time()
            #print 'Time(epupdate_parallel): %.8fs' % (t_stop-t_start)
            # Update EP parameters, and figure out where skips happened
            #t_start=time.time()
            if indok2.shape[0] < mm:
                # Complement of 'indok2'
                tarr = np.ones(mm,dtype=np.bool)
                tarr[indok2] = False
                indnok = np.nonzero(tarr)[0]
                rstat[indnok] = 0  # Filter out undef. cavity positions
                #print 'Time(indnok stuff): %.8fs' % (time.time()-t_start)
            indok2 = np.nonzero(rstat)[0]
            new_pi[:] = rep.ep_pi
            new_beta[:] = rep.ep_beta
            # 'nu0', 'alpha0' must remain full size
            if indok2.shape[0] < mm:
                nu = nu0[indok2].copy()
                alpha = alpha0[indok2].copy()
            else:
                nu = nu0
                alpha = alpha0
            # Just a sanity check (this should not fire)
            tvec = 1. - nu*crho[indok2]
            indok3 = np.nonzero(tvec >= 1e-7)[0]
            if indok3.shape[0] < indok2.shape[0]:
                print 'UUPS[EPCoupParallelInfDriver.inference]: On %d' % (indok2.shape[0]-indok3.shape[0])
                tvec = tvec[indok3].copy()
                nu = nu[indok3].copy()
                alpha = alpha[indok3].copy()
                indok2 = indok2[indok3].copy()
            if indok2.shape[0] < mm:
                indok = potman.updind[indok2].copy()
            else:
                indok = potman.updind
            new_pi[indok] = nu/tvec
            new_beta[indok] = (cmu[indok2]*nu + alpha)/tvec
            #t_stop=time.time()
            #print 'Time(inference:updpars1): %.8fs' % (t_stop-t_start)
            #t_start=time.time()
            # Damping
            if opts.damp > 0.:
                new_pi[indok] = (1.-opts.damp)*new_pi[indok] + \
                                opts.damp*rep.ep_pi[indok]
                new_beta[indok] = (1.-opts.damp)*new_beta[indok] + \
                                  opts.damp*rep.ep_beta[indok]
            nskip = mm - indok.shape[0]    # Number of skips
            # Recompute representation (refresh)
            # Posterior covariance is kept in 'rep.post_cov'
            indok = potman.updind
            old_margs[:mm] = rep.marg_means[indok]
            old_margs[mm:] = np.sqrt(rep.marg_vars[indok])
            rep.setpi(new_pi)
            rep.setbeta(new_beta)
            #t_stop=time.time()
            #print 'Time(inference:updpars2): %.8fs' % (t_stop-t_start)
            rep.refresh()
            new_margs[:mm] = rep.marg_means[indok]
            new_margs[mm:] = np.sqrt(rep.marg_vars[indok])
            res.delta = helpers.maxreldiff(old_margs,new_margs)
            # End of sweep: Write results
            res.nskip += nskip
            if opts.res_det:
                res_det.delta.append(res.delta)
                res_det.nskip.append(nskip)
            if opts.verbose>0:
                print 'It. %d: delta=%f, nskip=%d' % (res.nit,res.delta,nskip)
            if do_teststats:
                self._binclass_print_teststats(opts.bc_testmodel,targets,
                                               opts.imode)
            # Convergence?
            if res.delta < opts.deltaeps:
                res.rstat = 0
                break
            #t_stop1=time.time()
            #print 'Time(inference::sweep): %.8fs' % (t_stop1-t_start1)
        # Timing
        #t_stop0=time.time()
        #print 'Time(inference(ALL)): %.8fs' % (t_stop0-t_start0)
        # Return stuff
        if opts.res_det:
            return (res, res_det)
        else:
            return res

    def predict(self,pmodel,opts):
        # Make sure that 'opts.imode' is correct
        opts.imode = 'CoupParallel'
        return CoupledInfDriver.predict(self,pmodel,opts)

class EPCoupSequentialInfDriver(CoupledInfDriver):
    """
    EPCoupSequentialInfDriver
    =========================

    Implements sequential updating expectation propagation in 'inference'.
    The representation is updated after each local EP update, using a
    Cholesky update/downdate. This is much slower than parallel updating
    in general, but may converge more reliably.
    If 'rep.keep_margs'==True, the marginals are kept up-2-date at all
    times. Right now, this is a waste of time.
    TODO: Implement optimized update scheduling, based on forward scoring
    and marginal moments.

    """
    def __init__(self,model,rep):
        CoupledInfDriver.__init__(self,model,rep)

    def inference(self,opts):
        """
        Update representation by running sequential updating EP. In a sweep,
        we iterate over all potentials in model.potman.updind (non-Gaussians)
        in random ordering. 'opts' attributes:
        - maxit: Maximum number of sweeps
        - deltaeps: Threshold for convergence (statistic based on relative
          change of Gaussian means and stddevs.)
        - damp: Damping constant (def.: 0 -> no damping)
          On top of this, we apply selective damping to make sure that
            1 + (Delta pi_k) margvar_k >= caveps
        - caveps: Update on k is skipped if
            cavvar_k / margvar_k > 1/caveps,
          or if there is some numerical failure
        - skipeps: Update is skipped if absolute change in pi_k is smaller
          than 'skipeps'
        - refresh: Refresh representation after each sweep? Def.: True
        - upd_1stsweep: Optional. Set of str. If given, in the 1st sweep, we
          only update on potentials whose type name is contained in the set.
        - res_det: Return detailed results in 'res_det' (below)? Def.: False
        - verbose: Verbosity level (0: no messages, 1: some messages). Def.: 0
        - bc_testmodel: Optional. See EPCoupParallelInfDriver.inference.
        Returns 'res' or '(res, res_det)' (latter if 'opts.res_det'==True).
        Each update results in a skip status, summarized in 'nskip'
        histograms:
        - 0: Not skipped
        - 1: Skipped due to cavity marginal ('caveps') or local EP failure
        - 2: Skipped due to small change ('skipeps', selective damping
        - 3: Skipped due to Cholesky up/downdate error
        'res' has attributes:
        - rstat: Return status (0: Converged to 'deltaeps'; 1: Done
          'maxit' sweeps)
        - nit: Number of sweeps done
        - delta: Value convergence statistic after last sweep
        - nskip: Skip status histogram (vector of size 4), summed over all
          updates and sweeps
        'res_det' attributes (optional):
        - delta: Value after each sweep
        - nskip: Matrix, each row skip status histogram for a sweep
        """
        opts.imode = 'CoupSequential'
        self._infer_check_commonargs(opts)
        try:
            if not (isinstance(opts.skipeps,numbers.Real) and
                    opts.skipeps>0.):
                raise TypeError('OPTS.SKIPEPS wrong')
        except AttributeError:
            opts.skipeps = 1e-8
        try:
            if not isinstance(opts.upd_1stsweep,set):
                raise TypeError('OPTS.UPD_1STSWEEP wrong')
            do_1stsweep = True
        except AttributeError:
            do_1stsweep = False
        # Initialization
        res = helpers.Struct()
        res.rstat = 1
        res.nskip = np.zeros(4,dtype=np.int32)
        if opts.res_det:
            res_det = helpers.Struct()
            res_det.delta = []
            res_det.nskip = []
        bfact = self.model.bfact
        potman = self.model.potman
        rep = self.rep
        m, n = bfact.shape()
        potman.check_internal()
        if do_1stsweep:
            ind_swp1 = set(potman.filterpots(opts.upd_1stsweep))
        try:
            targets = self._binclass_assemble_targets(opts)
            do_teststats = True
        except AttributeError:
            do_teststats = False
        # Loop over sweeps
        vvec = np.empty(n)
        for res.nit in range(1,opts.maxit+1):
            updind = np.random.permutation(potman.updind)
            if do_1stsweep and res.nit==1:
                updind = [x for x in updind if x in ind_swp1]
            if len(updind)==0:
                raise IndexError('UPDIND empty: No potentials to update on?')
            # Loop over potentials in UPDIND
            nskip = [0]*4
            delta = 0.
            for j in updind:
                # np.int32 not instanceof numbers.Integral (sucks!)
                j = int(j)
                # Compute cavity marginals
                ep_pi = rep.ep_pi[j]
                ep_beta = rep.ep_beta[j]
                do_skip = 0
                # vvec = L^-1 B[j,:] required below
                (mu, rho) = rep.get_marg(j,vvec)
                tscal = 1. - ep_pi*rho
                if tscal >= opts.caveps:
                    crho = rho/tscal
                    cmu = (mu - ep_beta*rho)/tscal
                    # Local EP update
                    (rstat, alpha, nu, logz) \
                        = epx.epupdate_single_pman(potman.potids,potman.numpot,
                                                   potman.parvec,
                                                   potman.parshrd,
                                                   potman.annobj,j,cmu,crho)
                    if rstat == 0:
                        do_skip = 1  # Local EP update failed
                    else:
                        tscal = 1. - nu*crho;
                        if tscal>=1e-7:
                            new_pi = nu/tscal
                            new_beta = (cmu*nu + alpha)/tscal
                        else:
                            do_skip = 1  # Local EP update failed
                else:
                    do_skip = 1  # Cavity marginal invalid
                if do_skip == 0:
                    # Damping
                    dfl_pi = new_pi-ep_pi  # Full update
                    dfl_beta = new_beta-ep_beta
                    delpi = (1.-opts.damp)*dfl_pi
                    delbeta = (1.-opts.damp)*dfl_beta
                    delpi2 = delpi; delbeta2 = delbeta
                    # Selective damping
                    if delpi*rho + 1. < opts.caveps:
                        delpi = (opts.caveps-1.)/rho
                        delbeta = (delpi/dfl_pi)*delbeta
                    new_pi = ep_pi+delpi
                    new_beta = ep_beta+delbeta
                    if abs(delpi)>=opts.skipeps:
                        # Update representation
                        try:
                            rep.update_single(j,delpi,delbeta,vvec)
                        except sla.LinAlgError:
                            do_skip = 3  # Numerical error Cholesky up/down
                    else:
                        # Small |delpi| counted as skip only if due to
                        # selective damping
                        do_skip = 2 if abs(delpi2)>=opts.skipeps else 4
                nskip[do_skip if do_skip<4 else 0] += 1
                if do_skip == 0:
                    hrho = crho*(1. - nu*crho)
                    hmu = cmu + alpha*crho
                    delta = max(delta,
                                helpers.maxreldiff(np.array([hmu,
                                                             np.sqrt(hrho)]),
                                                   np.array([mu,
                                                             np.sqrt(rho)])))
            # Write back results
            res.nskip += np.array(nskip,dtype=np.int32)
            res.delta = delta
            if opts.res_det:
                res_det.delta.append(delta)
                res_det.nskip.append(nskip)
            if opts.refresh:
                rep.refresh()
            if opts.verbose>0:
                print 'It. %d: delta=%f, nnskip=%d' % (res.nit,res.delta,
                                                       sum(nskip[1:]))
                print '   nskip=', nskip
            if do_teststats:
                self._binclass_print_teststats(opts.bc_testmodel,targets,
                                               opts.imode)
            if res.delta < opts.deltaeps:
                res.rstat = 0
                break
        # Return stuff
        if opts.res_det:
            return (res, res_det)
        else:
            return res

class EPFactorizedInfDriver(InfDriver):
    """
    EPFactorizedInfDriver
    =====================

    Implements expectation propagation inference in factorized mode.

    """
    def __init__(self,model,rep):
        if not isinstance(model,ut.ModelFactorized):
            raise TypeError('MODEL must be instance of apbsint.ModelFactorized')
        if not isinstance(rep,ut.RepresentationFactorized):
            raise TypeError('REP must be instance of apbsint.RepresentationFactorized')
        InfDriver.__init__(self,model,rep)

    def init(self,mode,refresh=True,cav_var=1.):
        """
        Initialize EP parameters according to mode 'mode'. The representation
        is refreshed afterwards iff 'refresh'==True. Modes:
        - 'ADF': Parameters for all non-Gaussian potentials set to zero. For
          Gaussian potentials, we use a heuristic which depends on 'cav_var'
          (see technical report).
          NOTE: For a potential j with V_j = {i} and B[j,i] = 1, the Gaussian
          potential is represented exactly (independent of 'cav_var'), and the
          EP parameters remain fixed there.
        """
        if mode.upper() == 'ADF':
            bfact = self.model.bfact
            bmat = bfact.get_mat()
            potman = self.model.potman
            rep = self.rep
            m, n = bfact.shape()
            potman.check_internal()
            ep_pi = np.zeros(rep.size_pars())
            ep_beta = np.zeros(rep.size_pars())
            off = 0
            for el in potman.elem:
                numk = el.size
                if el.name == 'Gaussian':
                    # If potential is N(s | y_j,ssq_j) and cv=='cav_var':
                    #   pi_ji = b_ji^2 / ( (|V_j|-1) cv + ssq_j )
                    #   beta_ji = b_ji y_j / ( (|V_j|-1) cv + ssq_j )
                    # Offset into EP parameter vectors:
                    off2 = bmat[:off].getnnz() if off>0 else 0
                    mx_tmp = bfact.b2fact[off:off+numk].copy()
                    sz2 = mx_tmp.getnnz()
                    # Number of nonzeros per row minus 1:
                    vjsz = mx_tmp.indptr[1:] - mx_tmp.indptr[:-1] - 1
                    tvec = 1./(cav_var*vjsz + el.pars[1])
                    mx_dg = ssp.diags(tvec,0)
                    mx_tmp = mx_dg * mx_tmp
                    ep_pi[off2:off2+sz2] = mx_tmp.data
                    mx_tmp = bmat[off:off+numk].copy()
                    assert mx_tmp.getnnz() == sz2
                    # Some y_j's could be zero, which would change the sparsity
                    # pattern. Have to go a detour here
                    nzind = mx_tmp.nonzero()
                    tvec *= el.pars[0]
                    mx_dg = ssp.diags(tvec,0)
                    mx_tmp = mx_dg * mx_tmp
                    ep_beta[off2:off2+sz2] = mx_tmp[nzind[0],nzind[1]]
                off += numk
            rep.setpi(ep_pi)
            rep.setbeta(ep_beta)
            if refresh:
                rep.refresh()
        else:
            raise ValueError("Unknown mode '" + mode + "'")

    def predict(self,pmodel,opts):
        """
        Prediction on test model 'pmodel' (type apbsint.ModelFactorized).
        'opts' is a struct with attributes:
        - ptype: What predictive moments are returned?
          0: Gaussian means h_q
          1: Gaussian moments (h_q, rho_q)
          2: Predictive moments (logz, h_p, rho_p)
          3: Everything (h_q, rho_q, logz, h_p, rho_p)
          Here, the predictive marginal at a test point is
            p(s) = Z^-1 t(s) q(s),
          where t(s) is the potential, q(s) the Gaussian marginal. 'logz'
          returns the log Z values
        """
        model = self.model
        rep = self.rep
        if not isinstance(pmodel,ut.ModelFactorized):
            raise TypeError('PMODEL must be instance of apbsint.ModelFactorized')
        pbfact = pmodel.bfact
        pm, n = pbfact.shape()
        if n != model.bfact.shape(1):
            raise TypeError('PMODEL, MODEL: Different number of variables')
        if not (isinstance(opts.ptype,numbers.Integral) and opts.ptype>=0 and
                opts.ptype<=3):
            raise ValueError('opts.ptype wrong')
        # Compute Gaussian moments
        h_q = np.empty(pm)
        rho_q = np.empty(pm) if opts.ptype>0 else None
        rep.predict(pbfact,h_q,rho_q)
        if opts.ptype==0:
            return h_q
        elif opts.ptype==1:
            return (h_q, rho_q)
        if opts.ptype==2:
            res = ()
        else:
            res = (h_q, rho_q)
        return res + self._predict_epcomp(pmodel,h_q,rho_q)

    def inference(self,opts):
        """
        Update representation by running sweeps of EP updates. In one sweep,
        we iterate sequentially over all potentials in random ordering.
        Selective damping is used (and the SD representation updated) iff
        activated (see apbsint.RepresentationFactorized).
        'opts' attributes:
        - maxit: Maximum number of sweeps
        - deltaeps: Threshold for convergence (statistic based on relative
          change of Gaussian means and stddevs.)
        - damp: Damping constant (def.: 0 -> no damping)
        - piminthres: pi values for cavity moments must be > (.)/2 for an
          update not to fail. If selective damping is active, it aims to
          enforce that future cavity pi values are >= (.).
          Def.: 1e-8
        - refresh: If True, marginals are refreshed from messages after each
          sweep. Def.: True
        - skip_gauss: If True, EP updates are not done on potentials of type
          'Gaussian'. Def.: False
        - upd_1stsweep: See apbsint.EPCoupSequentialInfDriver.inference.
          Optional
        - res_det: Return detailed results in 'res_det' (below)? Def.: False
        - verbose: Verbosity level (0: no messages, 1: some messages). Def.: 0
        - bc_testmodel: See apbsint.EPCoupParallelInfDriver.inference.
          Optional
        Returns 'res' or '(res, res_det)' (latter if 'opts.res_det'==True).
        Each update results in a skip status, summarized in 'nskip'
        histograms:
        - 0: Not skipped
        - 1: Skipped due to invalid cavity marginal ('piminthres')
        - 2: Skipped due to local EP update error
        - 3: Skipped due to invalid new marginal ('piminthres')
        - 4: Skipped due to selective damping
        'res' attributes:
        - rstat: Return status (0: Converged to 'deltaeps'; 1: Done
          'maxit' sweeps)
        - nit: Number of sweeps done
        - delta: Value convergence statistic after last sweep
        - nskip: Skip status histogram (vector of size 5), summed over all
          updates and sweeps
        - nsdamp: Only if selective damping active. Number of non-skipped
          updates which were selectively damped
        'res_det' attributes (optional):
        - delta: Value after each sweep
        - nskip: Matrix, each row skip status histogram for a sweep
        - nsdamp: S.a. Value for each sweep
        """
        opts.imode = 'Factorized'
        self._infer_check_commonargs(opts)
        try:
            if not (isinstance(opts.piminthres,numbers.Real) and
                    opts.skipeps>0.):
                raise TypeError('OPTS.PIMINTHRES wrong')
        except AttributeError:
            opts.piminthres = 1e-8
        try:
            if not isinstance(opts.upd_1stsweep,set):
                raise TypeError('OPTS.UPD_1STSWEEP wrong')
            do_1stsweep = True
        except AttributeError:
            do_1stsweep = False
        try:
            if not isinstance(opts.skip_gauss,bool):
                raise TypeError('OPTS.SKIP_GAUSS wrong')
        except AttributeError:
            opts.skip_gauss = False
        # Initialization
        bfact = self.model.bfact
        potman = self.model.potman
        rep = self.rep
        m, n = bfact.shape()
        potman.check_internal()
        try:
            do_seldamp = (rep.sd_numk>0)
        except AttributeError:
            do_seldamp = False
        res = helpers.Struct()
        res.rstat = 1
        res.nskip = np.zeros(5,dtype=np.int32)
        if do_seldamp:
            res.nsdamp = 0
        if opts.res_det:
            res_det = helpers.Struct()
            res_det.delta = []
            res_det.nskip = []
            res_det.nsdamp = []
        if do_1stsweep:
            ind_swp1 = set(potman.filterpots(opts.upd_1stsweep))
        try:
            targets = self._binclass_assemble_targets(opts)
            do_teststats = True
        except AttributeError:
            do_teststats = False
        # DEBUG:
        # If 'opts.deb_matcomp_fname' is given, we directly compare against
        # intermediate results stored by Matlab. This is a file name string
        # with %d for 'res.nit' (sweep number; starting 1). The index
        # 'updind' is also loaded from that file
        try:
            if len(opts.deb_matcomp_fname)==0:
                raise ValueError('OPTS.DEB_MATCOMP_FNAME wrong')
            do_deb_matcomp = True
        except AttributeError:
            do_deb_matcomp = False
        # Loop over sweeps
        for res.nit in range(1,opts.maxit+1):
            if not do_deb_matcomp:
                if not opts.skip_gauss:
                    updind = np.int32(np.random.permutation(m))
                else:
                    updind = np.random.permutation(potman.updind)
                if do_1stsweep and res.nit==1:
                    # NOTE: This could be very slow...
                    updind = np.array([x for x in updind if x in ind_swp1],
                                      dtype=np.int32)
                    if updind.shape[0]==0:
                        raise IndexError('UPDIND empty: No potentials to update on?')
            else:
                deb_mc = scipy.io.loadmat(opts.deb_matcomp_fname % res.nit)
                updind = deb_mc['updind'].ravel()
            # Everything is done by epx.fact_sequpdates
            sz = updind.shape[0]
            rstat = np.empty(sz,dtype=np.int32)
            delta = np.empty(sz)
            if not do_seldamp:
                epx.fact_sequpdates(n,m,updind,potman.potids,potman.numpot,
                                    potman.parvec,potman.parshrd,potman.annobj,
                                    bfact.rowind,bfact.colind,bfact.bvals,
                                    rep.ep_pi,rep.ep_beta,rep.marg_pi,
                                    rep.marg_beta,opts.piminthres,opts.damp,
                                    rstat,delta)
            else:
                sd_dampfact = np.empty(sz)
                sd_nupd, sd_nrec = \
                    epx.fact_sequpdates(n,m,updind,potman.potids,potman.numpot,
                                        potman.parvec,potman.parshrd,
                                        potman.annobj,bfact.rowind,
                                        bfact.colind,bfact.bvals,rep.ep_pi,
                                        rep.ep_beta,rep.marg_pi,rep.marg_beta,
                                        opts.piminthres,opts.damp,rstat,delta,
                                        rep.sd_numvalid,rep.sd_topind,
                                        rep.sd_topval,rep.sd_subind,
                                        rep.sd_subexcl,sd_dampfact)
                # Among non-skipped updates, count those for which SD_DAMPFACT
                # larger than OPTS.DAMP
                nsdamp = np.sum(sd_dampfact[np.nonzero(rstat==0)] > opts.damp)
                res.nsdamp += nsdamp
                if opts.res_det:
                    res_det.nsdamp.append(nsdamp)
            res.delta = max(delta)
            nskip = [0]*5
            for k in xrange(5):
                nskip[k] = np.sum(rstat==k)
            res.nskip += np.array(nskip,dtype=np.int32)
            if opts.res_det:
                res_det.delta.append(res.delta)
                res_det.nskip.append(nskip)
            if opts.refresh:
                rep.refresh()
            if opts.verbose>0:
                print 'It. %d: delta=%f, nnskip=%d' % (res.nit,res.delta,
                                                       sum(nskip[1:]))
                if do_seldamp:
                    print '   nskip=', nskip, ', nsdamp=%d' % nsdamp
                else:
                    print '   nskip=', nskip
            if do_teststats:
                self._binclass_print_teststats(opts.bc_testmodel,targets,
                                               opts.imode)
            # DEBUG
            if do_deb_matcomp:
                deb_ep_pi = deb_mc['ep_pi'].ravel()
                deb_ep_beta = deb_mc['ep_beta'].ravel()
                deb_marg_pi = deb_mc['marg_pi'].ravel()
                deb_marg_beta = deb_mc['marg_beta'].ravel()
                print ('DEBUG[%d]: df(ep_pi)=%.4e, df(ep_beta)=%.4e, ' +
                       'df(m_ep)=%.4e, df(m_beta)=%.4e') % \
                       (res.nit, helpers.maxreldiff(rep.ep_pi,deb_ep_pi),
                        helpers.maxreldiff(rep.ep_beta,deb_ep_beta),
                        helpers.maxreldiff(rep.marg_pi,deb_marg_pi),
                        helpers.maxreldiff(rep.marg_beta,deb_marg_beta))
            # TODO: Plot absolute differences means, stddevs (as in Matlab)
            if res.delta < opts.deltaeps:
                res.rstat = 0
                break
        # Return stuff
        if opts.res_det:
            return (res, res_det)
        else:
            return res

# Testcode (really basic)

#if __name__ == "__main__":
