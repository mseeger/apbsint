"""
utilities
=========

Collection of classes and functions for components and services of ApBsInT:
- Potential manager
- Model (PM and B factor)
- Representation

"""

import numpy as np
import scipy.sparse as ssp
import scipy.linalg as sla
import numbers
import time  # For profiling

import apbsint.helpers as helpers
import apbsint.coup_fact as cf
import apbsint.eptools_ext as epx
import apbsint.ptannotate_ext as pta

__all__ = ['ElemPotManager', 'PotManager', 'Model', 'ModelCoupled',
           'ModelFactorized', 'Representation', 'RepresentationCoupled',
           'RepresentationFactorized']

# Potential manager classes

class ElemPotManager:
    """
    ElemPotManager
    ==============

    Elementary component type for PotManager, collects potential name,
    number of potentials and potential parameters.
    A potential can be annotated by an object of type
    pta.PotentialAnnotation. The default is None (no annotation).
    Attributes should be changed only via access methods: we maintain
    an up_date flag, based on which the internal representation is
    recomputed in PotManager. Otherwise, attributes are not controlled
    here, but only in PotManager.

    The potential type belongs to an argument group. Right now:
    - 0: Univariate t_j(s_j)
    - 1: Bivariate precision t_j(s_j,tau_k(j))
    For group 1, the index vector [k(j)] must be given as 'kind'. The
    validity of the index is checked in 'PotManager.check_internal'.
    """
    def __init__(self,name,size,pars,annobj=None,kind=None):
        self.setname(name,False)
        self.setsize(size,False)
        self.setpars(pars,False)
        self.setannobj(annobj,False)
        self.setkind(kind,False)
        self.up_date = False

    def setname(self,name,chk=True):
        if not isinstance(name,str):
            raise TypeError('NAME must be string')
        if chk:
            self.up_date = self.up_date and (name == self.name)
        self.name = name

    def setsize(self,size,chk=True):
        if not isinstance(size,numbers.Integral) or size<1:
            raise TypeError('SIZE must be positive integer')
        if chk:
            self.up_date = self.up_date and (size == self.size)
        self.size = size

    def setpars(self,pars,chk=True):
        if not isinstance(pars,tuple):
            pars = (pars,)
        # Must be scalar (float), 1D np.ndarray, or tuple of such entries
        for el in pars:
            if not (isinstance(el,float) or
                    (isinstance(el,np.ndarray) and el.ndim == 1)):
                raise TypeError('PARS must scalar, 1D numpy.ndarray, or list/tuple thereof')
        if chk:
            # Can't really check this, so assume it has changed
            self.up_date = False
        self.pars = pars

    def setannobj(self,annobj,chk=True):
        if not (annobj is None or isinstance(annobj,pta.PotentialAnnotation)):
            raise TypeError('ANNOBJ must be apbsint.ptannotate_ext.PotentialAnnotation')
        if chk:
            self.up_date = self.up_date and (annobj == self.annobj)
        self.annobj = annobj

    def setkind(self,kind,chk=True):
        if kind is not None:
            # 'kind' must be array of type np.int32. Neither length nor content
            # is checked here
            if not (isinstance(kind,np.ndarray) and kind.ndim == 1 and
                    kind.dtype == np.int32):
                raise TypeError('KIND must np.ndarray (type np.int32)')
        if chk:
            self.up_date = self.up_date and (kind == self.kind)
        self.kind = kind

# Mechanism to check whether internal representation has to be recomputed:
# - Recompute (in 'check_internal') if the up_date field of any ElemPotManager
#   object is False
# - Set all up_date fields to False upon construction here
# - Set all up_date fields to True in 'check_internal'
class PotManager:
    """
    PotManager
    ==========

    A potential manager consists of one or more (tuple) ElemPotManager
    objects, which are stacked. It also maintains an internal representation.
    Call 'check_internal' before accessing the internal representation: it is
    recomputed whenever any of the ElemPotManager objects has changed.

    If there are bivariate precision potentials (see 'ElemPotManager'), they
    must come last. This is checked in 'check_internal'. In this case, the
    'kind' indexes are joined to form a single flat [k(j)]. An internal
    representation of j <-> k is maintained in 'tauind'. 'num_bvprec' is the
    number of precision potentials, 'num_tau' the size of the tau vector.
    Internal consistency is checked in 'check_internal'.
    """
    def __init__(self,elem):
        if isinstance(elem,list):
            raise TypeError('ELEM must be tuple, not list')
        if not isinstance(elem,tuple):
            elem = (elem,)
        for el in elem:
            if not isinstance(el,ElemPotManager):
                raise TypeError('ELEM must be tuple of ElemPotManager')
        self.size = 0
        for el in elem:
            el.up_date = False
            self.size += el.size
        self.elem = elem

    # Internal representation consists of 'potids', 'numpot', 'parvec',
    # 'parshrd', 'annobj': all contiguous 1D np.ndarray, dtype np.int32
    # except 'parvec': np.float64, 'annobj': np.uint64.
    # 'annobj' stores a void* to the annotation object ('getptr' method),
    # or 0 if none.
    # See src/eptools/potentials/PotManagerFactory.h or documentation for
    # details.
    # Also, we compute 'updind' as index of all non-Gaussian potentials
    # (type other than 'Gaussian').
    def check_internal(self):
        # HIER: bvprec stuff!
        elem = self.elem
        do_recomp = False
        for el in elem:
            if not el.up_date:
                do_recomp = True
                break
        if do_recomp:
            # Loop 1: Everything except PARVEC, determine size
            nb = len(elem)
            self.potids = np.empty(nb,dtype=np.int32)
            self.numpot = np.empty(nb,dtype=np.int32)
            self.annobj = np.zeros(nb,dtype=np.uint64)
            parshrd = []
            pvsz = 0 # Size of PARVEC
            updind = []
            pid_gauss = epx.getpotid('Gaussian')
            off = 0
            self.num_bvprec = 0
            for k in xrange(nb):
                pid = epx.getpotid(elem[k].name)
                if pid == -1:
                    raise ValueError("Block {0}: Unknown potential name '{1}'".format(k,elem[k].name))
                self.potids[k] = pid
                numk = elem[k].size
                self.numpot[k] = numk
                agid = epx.getpotagroup(pid)
                if self.num_bvprec>0 and agid!=1:
                    raise ValueError('Block {0}: Bivariate precision potentials must come last'.format(k))
                if agid==1:
                    self.num_bvprec += numk
                if elem[k].annobj is not None:
                    self.annobj[k] = elem[k].annobj.getptr()
                pars = elem[k].pars
                nump = len(pars)
                for p in xrange(nump):
                    if isinstance(pars[p],float):
                        parshrd.append(1)
                        sz = 1
                    else:
                        sz = len(pars[p])
                        if sz == 1:
                            parshrd.append(1)
                        else:
                            parshrd.append(0)
                            if sz != numk:
                                raise ValueError('Block {0}, pars[{1}]: Wrong size'.format(k,p))
                    pvsz += sz
                if pid != pid_gauss:
                    updind.extend(range(off,off+numk))
                off += numk
            self.parshrd = np.array(parshrd,dtype=np.int32)
            self.updind = np.array(updind,dtype=np.int32)
            # Loop 2: Assemble PARVEC
            self.parvec = np.empty(pvsz,dtype=np.float64)
            off = 0
            for k in xrange(nb):
                pars = elem[k].pars
                nump = len(pars)
                for p in xrange(nump):
                    sz = 1 if isinstance(pars[p],float) else len(pars[p])
                    self.parvec[off:off+sz] = pars[p]
                    off += sz
            # Test whether all parameter values are valid
            msg = epx.potmanager_isvalid(self.potids,self.numpot,self.parvec,
                                         self.parshrd,self.annobj)
            if len(msg)>0:
                raise ValueError(msg)
            # Set all up_date flags
            for el in elem:
                el.up_date = True

    def filterpots(self,nameset):
        """
        Return index of potential positions corresponding to those with
        type names in 'nameset' (set of str)
        """
        if not isinstance(nameset,set):
            raise TypeError('NAMESET must be set of strings')
        off = 0
        res = []
        for el in self.elem:
            numk = el.size
            if el.name in nameset:
                res.extend(range(off,off+numk))
            off += numk
        return np.array(res,dtype=np.int32)

class Model:
    """
    Model
    =====

    Collects B coupling factor and potential manager (type apbsint.PotManager).
    """
    def __init__(self,bfact,potman):
        if not isinstance(potman,PotManager):
            raise TypeError('POTMAN must be instance of apbsint.PotManager')
        self.bfact = bfact
        self.potman = potman

class ModelCoupled(Model):
    """
    ModelCoupled
    ============

    Model for coupled mode. The B factor must be of type apbsint.Mat.
    """
    def __init__(self,bfact,potman):
        if not isinstance(bfact,cf.Mat):
            raise TypeError('BFACT must be instance of apbsint.Mat')
        if bfact.transp:
            raise TypeError('BFACT must not be transposed')
        Model.__init__(self,bfact,potman)
        if bfact.shape(0) != potman.size:
            raise TypeError('BFACT, POTMAN must have same size')

class ModelFactorized(Model):
    """
    ModelFactorized
    ===============

    Model for factorized mode. The B factor must be of type
    apbsint.MatFactorizedInf.
    """
    def __init__(self,bfact,potman):
        if not isinstance(bfact,cf.MatFactorizedInf):
            raise TypeError('BFACT must be instance of apbsint.MatFactorizedInf')
        Model.__init__(self,bfact,potman)
        if bfact.shape(0) != potman.size:
            raise TypeError('BFACT, POTMAN must have same size')

# Representation classes (coupled mode for now)

class Representation:
    """
    Representation
    ==============

    Base class for EP posterior representations. The EP (message)
    parameters are also maintained here.
    """
    def __init__(self,bfact,ep_pi=None,ep_beta=None):
        self.ep_pi = None
        if ep_pi is not None:
            self.setpi(ep_pi)
        self.ep_beta = None
        if ep_beta is not None:
            self.setbeta(ep_beta)
        self.bfact = bfact

    def setpi(self,ep_pi):
        sz = self.size_pars()
        if not helpers.check_vecsize(ep_pi,sz):
            raise TypeError('EP_PI must be vector of size {0}'.format(sz))
        if self.ep_pi is None:
            self.ep_pi = np.empty(sz)
        self.ep_pi[:] = ep_pi

    def setbeta(self,ep_beta):
        sz = self.size_pars()
        if not helpers.check_vecsize(ep_beta,sz):
            raise TypeError('EP_BETA must be vector of size {0}'.format(sz))
        if self.ep_beta is None:
            self.ep_beta = np.empty(sz)
        self.ep_beta[:] = ep_beta

    # Internal methods

    def size_pars(self):
        """
        Returns size of EP parameter vectors ep_pi, ep_beta
        """
        raise NotImplementedError('SIZE_PARS must be implemented')

class RepresentationCoupled(Representation):
    """
    RepresentationCoupled
    =====================

    EP posterior representation in coupled mode. Used for parallel and
    sequential updating EP. We maintain the Cholesky factor L of the
    posterior inverse covariance matrix A: A = L L^T, as well as
      c = L^-1 B^T beta,
    B the factor given in 'bfact' (must be type 'Mat'). If 'keep_margs'
    is True, we also keep marginal moments in 'marg_means', 'marg_vars'.
    """
    def __init__(self,bfact,ep_pi=None,ep_beta=None,keep_margs=False):
        if not isinstance(bfact,cf.Mat):
            raise TypeError('BFACT must be instance of apbsint.Mat')
        Representation.__init__(self,bfact,ep_pi,ep_beta)
        self.keep_margs = keep_margs

    def size_pars(self):
        return self.bfact.shape(0)

    def refresh(self):
        """
        Recompute representation from scratch, given EP paraemeters and B
        coupling factor. If 'keep_margs'==True, the covariance A^-1 is
        computed as byproduct. In this case, A^-1 is kept as attribute
        'post_cov'.
        ATTENTION: 'post_cov' is valid only directly after a call of
        'refresh'. It is not kept up-2-date, and may even be overwritten
        by other methods.
        """
        #t_start0=time.time()
        bfact = self.bfact
        m, n = bfact.shape()
        # Cholesky factor L and c vector
        # We build the A matrix in 'self.lfact'. 'sla.cholesky' overwrites A
        # directly by L.
        # NOTE: 'resize' method does not work as documented. Even with no
        # references to the object, it raises an error, so have to call with
        # 'refcheck=False'
        try:
            self.lfact.resize((n,n),refcheck=False)
        except AttributeError:
            self.lfact = np.empty((n,n))
        #t_start=time.time()
        bfact.mat_btdb(self.ep_pi,self.lfact)
        #t_stop=time.time()
        #print 'Time(refresh::mat_btdb): %.8fs' % (t_stop-t_start)
        self.lfact = sla.cholesky(self.lfact,lower=True,overwrite_a=True)
        self.cvec = sla.solve_triangular(self.lfact,bfact.T().mvm(self.ep_beta),
                                         lower=True,trans='N')
        if self.keep_margs:
            # Recompute marginal moments
            self.marg_means = bfact.mvm(sla.solve_triangular(self.lfact,
                                                             self.cvec,
                                                             lower=True,
                                                             trans='T'))
            # We need A^-1 for the marginal variances. It is written into
            # 'self.post_cov'
            try:
                self.post_cov.resize((n,n),refcheck=False)
            except AttributeError:
                self.post_cov = np.empty((n,n))
            self._comp_inva(self.post_cov)
            #t_start=time.time()
            self.marg_vars = bfact.diag_bsbt(self.post_cov)
            #t_stop=time.time()
            #print 'Time(refresh::diag_bsbt): %.8fs' % (t_stop-t_start)
        # Cholesky up/downdate functions require F contiguous
        # (do this down here, in case other methods called here prefer C
        # contiguous)
        if not self.lfact.flags['F_CONTIGUOUS']:
            # NOTE: Very silly! 'sla.cholesky' produces C contiguous, so here
            # a useless copy is done!
            # We could use 'sla.cholesky' with 'lower=false', then assign
            # the transpose. This would avoid the copy (dodgy?)
            self.lfact = np.asfortranarray(self.lfact)
        #t_stop0=time.time()
        #print 'Time(refresh(ALL)): %.8fs' % (t_stop0-t_start0)

    def update_single(self,j,delpi,delbeta,vvec=None):
        """
        Change of EP parameters:
          ep_pi[j] += delpi; ep_beta[j] += delbeta
        The representation is updated accordingly. In particular, the
        Cholesky factor 'lfact' is updated ('delpi'>0) or downdated
        ('delpi'<0). If 'keep_margs'==True, the marginal moments are
        updated as well.
        In 'vvec', the vector L^-1 B[j,:] can be passed. If not, it is
        recomputed here.
        NOTE: 'post_cov' (if given) is not updated!
        """
        bfact = self.bfact
        m, n = bfact.shape()
        if not (isinstance(j,numbers.Integral) and j>=0 and j<m and
                isinstance(delpi,numbers.Real) and
                isinstance(delbeta,numbers.Real)):
            raise ValueError('J, DELPI or DELBETA wrong')
        if not (vvec is None or helpers.check_vecsize(vvec,n)):
            raise TypeError('VVEC wrong')
        # Scratch variables. We keep them as members, to avoid having to
        # allocate them in every call
        try:
            self.cup_c.resize(n,refcheck=False)
            self.cup_s.resize(n,refcheck=False)
            self.cup_wk.resize(n,refcheck=False)
            self.cup_z.resize((1,n),refcheck=False)
            self.us_bvec.resize(n,refcheck=False)
            if self.keep_margs:
                self.us_wvec.resize(m,refcheck=False)
                self.us_w2vec.resize(m,refcheck=False)
        except AttributeError:
            self.cup_c = np.empty(n)
            self.cup_s = np.empty(n)
            self.cup_wk = np.empty(n)
            self.cup_z = np.empty((1,n),order='F')
            self.us_bvec = np.empty(n)
            if self.keep_margs:
                self.us_wvec = np.empty(m)
                self.us_w2vec = np.empty(m)
        bvec = self.us_bvec
        if self.keep_margs:
            wvec = self.us_wvec
            w2vec = self.us_w2vec
        if delpi>0.:
            # Cholesky update
            tscal = np.sqrt(delpi)
            bfact.T().getcol(j,bvec)
            if self.keep_margs:
                if vvec is None:
                    # Need 'vvec' below, so compute it here
                    vvec = sla.solve_triangular(self.lfact,bvec,lower=True,
                                                trans='N')
                mu = np.inner(vvec,self.cvec)
                rho = np.inner(vvec,vvec)
            bvec *= tscal
            yscal = np.empty(1)
            yscal[0] = delbeta/tscal
            self.cup_z[0] = self.cvec
            stat = epx.choluprk1(self.lfact,'L',bvec,self.cup_c,self.cup_s,
                                 self.cup_wk,self.cup_z,yscal)
            if stat != 0:
                raise sla.LinAlgError("Numerical error in 'choluprk1' (external)")
            self.cvec[:] = self.cup_z.ravel()
        else:
            # Cholesky downdate
            tscal = np.sqrt(-delpi)
            if vvec is None:
                bfact.T().getcol(j,bvec)
                vvec = sla.solve_triangular(self.lfact,bvec,lower=True,
                                            trans='N')
            if self.keep_margs:
                mu = np.inner(vvec,self.cvec)
                rho = np.inner(vvec,vvec)
            yscal = np.empty(1)
            yscal[0] = -delbeta/tscal
            self.cup_z[0] = self.cvec
            bvec[:] = vvec; bvec *= tscal
            stat = epx.choldnrk1(self.lfact,'L',bvec,self.cup_c,self.cup_s,
                                 self.cup_wk,self.cup_z,yscal)
            if stat != 0:
                raise sla.LinAlgError("Numerical error in 'choldnrk1' (external)")
            self.cvec[:] = self.cup_z.ravel()
        self.ep_pi[j] += delpi
        self.ep_beta[j] += delbeta
        if self.keep_margs:
            # Update marginal moments
            assert vvec is not None
            bfact.mvm(sla.solve_triangular(self.lfact,vvec,lower=True,
                                           trans='T'),wvec)
            tscal = 1./(delpi*rho+1.);
            w2vec[:] = wvec; w2vec *= ((delbeta-delpi*mu)*tscal)
            self.marg_means += w2vec
            wvec *= wvec; wvec *= (delpi*tscal)
            self.marg_vars -= wvec

    def get_marg(self,j,vvec=None):
        """
        Returns (mu, rho), mu marginal mean, rho marginal variance at potential
        j (Gaussian marginal, not tilted marginal).
        If 'vvec' is given, L^-1 B[j,:] is written there. In this case, the
        marginal is always computed from scratch. Otherwise, if
        'keep_margs'==True, we use 'marg_XXX'.
        If (mu, rho) are computed from scratch and 'keep_margs'==True, the
        corr. entries of 'marg_XXX' are refreshed.
        """
        bfact = self.bfact
        m, n = bfact.shape()
        if not (isinstance(j,numbers.Integral) and j>=0 and j<m):
            raise ValueError('J wrong')
        if vvec is None:
            if self.keep_margs:
                return (self.marg_means[j], self.marg_vars[j])
            vvec = np.empty(n)
        else:
            if not helpers.check_vecsize(vvec,n):
                raise TypeError('VVEC wrong')
        try:
            self.us_bvec.resize(n,refcheck=False)
        except AttributeError:
            self.us_bvec = np.empty(n)
        bfact.T().getcol(j,self.us_bvec)
        vvec[:] = sla.solve_triangular(self.lfact,self.us_bvec,lower=True,
                                       trans='N')
        mu = np.inner(vvec,self.cvec)
        rho = np.inner(vvec,vvec)
        if self.keep_margs:
            # Refresh entries
            self.marg_means[j] = mu
            self.marg_vars[j] = rho
        return (mu, rho)

    def predict(self,pbfact,pmeans,pvars=None,use_cov=False):
        """
        Compute predictive means (and variances, optional), given test set
        coupling factor B_p (in 'pbfact').
        Predictive variances require A^-1. If 'post_cov' is defined and
        'use_cov'==True, A^-1 is taken from 'post_cov'. Otherwise, it is
        computed here (and written into 'post_cov').
        NOTE: Use 'use_cov'=True if 'refresh' with 'keep_margs'=True has
        been called just before.
        """
        if not isinstance(pbfact,cf.Mat):
            raise TypeError('PBFACT must be instance of apbsint.Mat')
        pm, n = pbfact.shape()
        if n != self.bfact.shape(1):
            raise TypeError('PBFACT has wrong size')
        if not (helpers.check_vecsize(pmeans,pm) and
                (pvars is None or helpers.check_vecsize(pvars,pm))):
            raise TypeError('PMEANS or PVARS wrong')
        # Predictive means
        pbfact.mvm(sla.solve_triangular(self.lfact,self.cvec,lower=True,
                                        trans='T'),pmeans)
        if pvars is not None:
            # Predictive variances: Need inverse A^-1
            try:
                if self.post_cov.shape != (n,n):
                    raise TypeError('Internal error: POST_COV attribute has wrong size')
            except AttributeError:
                if use_cov:
                    raise ValueError('POST_COV is not defined')
                self.post_cov = np.empty((n,n))
            amat = self.post_cov
            if not use_cov:
                self._comp_inva(amat)
            pbfact.diag_bsbt(amat,pvars)

    # Internal methods

    def _comp_inva(self,amat):
        """
        Compute inverse of A and write into 'amat' (must be right size and
        C contiguous). The Cholesky factor of A is in 'lfact'.
        """
        # Set 'amat' to identity matrix
        amat.fill(0.)
        np.fill_diagonal(amat,1.)
        #if not (amat.flags['C_CONTIGUOUS'] or amat.flags['F_CONTIGUOUS']):
        #    raise TypeError('Internal error: AMAT should be contiguous!')
        #aflat = amat.ravel()
        #if aflat.flags['OWNDATA']:
        #    raise TypeError('Internal error: Need view here!')
        #aflat.fill(0.)
        #n = amat.shape[0]
        #aflat[0::n+1] = 1.
        amat[:] = sla.cho_solve((self.lfact,True),amat,overwrite_b=True)

class RepresentationFactorized(Representation):
    """
    RepresentationFactorized
    ========================

    EP posterior representation in factorized mode. B (in 'bfact') is a
    sparse matrix (see apbsint.MatFactorizedInf).
    'ep_pi', 'ep_beta' are the message parameters (flat, size 'bfact.nnz()').
    'marg_pi', 'marg_beta' are the marginals (natural parameters).

    Selective damping is supported if 'sd_numk' is given. The SD
    representation tracks max_k pi_{k,i} for each variable, it is
    initialized/recomputed by 'seldamp_reset'.
    """
    def __init__(self,bfact,ep_pi=None,ep_beta=None):
        if not isinstance(bfact,cf.MatFactorizedInf):
            raise TypeError('BFACT must be apbsint.MatFactorizedInf')
        Representation.__init__(self,bfact,ep_pi,ep_beta)

    def size_pars(self):
        return self.bfact.nnz()

    def refresh(self):
        """
        Recomputes marginals 'marg_pi', 'marg_beta' from message parameters
        'ep_pi', 'ep_beta'.
        """
        bf = self.bfact
        m, n = bf.shape()
        if self.ep_pi is None or self.ep_beta is None:
            raise ValueError('EP parameters must be initialized')
        try:
            self.marg_pi.resize(n,refcheck=False)
            self.marg_beta.resize(n,refcheck=False)
        except AttributeError:
            self.marg_pi = np.empty(n)
            self.marg_beta = np.empty(n)
        epx.fact_compmarginals(n,m,bf.rowind,bf.colind,bf.bvals,self.ep_pi,
                               self.ep_beta,self.marg_pi,self.marg_beta)

    def predict(self,pbfact,pmeans,pvars=None):
        if not isinstance(pbfact,cf.MatFactorizedInf):
            raise TypeError('PBFACT must be apbsint.MatFactorizedInf')
        pm, n = pbfact.shape()
        if n != self.bfact.shape(1):
            raise TypeError('PBFACT has wrong size')
        if not (helpers.check_vecsize(pmeans,pm) and
                (pvars is None or helpers.check_vecsize(pvars,pm))):
            raise TypeError('PMEANS or PVARS wrong')
        tvec = 1./self.marg_pi
        if pvars is not None:
            # 'pbfact.b2fact' is B_test**2
            pvars[:] = pbfact.b2fact.dot(tvec)
        tvec *= self.marg_beta
        if pmeans.flags['C_CONTIGUOUS']:
            pbfact.mvm(tvec,pmeans)
        else:
            pmeans[:] = pbfact.mvm(tvec)

    def seldamp_reset(self,numk,subind=None,subexcl=False):
        """
        Initializes or resets the selective damping (SD) representation. SD
        ensures that cavity marginals are well-defined after each EP update.
        The SD representation tracks max_k pi_{k,i} for each variable, by
        storing the 'numk' largest pi values for each i. The larger 'numk',
        the less often maxima have to be recomputed.
        If 'subind' is given, max_k runs only over this index (if
        'subexcl'==False) or over its complement (if 'subexcl'==True).
        """
        bf = self.bfact
        m, n = bf.shape()
        if not isinstance(numk,numbers.Integral) or numk<2:
            raise TypeError('NUMK must be integer > 1')
        if not (subind is None or
                (helpers.check_vecsize(subind) and subind.dtype == np.int32)):
            raise TypeError('SUBIND must be numpy.ndarray with dtype numpy.int32')
        (self.sd_numvalid, self.sd_topind, self.sd_topval) \
            = epx.fact_compmaxpi(n,m,bf.rowind,bf.colind,bf.bvals,self.ep_pi,
                                 self.ep_beta,numk,subind,subexcl)
        self.sd_subind = subind
        self.sd_subexcl = subexcl
        self.sd_numk = numk

# Testcode (really basic)

if __name__ == "__main__":
    pelem1 = ElemPotManager('Laplace',100,(0., 1.2))
    pelem2 = ElemPotManager('Gaussian',200,(np.random.randn(200), 1.5))
    pelem3 = ElemPotManager('Probit',7,(np.array([1.,-1.,1.,1.,1.,1.,-1.]),
                                        0.))
    pman = PotManager((pelem1, pelem2, pelem3))
    pman.check_internal()
    assert pman.updind == range(100) + range(300,307), \
        'PotManager.check_internal: updind is wrong'
    print('PotManager.check_internal seems OK')
    # HIER: Test code for RepresentationCoupled!
