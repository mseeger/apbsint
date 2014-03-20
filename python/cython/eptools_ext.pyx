# -------------------------------------------------------------------
# EPTOOLS_EXT
# -------------------------------------------------------------------
# Cython code wrapping external C++ functions.
# Calls functions from src/eptools/wrap (eptools part)
# Author: Matthias Seeger
# -------------------------------------------------------------------

# TODO:
# - Add docstring's (best: use reStructured text)

import cython
import numpy as np
import apbsint.exceptions as exc

cimport numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from ceptools_ext cimport *

# We get pointers to BLAS functions from scipy:
# See mail.scipy.org/pipermail/numpy-discussion/2013-February/065576.html
# ATTENTION: According to Python 2.7.6 documentation, the CObject API is
# deprecated from Python 2.7 on, and one should use Capsules.
# But would this work in the specific SciPy context here?
import scipy
from cpython cimport PyCObject_AsVoidPtr
__import__('scipy.linalg.blas')

# Helper functions

# Converts NumPy array with uint64 entries into Cython array with void*
# entries. The returned array has to be deallocated with 'PyMem_Free'
cdef void** make_voidptr_array(np.ndarray[np.uint64_t,ndim=1] arr):
    cdef int i
    cdef void** ret = <void**>PyMem_Malloc(arr.shape[0]*sizeof(void*))
    for i in range(arr.shape[0]):
        ret[i] = <void*>arr[i]
    return ret

cdef check_contiguous_array(np.ndarray a,bytes nma):
    if not a.flags.c_contiguous:
        raise TypeError('%s must be contiguous array' % nma.upper())

cdef check_contiguous_array_size(np.ndarray a,bytes nma,int sz):
    if not (a.flags.c_contiguous and a.shape[0]==sz):
        raise TypeError('%s must be contiguous array of size %d' %
                        (nma.upper(),sz))

# Cython functions

# rstat, alpha, nu, logz (optional) are return arguments (contiguous vectors
# of same size as cmu, rstat is int32, others are double).
@cython.boundscheck(False)
@cython.wraparound(False)
def epupdate_parallel(np.ndarray[int,ndim=1] potids not None,
                      np.ndarray[int,ndim=1] numpot not None,
                      np.ndarray[np.double_t,ndim=1] parvec not None,
                      np.ndarray[int,ndim=1] parshrd not None,
                      np.ndarray[np.uint64_t,ndim=1] annobj not None,
                      np.ndarray[np.double_t,ndim=1] cmu not None,
                      np.ndarray[np.double_t,ndim=1] crho not None,
                      np.ndarray[int,ndim=1] rstat not None,
                      np.ndarray[np.double_t,ndim=1] alpha not None,
                      np.ndarray[np.double_t,ndim=1] nu not None,
                      np.ndarray[np.double_t,ndim=1] logz = None,
                      np.ndarray[int,ndim=1] updind = None):
    cdef int rsz, errcode, ain, aout
    cdef char errstr[512]
    cdef void** annobj_p
    cdef int logz_n, updind_n
    cdef double* logz_p
    cdef int* updind_p
    # Ensure that input/output arguments are contiguous
    rsz = cmu.shape[0]
    check_contiguous_array_size(rstat,'RSTAT',rsz)
    check_contiguous_array_size(alpha,'ALPHA',rsz)
    check_contiguous_array_size(nu,'NU',rsz)
    if logz is not None:
        check_contiguous_array_size(logz,'LOGZ',rsz)
    potids = np.ascontiguousarray(potids)
    numpot = np.ascontiguousarray(numpot)
    parvec = np.ascontiguousarray(parvec)
    parshrd = np.ascontiguousarray(parshrd)
    cmu = np.ascontiguousarray(cmu)
    crho = np.ascontiguousarray(crho)
    annobj_p = make_voidptr_array(annobj)  # Convert to void* array
    # Create return arguments
    # NOTE: They are now passed as arguments (more efficient)
    #rsz = cmu.shape[0]
    #cdef np.ndarray[int,ndim=1] rstat = np.zeros(rsz,dtype=np.int32)
    #cdef np.ndarray[np.double_t,ndim=1] alpha = np.zeros(rsz,dtype=np.double)
    #cdef np.ndarray[np.double_t,ndim=1] nu = np.zeros(rsz,dtype=np.double)
    #cdef np.ndarray[np.double_t,ndim=1] logz = np.zeros(rsz,dtype=np.double)
    # Call C function
    if updind is None:
        updind_n = 0
        updind_p = NULL
        ain = 7
    else:
        updind = np.ascontiguousarray(updind)
        updind_n = updind.shape[0]
        updind_p = &updind[0]
        ain = 8
    if logz is None:
        logz_n = 0
        logz_p = NULL
        aout = 3
    else:
        logz_n = logz.shape[0]
        logz_p = &logz[0]
        aout = 4
    eptwrap_epupdate_parallel(ain,aout,&potids[0],potids.shape[0],&numpot[0],
                              numpot.shape[0],&parvec[0],parvec.shape[0],
                              &parshrd[0],parshrd.shape[0],annobj_p,
                              annobj.shape[0],&cmu[0],cmu.shape[0],&crho[0],
                              crho.shape[0],updind_p,updind_n,&rstat[0],
                              rstat.shape[0],&alpha[0],alpha.shape[0],&nu[0],
                              nu.shape[0],logz_p,logz_n,&errcode,errstr)
    PyMem_Free(annobj_p)  # Free temp. void* array
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)

@cython.boundscheck(False)
@cython.wraparound(False)
def epupdate_parallel_bvprec(np.ndarray[int,ndim=1] potids not None,
                             np.ndarray[int,ndim=1] numpot not None,
                             np.ndarray[np.double_t,ndim=1] parvec not None,
                             np.ndarray[int,ndim=1] parshrd not None,
                             np.ndarray[np.uint64_t,ndim=1] annobj not None,
                             np.ndarray[np.double_t,ndim=1] cmu not None,
                             np.ndarray[np.double_t,ndim=1] crho not None,
                             np.ndarray[np.double_t,ndim=1] ca not None,
                             np.ndarray[np.double_t,ndim=1] cc not None,
                             np.ndarray[int,ndim=1] rstat not None,
                             np.ndarray[np.double_t,ndim=1] alpha not None,
                             np.ndarray[np.double_t,ndim=1] nu not None,
                             np.ndarray[np.double_t,ndim=1] hata not None,
                             np.ndarray[np.double_t,ndim=1] hatc not None,
                             np.ndarray[np.double_t,ndim=1] logz = None,
                             np.ndarray[int,ndim=1] updind = None):
    cdef int rsz, errcode, ain, aout
    cdef char errstr[512]
    cdef void** annobj_p
    cdef int logz_n, updind_n
    cdef double* logz_p
    cdef int* updind_p
    # Ensure that input/output arguments are contiguous
    rsz = cmu.shape[0]
    check_contiguous_array_size(rstat,'RSTAT',rsz)
    check_contiguous_array_size(alpha,'ALPHA',rsz)
    check_contiguous_array_size(nu,'NU',rsz)
    if logz is not None:
        check_contiguous_array_size(logz,'LOGZ',rsz)
    rsz = ca.shape[0]
    check_contiguous_array_size(hata,'HATA',rsz)
    check_contiguous_array_size(hatc,'HATC',rsz)
    potids = np.ascontiguousarray(potids)
    numpot = np.ascontiguousarray(numpot)
    parvec = np.ascontiguousarray(parvec)
    parshrd = np.ascontiguousarray(parshrd)
    cmu = np.ascontiguousarray(cmu)
    crho = np.ascontiguousarray(crho)
    ca = np.ascontiguousarray(ca)
    cc = np.ascontiguousarray(cc)
    annobj_p = make_voidptr_array(annobj)  # Convert to void* array
    # Call C function
    if updind is None:
        updind_n = 0
        updind_p = NULL
        ain = 9
    else:
        updind = np.ascontiguousarray(updind)
        updind_n = updind.shape[0]
        updind_p = &updind[0]
        ain = 10
    if logz is None:
        logz_n = 0
        logz_p = NULL
        aout = 5
    else:
        logz_n = logz.shape[0]
        logz_p = &logz[0]
        aout = 6
    eptwrap_epupdate_parallel_bvprec(ain,aout,&potids[0],potids.shape[0],
                                     &numpot[0],numpot.shape[0],&parvec[0],
                                     parvec.shape[0],&parshrd[0],
                                     parshrd.shape[0],annobj_p,annobj.shape[0],
                                     &cmu[0],cmu.shape[0],&crho[0],
                                     crho.shape[0],&ca[0],ca.shape[0],&cc[0],
                                     cc.shape[0],updind_p,updind_n,&rstat[0],
                                     rstat.shape[0],&alpha[0],alpha.shape[0],
                                     &nu[0],nu.shape[0],&hata[0],hata.shape[0],
                                     &hatc[0],hatc.shape[0],logz_p,logz_n,
                                     &errcode,errstr)
    PyMem_Free(annobj_p)  # Free temp. void* array
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)

def getpotid(bytes name not None):
    cdef int errcode, pid
    cdef char errstr[512]
    # Call C function
    eptwrap_getpotid(1,1,<char*>name,&pid,&errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return pid

def getpotname(int pid):
    cdef int errcode
    cdef char errstr[512]
    cdef char* name
    # Call C function
    eptwrap_getpotname(1,1,pid,&name,&errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return <bytes>name

def getpotagroup(int pid):
    cdef int errcode, agid
    cdef char errstr[512]
    # Call C function
    eptwrap_getpotagroup(1,1,pid,&agid,&errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return agid

@cython.boundscheck(False)
@cython.wraparound(False)
def fact_compmarginals(int n,int m,np.ndarray[int,ndim=1] rp_rowind not None,
                       np.ndarray[int,ndim=1] rp_colind not None,
                       np.ndarray[np.double_t,ndim=1] rp_bvals not None,
                       np.ndarray[np.double_t,ndim=1] rp_pi not None,
                       np.ndarray[np.double_t,ndim=1] rp_beta not None,
                       np.ndarray[np.double_t,ndim=1] margpi not None,
                       np.ndarray[np.double_t,ndim=1] margbeta not None):
    cdef int errcode
    cdef char errstr[512]
    # Ensure that input/output arguments are contiguous
    check_contiguous_array(margpi,'MARGPI')
    check_contiguous_array(margbeta,'MARGBETA')
    rp_rowind = np.ascontiguousarray(rp_rowind)
    rp_colind = np.ascontiguousarray(rp_colind)
    rp_bvals = np.ascontiguousarray(rp_bvals)
    rp_pi = np.ascontiguousarray(rp_pi)
    rp_beta = np.ascontiguousarray(rp_beta)
    # Call C function
    eptwrap_fact_compmarginals(9,0,n,m,&rp_rowind[0],rp_rowind.shape[0],
                               &rp_colind[0],rp_colind.shape[0],&rp_bvals[0],
                               rp_bvals.shape[0],&rp_pi[0],rp_pi.shape[0],
                               &rp_beta[0],rp_beta.shape[0],&margpi[0],
                               margpi.shape[0],&margbeta[0],margbeta.shape[0],
                               &errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)

@cython.boundscheck(False)
@cython.wraparound(False)
def fact_compmarginals_bvprec(int n,int m,
                              np.ndarray[int,ndim=1] rp_rowind not None,
                              np.ndarray[int,ndim=1] rp_colind not None,
                              np.ndarray[np.double_t,ndim=1] rp_bvals not None,
                              np.ndarray[np.double_t,ndim=1] rp_pi not None,
                              np.ndarray[np.double_t,ndim=1] rp_beta not None,
                              np.ndarray[int,ndim=1] rp_tauind not None,
                              np.ndarray[np.double_t,ndim=1] rp_a not None,
                              np.ndarray[np.double_t,ndim=1] rp_c not None,
                              np.ndarray[np.double_t,ndim=1] margpi not None,
                              np.ndarray[np.double_t,ndim=1] margbeta not None,
                              np.ndarray[np.double_t,ndim=1] marga not None,
                              np.ndarray[np.double_t,ndim=1] margc not None):
    cdef int errcode
    cdef char errstr[512]
    # Ensure that input/output arguments are contiguous
    check_contiguous_array(margpi,'MARGPI')
    check_contiguous_array(margbeta,'MARGBETA')
    check_contiguous_array(marga,'MARGA')
    check_contiguous_array(margc,'MARGC')
    rp_rowind = np.ascontiguousarray(rp_rowind)
    rp_colind = np.ascontiguousarray(rp_colind)
    rp_bvals = np.ascontiguousarray(rp_bvals)
    rp_pi = np.ascontiguousarray(rp_pi)
    rp_beta = np.ascontiguousarray(rp_beta)
    rp_tauind = np.ascontiguousarray(rp_tauind)
    rp_a = np.ascontiguousarray(rp_a)
    rp_c = np.ascontiguousarray(rp_c)
    # Call C function
    eptwrap_fact_compmarginals_bvprec(15,0,n,m,&rp_rowind[0],rp_rowind.shape[0],
                                      &rp_colind[0],rp_colind.shape[0],
                                      &rp_bvals[0],rp_bvals.shape[0],&rp_pi[0],
                                      rp_pi.shape[0],&rp_beta[0],
                                      rp_beta.shape[0],&rp_tauind[0],
                                      rp_tauind.shape[0],&rp_a[0],
                                      rp_a.shape[0],&rp_c[0],rp_c.shape[0],
                                      &margpi[0],margpi.shape[0],&margbeta[0],
                                      margbeta.shape[0],&marga[0],
                                      marga.shape[0],&margc[0],margc.shape[0],
                                      &errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)

@cython.boundscheck(False)
@cython.wraparound(False)
def fact_compmaxpi(int n,int m,np.ndarray[int,ndim=1] rp_rowind not None,
                   np.ndarray[int,ndim=1] rp_colind not None,
                   np.ndarray[np.double_t,ndim=1] rp_bvals not None,
                   np.ndarray[np.double_t,ndim=1] rp_pi not None,
                   np.ndarray[np.double_t,ndim=1] rp_beta not None,
                   int sd_k,np.ndarray[int,ndim=1] sd_subind = None,
                   int sd_subexcl = 0):
    cdef int errcode, rsz, subind_n, ain
    cdef char errstr[512]
    cdef int* subind_p
    # Ensure that input arguments are contiguous
    rp_rowind = np.ascontiguousarray(rp_rowind)
    rp_colind = np.ascontiguousarray(rp_colind)
    rp_bvals = np.ascontiguousarray(rp_bvals)
    rp_pi = np.ascontiguousarray(rp_pi)
    rp_beta = np.ascontiguousarray(rp_beta)
    # Create return arguments
    if sd_k<2:
        raise ValueError('SD_K: Must be >1')
    if n<1:
        raise ValueError('N: Must be positive')
    rsz = n*(sd_k+1)
    cdef np.ndarray[int,ndim=1] sd_numvalid = np.zeros(n,dtype=np.int32)
    cdef np.ndarray[int,ndim=1] sd_topind = np.zeros(rsz,dtype=np.int32)
    cdef np.ndarray[np.double_t,ndim=1] sd_topval = \
        np.zeros(rsz,dtype=np.double)
    # Call C function
    if sd_subind is None:
        subind_n = 0
        subind_p = NULL
        ain = 8
    else:
        sd_subind = np.ascontiguousarray(sd_subind)
        subind_n = sd_subind.shape[0]
        subind_p = &sd_subind[0]
        ain = 10
    eptwrap_fact_compmaxpi(ain,3,n,m,&rp_rowind[0],rp_rowind.shape[0],
                           &rp_colind[0],rp_colind.shape[0],&rp_bvals[0],
                           rp_bvals.shape[0],&rp_pi[0],rp_pi.shape[0],
                           &rp_beta[0],rp_beta.shape[0],sd_k,subind_p,subind_n,
                           sd_subexcl,&sd_numvalid[0],sd_numvalid.shape[0],
                           &sd_topind[0],sd_topind.shape[0],&sd_topval[0],
                           sd_topval.shape[0],&errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return (sd_numvalid,sd_topind,sd_topval)

@cython.boundscheck(False)
@cython.wraparound(False)
def fact_compmaxac(int n,int m,int numk,
                   np.ndarray[int,ndim=1] rp_rowind not None,
                   np.ndarray[int,ndim=1] rp_colind not None,
                   np.ndarray[np.double_t,ndim=1] rp_bvals not None,
                   np.ndarray[np.double_t,ndim=1] rp_pi not None,
                   np.ndarray[np.double_t,ndim=1] rp_beta not None,
                   np.ndarray[int,ndim=1] rp_tauind not None,
                   np.ndarray[np.double_t,ndim=1] rp_a not None,
                   np.ndarray[np.double_t,ndim=1] rp_c not None,
                   int sda_k,int sdc_k):
    cdef int errcode, rsz
    cdef char errstr[512]
    # Ensure that input arguments are contiguous
    rp_rowind = np.ascontiguousarray(rp_rowind)
    rp_colind = np.ascontiguousarray(rp_colind)
    rp_bvals = np.ascontiguousarray(rp_bvals)
    rp_pi = np.ascontiguousarray(rp_pi)
    rp_beta = np.ascontiguousarray(rp_beta)
    rp_tauind = np.ascontiguousarray(rp_tauind)
    rp_a = np.ascontiguousarray(rp_a)
    rp_c = np.ascontiguousarray(rp_c)
    # Create return arguments
    if sda_k<2:
        raise ValueError('SDA_K: Must be >1')
    if sdc_k<2:
        raise ValueError('SDC_K: Must be >1')
    if numk<1:
        raise ValueError('NUMK: Must be positive')
    rsz = numk*(sda_k+1)
    cdef np.ndarray[int,ndim=1] sda_numvalid = np.zeros(numk,dtype=np.int32)
    cdef np.ndarray[int,ndim=1] sda_topind = np.zeros(rsz,dtype=np.int32)
    cdef np.ndarray[np.double_t,ndim=1] sda_topval = \
        np.zeros(rsz,dtype=np.double)
    rsz = numk*(sdc_k+1)
    cdef np.ndarray[int,ndim=1] sdc_numvalid = np.zeros(numk,dtype=np.int32)
    cdef np.ndarray[int,ndim=1] sdc_topind = np.zeros(rsz,dtype=np.int32)
    cdef np.ndarray[np.double_t,ndim=1] sdc_topval = \
        np.zeros(rsz,dtype=np.double)
    # Call C function
    eptwrap_fact_compmaxac(12,6,n,m,&rp_rowind[0],rp_rowind.shape[0],
                           &rp_colind[0],rp_colind.shape[0],&rp_bvals[0],
                           rp_bvals.shape[0],&rp_pi[0],rp_pi.shape[0],
                           &rp_beta[0],rp_beta.shape[0],&rp_tauind[0],
                           rp_tauind.shape[0],&rp_a[0],rp_a.shape[0],&rp_c[0],
                           rp_c.shape[0],sda_k,sdc_k,&sda_numvalid[0],
                           sda_numvalid.shape[0],&sda_topind[0],
                           sda_topind.shape[0],&sda_topval[0],
                           sda_topval.shape[0],&sdc_numvalid[0],
                           sdc_numvalid.shape[0],&sdc_topind[0],
                           sdc_topind.shape[0],&sdc_topval[0],
                           sdc_topval.shape[0],&errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return (sda_numvalid,sda_topind,sda_topval,sdc_numvalid,sdc_topind,
            sdc_topval)

# NOTE: sd_nupd, sd_nrec are returned only if rstat, delta, sd_dampfact and
# sd_numvalid are all given
@cython.boundscheck(False)
@cython.wraparound(False)
def fact_sequpdates(int n,int m,np.ndarray[int,ndim=1] updjind not None,
                    np.ndarray[int,ndim=1] pm_potids not None,
                    np.ndarray[int,ndim=1] pm_numpot not None,
                    np.ndarray[np.double_t,ndim=1] pm_parvec not None,
                    np.ndarray[int,ndim=1] pm_parshrd not None,
                    np.ndarray[np.uint64_t,ndim=1] pm_annobj not None,
                    np.ndarray[int,ndim=1] rp_rowind not None,
                    np.ndarray[int,ndim=1] rp_colind not None,
                    np.ndarray[np.double_t,ndim=1] rp_bvals not None,
                    np.ndarray[np.double_t,ndim=1] rp_pi not None,
                    np.ndarray[np.double_t,ndim=1] rp_beta not None,
                    np.ndarray[np.double_t,ndim=1] margpi not None,
                    np.ndarray[np.double_t,ndim=1] margbeta not None,
                    double piminthres,double dampfact = 0.,
                    np.ndarray[int,ndim=1] rstat = None,
                    np.ndarray[np.double_t,ndim=1] delta = None,
                    np.ndarray[int,ndim=1] sd_numvalid = None,
                    np.ndarray[int,ndim=1] sd_topind = None,
                    np.ndarray[np.double_t,ndim=1] sd_topval = None,
                    np.ndarray[int,ndim=1] sd_subind = None,
                    int sd_subexcl = 0,
                    np.ndarray[np.double_t,ndim=1] sd_dampfact = None):
    cdef int errcode, rsz, sd_nupd, sd_nrec, aout, ain
    cdef char errstr[512]
    cdef void** annobj_p
    cdef int rstat_n, delta_n, numvalid_n, topind_n, topval_n, subind_n
    cdef int dampfact_n
    cdef int* rstat_p
    cdef double* delta_p
    cdef int* numvalid_p
    cdef int* topind_p
    cdef double* topval_p
    cdef int* subind_p
    cdef double* dampfact_p
    # Ensure that input/output arguments are contiguous
    updjind = np.ascontiguousarray(updjind)
    pm_potids = np.ascontiguousarray(pm_potids)
    pm_numpot = np.ascontiguousarray(pm_numpot)
    pm_parvec = np.ascontiguousarray(pm_parvec)
    pm_parshrd = np.ascontiguousarray(pm_parshrd)
    rp_rowind = np.ascontiguousarray(rp_rowind)
    rp_colind = np.ascontiguousarray(rp_colind)
    rp_bvals = np.ascontiguousarray(rp_bvals)
    check_contiguous_array(rp_pi,'RP_PI')
    check_contiguous_array(rp_beta,'RP_BETA')
    check_contiguous_array(margpi,'MARGPI')
    check_contiguous_array(margbeta,'MARGBETA')
    check_contiguous_array(rstat,'RSTAT')
    check_contiguous_array(delta,'DELTA')
    if sd_numvalid is not None:
        check_contiguous_array(sd_numvalid,'SD_NUMVALID')
        if sd_topind is None or sd_topval is None:
            raise ValueError('SD_TOPIND, SD_TOPVAL must be given')
        check_contiguous_array(sd_topind,'SD_TOPIND')
        check_contiguous_array(sd_topval,'SD_TOPVAL')
        if sd_dampfact is not None:
            check_contiguous_array(sd_dampfact,'SD_DAMPFACT')
    # Call C function
    rsz = updjind.shape[0]
    if rsz<1:
        raise ValueError('UPDJIND must not be empty')
    aout = 0
    ain = 17
    rstat_n = 0
    rstat_p = NULL
    delta_n = 0
    delta_p = NULL
    numvalid_n = 0
    numvalid_p = NULL
    topind_n = 0
    topind_p = NULL
    topval_n = 0
    topval_p = NULL
    subind_n = 0
    subind_p = NULL
    dampfact_n = 0
    dampfact_p = NULL
    if rstat is not None:
        rstat_n = rstat.shape[0]
        rstat_p = &rstat[0]
        aout += 1
        if delta is not None:
            delta_n = delta.shape[0]
            delta_p = &delta[0]
            aout += 1
    if sd_numvalid is not None:
        numvalid_n = sd_numvalid.shape[0]
        numvalid_p = &sd_numvalid[0]
        topind_n = sd_topind.shape[0]
        topind_p = &sd_topind[0]
        topval_n = sd_topval.shape[0]
        topval_p = &sd_topval[0]
        ain += 3
        if sd_subind is not None:
            sd_subind = np.ascontiguousarray(sd_subind)
            subind_n = sd_subind.shape[0]
            subind_p = &sd_subind[0]
            ain += 2
        if sd_dampfact is not None:
            dampfact_n = sd_dampfact.shape[0]
            dampfact_p = &sd_dampfact[0]
            if aout==2:
                aout = 5
    annobj_p = make_voidptr_array(pm_annobj)  # Convert to void* array
    eptwrap_fact_sequpdates(ain,aout,n,m,&updjind[0],updjind.shape[0],
                            &pm_potids[0],pm_potids.shape[0],&pm_numpot[0],
                            pm_numpot.shape[0],&pm_parvec[0],
                            pm_parvec.shape[0],&pm_parshrd[0],
                            pm_parshrd.shape[0],annobj_p,pm_annobj.shape[0],
                            &rp_rowind[0],rp_rowind.shape[0],&rp_colind[0],
                            rp_colind.shape[0],&rp_bvals[0],rp_bvals.shape[0],
                            &rp_pi[0],rp_pi.shape[0],&rp_beta[0],
                            rp_beta.shape[0],&margpi[0],margpi.shape[0],
                            &margbeta[0],margbeta.shape[0],piminthres,dampfact,
                            numvalid_p,numvalid_n,topind_p,topind_n,topval_p,
                            topval_n,subind_p,subind_n,sd_subexcl,rstat_p,
                            rstat_n,delta_p,delta_n,dampfact_p,dampfact_n,
                            &sd_nupd,&sd_nrec,&errcode,errstr)
    PyMem_Free(annobj_p)  # Free temp. void* array
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    if aout>2:
        return (sd_nupd,sd_nrec)

# NOTE: sd_nupd, sd_nrec are returned only if rstat, delta, sd_dampfact and
# sd_numvalid are all given
@cython.boundscheck(False)
@cython.wraparound(False)
def fact_sequpdates_bvprec(int n,int m,np.ndarray[int,ndim=1] updjind not None,
                           np.ndarray[int,ndim=1] pm_potids not None,
                           np.ndarray[int,ndim=1] pm_numpot not None,
                           np.ndarray[np.double_t,ndim=1] pm_parvec not None,
                           np.ndarray[int,ndim=1] pm_parshrd not None,
                           np.ndarray[np.uint64_t,ndim=1] pm_annobj not None,
                           np.ndarray[int,ndim=1] rp_rowind not None,
                           np.ndarray[int,ndim=1] rp_colind not None,
                           np.ndarray[np.double_t,ndim=1] rp_bvals not None,
                           np.ndarray[np.double_t,ndim=1] rp_pi not None,
                           np.ndarray[np.double_t,ndim=1] rp_beta not None,
                           np.ndarray[int,ndim=1] rp_tauind not None,
                           np.ndarray[np.double_t,ndim=1] rp_a not None,
                           np.ndarray[np.double_t,ndim=1] rp_c not None,
                           np.ndarray[np.double_t,ndim=1] margpi not None,
                           np.ndarray[np.double_t,ndim=1] margbeta not None,
                           np.ndarray[np.double_t,ndim=1] marga not None,
                           np.ndarray[np.double_t,ndim=1] margc not None,
                           double piminthres,double aminthres,double cminthres,
                           double dampfact = 0.,
                           np.ndarray[int,ndim=1] rstat = None,
                           np.ndarray[np.double_t,ndim=1] delta = None,
                           np.ndarray[int,ndim=1] sd_numvalid = None,
                           np.ndarray[int,ndim=1] sd_topind = None,
                           np.ndarray[np.double_t,ndim=1] sd_topval = None,
                           np.ndarray[int,ndim=1] sda_numvalid = None,
                           np.ndarray[int,ndim=1] sda_topind = None,
                           np.ndarray[np.double_t,ndim=1] sda_topval = None,
                           np.ndarray[int,ndim=1] sdc_numvalid = None,
                           np.ndarray[int,ndim=1] sdc_topind = None,
                           np.ndarray[np.double_t,ndim=1] sdc_topval = None,
                           np.ndarray[int,ndim=1] sd_subind = None,
                           int sd_subexcl = 0,
                           np.ndarray[np.double_t,ndim=1] sd_dampfact = None):
    cdef int errcode, rsz, sd_nupd, sd_nrec, aout, ain
    cdef char errstr[512]
    cdef void** annobj_p
    cdef int rstat_n, delta_n, numvalid_n, topind_n, topval_n, subind_n
    cdef int dampfact_n, anumvalid_n, atopind_n, atopval_n, cnumvalid_n
    cdef int ctopind_n, ctopval_n
    cdef int* rstat_p
    cdef double* delta_p
    cdef int* numvalid_p
    cdef int* topind_p
    cdef double* topval_p
    cdef int* anumvalid_p
    cdef int* atopind_p
    cdef double* atopval_p
    cdef int* cnumvalid_p
    cdef int* ctopind_p
    cdef double* ctopval_p
    cdef int* subind_p
    cdef double* dampfact_p
    # Ensure that input/output arguments are contiguous
    updjind = np.ascontiguousarray(updjind)
    pm_potids = np.ascontiguousarray(pm_potids)
    pm_numpot = np.ascontiguousarray(pm_numpot)
    pm_parvec = np.ascontiguousarray(pm_parvec)
    pm_parshrd = np.ascontiguousarray(pm_parshrd)
    rp_rowind = np.ascontiguousarray(rp_rowind)
    rp_colind = np.ascontiguousarray(rp_colind)
    rp_bvals = np.ascontiguousarray(rp_bvals)
    rp_tauind = np.ascontiguousarray(rp_tauind)
    check_contiguous_array(rp_pi,'RP_PI')
    check_contiguous_array(rp_beta,'RP_BETA')
    check_contiguous_array(rp_a,'RP_A')
    check_contiguous_array(rp_c,'RP_C')
    check_contiguous_array(margpi,'MARGPI')
    check_contiguous_array(margbeta,'MARGBETA')
    check_contiguous_array(margpi,'MARGA')
    check_contiguous_array(margpi,'MARGC')
    check_contiguous_array(rstat,'RSTAT')
    check_contiguous_array(delta,'DELTA')
    if sd_numvalid is not None:
        check_contiguous_array(sd_numvalid,'SD_NUMVALID')
        if sd_topind is None or sd_topval is None:
            raise ValueError('SD_TOPIND, SD_TOPVAL must be given')
        check_contiguous_array(sd_topind,'SD_TOPIND')
        check_contiguous_array(sd_topval,'SD_TOPVAL')
        if sd_dampfact is not None:
            check_contiguous_array(sd_dampfact,'SD_DAMPFACT')
        if sda_numvalid is not None:
            check_contiguous_array(sda_numvalid,'SDA_NUMVALID')
            if sda_topind is None or sda_topval is None:
                raise ValueError('SDA_TOPIND, SDA_TOPVAL must be given')
            check_contiguous_array(sda_topind,'SDA_TOPIND')
            check_contiguous_array(sda_topval,'SDA_TOPVAL')
        if sdc_numvalid is not None:
            check_contiguous_array(sdc_numvalid,'SDC_NUMVALID')
            if sdc_topind is None or sdc_topval is None:
                raise ValueError('SDC_TOPIND, SDC_TOPVAL must be given')
            check_contiguous_array(sdc_topind,'SDC_TOPIND')
            check_contiguous_array(sdc_topval,'SDC_TOPVAL')
    # Call C function
    rsz = updjind.shape[0]
    if rsz<1:
        raise ValueError('UPDJIND must not be empty')
    aout = 0
    ain = 24
    rstat_n = 0
    rstat_p = NULL
    delta_n = 0
    delta_p = NULL
    numvalid_n = 0
    numvalid_p = NULL
    topind_n = 0
    topind_p = NULL
    topval_n = 0
    topval_p = NULL
    anumvalid_n = 0
    anumvalid_p = NULL
    atopind_n = 0
    atopind_p = NULL
    atopval_n = 0
    atopval_p = NULL
    cnumvalid_n = 0
    cnumvalid_p = NULL
    ctopind_n = 0
    ctopind_p = NULL
    ctopval_n = 0
    ctopval_p = NULL
    subind_n = 0
    subind_p = NULL
    dampfact_n = 0
    dampfact_p = NULL
    if rstat is not None:
        rstat_n = rstat.shape[0]
        rstat_p = &rstat[0]
        aout += 1
        if delta is not None:
            delta_n = delta.shape[0]
            delta_p = &delta[0]
            aout += 1
    if sd_numvalid is not None:
        numvalid_n = sd_numvalid.shape[0]
        numvalid_p = &sd_numvalid[0]
        topind_n = sd_topind.shape[0]
        topind_p = &sd_topind[0]
        topval_n = sd_topval.shape[0]
        topval_p = &sd_topval[0]
        ain += 3
        if sda_numvalid is not None:
            anumvalid_n = sda_numvalid.shape[0]
            anumvalid_p = &sda_numvalid[0]
            atopind_n = sda_topind.shape[0]
            atopind_p = &sda_topind[0]
            atopval_n = sda_topval.shape[0]
            atopval_p = &sda_topval[0]
            ain += 3
        if sdc_numvalid is not None:
            cnumvalid_n = sdc_numvalid.shape[0]
            cnumvalid_p = &sdc_numvalid[0]
            ctopind_n = sdc_topind.shape[0]
            ctopind_p = &sdc_topind[0]
            ctopval_n = sdc_topval.shape[0]
            ctopval_p = &sdc_topval[0]
            ain += 3
        if sd_subind is not None:
            sd_subind = np.ascontiguousarray(sd_subind)
            subind_n = sd_subind.shape[0]
            subind_p = &sd_subind[0]
            ain += 2
        if sd_dampfact is not None:
            dampfact_n = sd_dampfact.shape[0]
            dampfact_p = &sd_dampfact[0]
            if aout==2:
                aout = 5
    annobj_p = make_voidptr_array(pm_annobj)  # Convert to void* array
    eptwrap_fact_sequpdates_bvprec(ain,aout,n,m,&updjind[0],updjind.shape[0],
                                   &pm_potids[0],pm_potids.shape[0],
                                   &pm_numpot[0],pm_numpot.shape[0],
                                   &pm_parvec[0],pm_parvec.shape[0],
                                   &pm_parshrd[0],pm_parshrd.shape[0],
                                   annobj_p,pm_annobj.shape[0],&rp_rowind[0],
                                   rp_rowind.shape[0],&rp_colind[0],
                                   rp_colind.shape[0],&rp_bvals[0],
                                   rp_bvals.shape[0],&rp_pi[0],rp_pi.shape[0],
                                   &rp_beta[0],rp_beta.shape[0],&rp_tauind[0],
                                   rp_tauind.shape[0],&rp_a[0],rp_a.shape[0],
                                   &rp_c[0],rp_c.shape[0],&margpi[0],
                                   margpi.shape[0],&margbeta[0],
                                   margbeta.shape[0],&marga[0],marga.shape[0],
                                   &margc[0],margc.shape[0],piminthres,
                                   aminthres,cminthres,dampfact,numvalid_p,
                                   numvalid_n,topind_p,topind_n,topval_p,
                                   topval_n,anumvalid_p,anumvalid_n,atopind_p,
                                   atopind_n,atopval_p,atopval_n,cnumvalid_p,
                                   cnumvalid_n,ctopind_p,ctopind_n,ctopval_p,
                                   ctopval_n,subind_p,subind_n,sd_subexcl,
                                   rstat_p,rstat_n,delta_p,delta_n,dampfact_p,
                                   dampfact_n,&sd_nupd,&sd_nrec,&errcode,
                                   errstr)
    PyMem_Free(annobj_p)  # Free temp. void* array
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    if aout>2:
        return (sd_nupd,sd_nrec)

@cython.boundscheck(False)
@cython.wraparound(False)
def potmanager_isvalid(np.ndarray[int,ndim=1] potids not None,
                       np.ndarray[int,ndim=1] numpot not None,
                       np.ndarray[np.double_t,ndim=1] parvec not None,
                       np.ndarray[int,ndim=1] parshrd not None,
                       np.ndarray[np.uint64_t,ndim=1] annobj not None,
                       int posoff = 0):
    cdef int errcode
    cdef char errstr[512]
    cdef char* retstr
    cdef void** annobj_p
    # Ensure that input arguments are contiguous
    potids = np.ascontiguousarray(potids)
    numpot = np.ascontiguousarray(numpot)
    parvec = np.ascontiguousarray(parvec)
    parshrd = np.ascontiguousarray(parshrd)
    annobj_p = make_voidptr_array(annobj)  # Convert to void* array
    # Call C function
    eptwrap_potmanager_isvalid(6,1,&potids[0],potids.shape[0],&numpot[0],
                               numpot.shape[0],&parvec[0],parvec.shape[0],
                               &parshrd[0],parshrd.shape[0],annobj_p,
                               annobj.shape[0],posoff,&retstr,&errcode,errstr)
    PyMem_Free(annobj_p)  # Free temp. void* array
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return <bytes>retstr

def epupdate_single(pid,np.ndarray[np.double_t,ndim=1] pars not None,
                    np.uint64_t annobj,double cmu,double crho):
    cdef int errcode, rstat
    cdef char errstr[512]
    cdef double alpha, nu, logz
    # Ensure that input arguments are contiguous
    pars = np.ascontiguousarray(pars)
    # Call C function
    if isinstance(pid,str):
        eptwrap_epupdate_single2(5,4,<char*>pid,&pars[0],pars.shape[0],
                                 <void*>annobj,cmu,crho,&rstat,&alpha,&nu,
                                 &logz,&errcode,errstr)
    else:
        eptwrap_epupdate_single1(5,4,<int>pid,&pars[0],pars.shape[0],
                                 <void*>annobj,cmu,crho,&rstat,&alpha,&nu,
                                 &logz,&errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return (rstat,alpha,nu,logz)

@cython.boundscheck(False)
@cython.wraparound(False)
def epupdate_single_pman(np.ndarray[int,ndim=1] potids not None,
                         np.ndarray[int,ndim=1] numpot not None,
                         np.ndarray[np.double_t,ndim=1] parvec not None,
                         np.ndarray[int,ndim=1] parshrd not None,
                         np.ndarray[np.uint64_t,ndim=1] annobj not None,
                         int pind,double cmu,double crho):
    cdef int errcode, rstat
    cdef char errstr[512]
    cdef void** annobj_p
    cdef double alpha, nu, logz
    # Ensure that input arguments are contiguous
    potids = np.ascontiguousarray(potids)
    numpot = np.ascontiguousarray(numpot)
    parvec = np.ascontiguousarray(parvec)
    parshrd = np.ascontiguousarray(parshrd)
    annobj_p = make_voidptr_array(annobj)  # Convert to void* array
    # Call C function
    eptwrap_epupdate_single3(8,4,&potids[0],potids.shape[0],&numpot[0],
                             numpot.shape[0],&parvec[0],parvec.shape[0],
                             &parshrd[0],parshrd.shape[0],annobj_p,
                             annobj.shape[0],pind,cmu,crho,&rstat,&alpha,&nu,
                             &logz,&errcode,errstr)
    PyMem_Free(annobj_p)  # Free temp. void* array
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return (rstat,alpha,nu,logz)

@cython.boundscheck(False)
@cython.wraparound(False)
def choluprk1(np.ndarray[np.double_t,ndim=2] l not None,bytes luplo not None,
              np.ndarray[np.double_t,ndim=1] vec not None,
              np.ndarray[np.double_t,ndim=1] cvec not None,
              np.ndarray[np.double_t,ndim=1] svec not None,
              np.ndarray[np.double_t,ndim=1] workv not None,
              np.ndarray[np.double_t,ndim=2] z = None,
              np.ndarray[np.double_t,ndim=1] y = None):
    cdef int errcode, stat
    cdef char errstr[512]
    # Ensure that input/output arguments are contiguous
    if not l.flags.f_contiguous:
        raise TypeError('L must be Fortran contiguous (column-major)')
    if not vec.flags.c_contiguous:
        raise TypeError('VEC must be contiguous array')
    if not cvec.flags.c_contiguous:
        raise TypeError('CVEC must be contiguous array')
    if not svec.flags.c_contiguous:
        raise TypeError('SVEC must be contiguous array')
    if not workv.flags.c_contiguous:
        raise TypeError('WORKV must be contiguous array')
    # We use fst_matrix transfer type
    cdef fst_matrix lmat, zmat
    lmat.buff = &l[0,0]
    lmat.m, lmat.n = l.shape[0], l.shape[1]
    lmat.stride = l.shape[0]
    lmat.strcode[0] = luplo[0]; lmat.strcode[1] = 0
    lmat.strcode[2] = 'N'; lmat.strcode[3] = 0
    # Get BLAS function pointers from scipy
    cdef dcopy_type f_dcopy = <dcopy_type>PyCObject_AsVoidPtr( \
        scipy.linalg.blas.cblas.dcopy._cpointer)
    cdef drotg_type f_drotg = <drotg_type>PyCObject_AsVoidPtr( \
        scipy.linalg.blas.cblas.drotg._cpointer)
    cdef drot_type f_drot   = <drot_type>PyCObject_AsVoidPtr( \
        scipy.linalg.blas.cblas.drot._cpointer)
    # Call C function
    if z is None:
        eptwrap_choluprk1(5,1,&lmat,&vec[0],vec.shape[0],&cvec[0],
                          cvec.shape[0],&svec[0],svec.shape[0],&workv[0],
                          workv.shape[0],NULL,NULL,0,&stat,f_dcopy,f_drotg,
                          f_drot,&errcode,errstr)
    else:
        if not z.flags.f_contiguous:
            raise TypeError('Z must be Fortran contiguous (column-major)')
        if not y.flags.c_contiguous:
            raise TypeError('Y must be contiguous array')
        zmat.buff = &z[0,0]
        zmat.m, zmat.n = z.shape[0], z.shape[1]
        zmat.stride = z.shape[0]
        # Not used:
        zmat.strcode[0] = ' '; zmat.strcode[1] = 0
        zmat.strcode[2] = ' '; zmat.strcode[3] = 0
        eptwrap_choluprk1(7,1,&lmat,&vec[0],vec.shape[0],&cvec[0],
                          cvec.shape[0],&svec[0],svec.shape[0],&workv[0],
                          workv.shape[0],&zmat,&y[0],y.shape[0],&stat,
                          f_dcopy,f_drotg,f_drot,&errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return stat

# NOTE: The ISP argument, which can be used for the MEX version, is not
# exposed here, but is fixed to 1 (so VEC must contain the p vector
# already). This is because BLAS dtrsv cannot be accessed via scipy
@cython.boundscheck(False)
@cython.wraparound(False)
def choldnrk1(np.ndarray[np.double_t,ndim=2] l not None,bytes luplo not None,
              np.ndarray[np.double_t,ndim=1] vec not None,
              np.ndarray[np.double_t,ndim=1] cvec not None,
              np.ndarray[np.double_t,ndim=1] svec not None,
              np.ndarray[np.double_t,ndim=1] workv not None,
              np.ndarray[np.double_t,ndim=2] z = None,
              np.ndarray[np.double_t,ndim=1] y = None):
    cdef int errcode, stat, isp
    cdef char errstr[512]
    isp = 1
    # Ensure that input/output arguments are contiguous
    if not l.flags.f_contiguous:
        raise TypeError('L must be Fortran contiguous (column-major)')
    if not vec.flags.c_contiguous:
        raise TypeError('VEC must be contiguous array')
    if not cvec.flags.c_contiguous:
        raise TypeError('CVEC must be contiguous array')
    if not svec.flags.c_contiguous:
        raise TypeError('SVEC must be contiguous array')
    if not workv.flags.c_contiguous:
        raise TypeError('WORKV must be contiguous array')
    # We use fst_matrix transfer type
    cdef fst_matrix lmat, zmat
    lmat.buff = &l[0,0]
    lmat.m, lmat.n = l.shape[0], l.shape[1]
    lmat.stride = l.shape[0]
    lmat.strcode[0] = luplo[0]; lmat.strcode[1] = 0
    lmat.strcode[2] = 'N'; lmat.strcode[3] = 0
    # Get BLAS function pointers from scipy
    cdef dcopy_type f_dcopy = <dcopy_type>PyCObject_AsVoidPtr( \
        scipy.linalg.blas.cblas.dcopy._cpointer)
    cdef ddot_type f_ddot   = <ddot_type>PyCObject_AsVoidPtr( \
        scipy.linalg.blas.cblas.ddot._cpointer)
    cdef drotg_type f_drotg = <drotg_type>PyCObject_AsVoidPtr( \
        scipy.linalg.blas.cblas.drotg._cpointer)
    cdef drot_type f_drot   = <drot_type>PyCObject_AsVoidPtr( \
        scipy.linalg.blas.cblas.drot._cpointer)
    cdef dscal_type f_dscal = <dscal_type>PyCObject_AsVoidPtr( \
        scipy.linalg.blas.cblas.dscal._cpointer)
    cdef daxpy_type f_daxpy = <daxpy_type>PyCObject_AsVoidPtr( \
        scipy.linalg.blas.cblas.daxpy._cpointer)
    # Call C function
    if z is None:
        eptwrap_choldnrk1(6,1,&lmat,&vec[0],vec.shape[0],&cvec[0],
                          cvec.shape[0],&svec[0],svec.shape[0],&workv[0],
                          workv.shape[0],isp,NULL,NULL,0,&stat,f_dcopy,NULL,
                          f_ddot,f_drotg,f_drot,f_dscal,f_daxpy,&errcode,
                          errstr)
    else:
        if not z.flags.f_contiguous:
            raise TypeError('Z must be Fortran contiguous (column-major)')
        if not y.flags.c_contiguous:
            raise TypeError('Y must be contiguous array')
        zmat.buff = &z[0,0]
        zmat.m, zmat.n = z.shape[0], z.shape[1]
        zmat.stride = z.shape[0]
        # Not used:
        zmat.strcode[0] = ' '; zmat.strcode[1] = 0
        zmat.strcode[2] = ' '; zmat.strcode[3] = 0
        eptwrap_choldnrk1(8,1,&lmat,&vec[0],vec.shape[0],&cvec[0],
                          cvec.shape[0],&svec[0],svec.shape[0],&workv[0],
                          workv.shape[0],isp,&zmat,&y[0],y.shape[0],&stat,
                          f_dcopy,NULL,f_ddot,f_drotg,f_drot,f_dscal,f_daxpy,
                          &errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
    return stat

def debug_castannobj(np.uint64_t annobj):
    cdef int errcode
    cdef char errstr[512]
    # Call C function
    eptwrap_debug_castannobj(<void*>annobj,&errcode,errstr)
    # Check for error, raise exception
    if errcode != 0:
        raise exc.ApBsWrapError(<bytes>errstr)
