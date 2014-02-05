# -------------------------------------------------------------------
# CEPTOOLS_EXT
# -------------------------------------------------------------------
# Extern declarations for eptools_ext extension module.
# Author: Matthias Seeger
# -------------------------------------------------------------------

# Declarations: Pointer_to_function types for BLAS functions. Required by
# eptwrap_choluprk1, eptwrap_choldnrk1

cdef extern from "src/eptools/wrap/matrix_types.h":
    ctypedef int blasint_t

    ctypedef void (* dcopy_type) (blasint_t* n,double* x,blasint_t* incx,
                                  double* y,blasint_t* incy)

    ctypedef void (* dscal_type) (blasint_t* n,double* alpha,double* x,
                                  blasint_t* incx)

    ctypedef double (* ddot_type) (blasint_t* n,double* a,blasint_t* lda,
                                   double* b,blasint_t* ldb)

    ctypedef void (* daxpy_type) (blasint_t *n,double* alpha,double* x,
                                  blasint_t* incx,double* y,blasint_t* incy)

    ctypedef void (* drotg_type) (double* a,double* b,double* c,double* s)

    ctypedef void (* drot_type) (blasint_t* n,double* x,blasint_t* incx,
                                 double* y,blasint_t* incy,double* c,double* s)

    ctypedef void (* dtrsv_type) (char* uplo,char* trans,char* diag,
                                  blasint_t* n,double* a,blasint_t* lda,
                                  double* x,blasint_t* incx)

    ctypedef struct fst_matrix:
        double* buff
        int m,n
        int stride
        char strcode[4]

# Declarations: C wrapper functions

cdef extern from "src/eptools/wrap/eptwrap_epupdate_parallel.h":
    void eptwrap_epupdate_parallel(int ain,int aout,int* potids,int npotids,
                                   int* numpot,int nnumpot,double* parvec,
                                   int nparvec,int* parshrd,int nparshrd,
                                   void** annobj,int nannobj,
                                   double* cmu,int ncmu,double* crho,int ncrho,
                                   int* updind,int nupdind,int* rstat,
                                   int nrstat,double* alpha,int nalpha,
                                   double* nu,int nnu,double* logz,int nlogz,
                                   int* errcode,char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_getpotid.h":
    void eptwrap_getpotid(int ain,int aout,char* name,int* pid,int* errcode,
                          char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_getpotname.h":
    void eptwrap_getpotname(int ain,int aout,int pid,char** name,int* errcode,
                            char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_fact_compmarginals.h":
    void eptwrap_fact_compmarginals(int ain,int aout,int n,int m,int* rp_rowind,
                                    int nrp_rowind,int* rp_colind,
                                    int nrp_colind,double* rp_bvals,
                                    int nrp_bvals,double* rp_pi,int nrp_pi,
                                    double* rp_beta,int nrp_beta,
                                    double* margpi,int nmargpi,
                                    double* margbeta,int nmargbeta,
                                    int* errcode,char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_fact_compmaxpi.h":
    void eptwrap_fact_compmaxpi(int ain,int aout,int n,int m,int* rp_rowind,
                                int nrp_rowind,int* rp_colind,int nrp_colind,
                                double* rp_bvals,int nrp_bvals,double* rp_pi,
                                int nrp_pi,double* rp_beta,int nrp_beta,
                                int sd_k,int* sd_subind,int nsd_subind,
                                int sd_subexcl,int* sd_numvalid,
                                int nsd_numvalid,int* sd_topind,int nsd_topind,
                                double* sd_topval,int nsd_topval,int* errcode,
                                char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_fact_sequpdates.h":
    void eptwrap_fact_sequpdates(int ain,int aout,int n,int m,int* updjind,
                                 int nupdjind,int* pm_potids,int npm_potids,
				 int* pm_numpot,int npm_numpot,
                                 double* pm_parvec,int npm_parvec,
                                 int* pm_parshrd,int npm_parshrd,
                                 void** pm_annobj,int npm_annobj,
                                 int* rp_rowind,int nrp_rowind,int* rp_colind,
                                 int nrp_colind,double* rp_bvals,int nrp_bvals,
                                 double* rp_pi,int nrp_pi,double* rp_beta,
                                 int nrp_beta,double* margpi,int nmargpi,
                                 double* margbeta,int nmargbeta,
                                 double piminthres,double dampfact,
                                 int* sd_numvalid,int nsd_numvalid,
                                 int* sd_topind,int nsd_topind,
                                 double* sd_topval,int nsd_topval,
                                 int* sd_subind,int nsd_subind,int sd_subexcl,
                                 int* rstat,int nrstat,double* delta,
                                 int ndelta,double* sd_dampfact,
                                 int nsd_dampfact,int* sd_nupd,int* sd_nrec,
                                 int* errcode,char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_potmanager_isvalid.h":
    void eptwrap_potmanager_isvalid(int ain,int aout,int* potids,int npotids,
                                    int* numpot,int nnumpot,double* parvec,
                                    int nparvec,int* parshrd,int nparshrd,
                                    void** annobj,int nannobj,
                                    int posoff,char** retstr,int* errcode,
                                    char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_epupdate_single.h":
    void eptwrap_epupdate_single1(int ain,int aout,int pid,double* pars,
                                  int npars,void* annobj,double cmu,
                                  double crho,int* rstat,
                                  double* alpha,double* nu,double* logz,
                                  int* errcode,char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_epupdate_single.h":
    void eptwrap_epupdate_single2(int ain,int aout,char* pname,double* pars,
                                  int npars,void* annobj,double cmu,
                                  double crho,int* rstat,
                                  double* alpha,double* nu,double* logz,
                                  int* errcode,char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_epupdate_single.h":
    void eptwrap_epupdate_single3(int ain,int aout,int* potids,int npotids,
                                  int* numpot,int nnumpot,double* parvec,
                                  int nparvec,int* parshrd,int nparshrd,
                                  void** annobj,int nannobj,
                                  int pind,double cmu,double crho,int* rstat,
                                  double* alpha,double* nu,double* logz,
                                  int* errcode,char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_choluprk1.h":
    void eptwrap_choluprk1(int ain,int aout,fst_matrix* lmat,double* vvec,
                           int nvvec,double* cvec,int ncvec,double* svec,
                           int nsvec,double* wkvec,int nwkvec,fst_matrix* zmat,
                           double* yvec,int nyvec,int* stat,dcopy_type f_dcopy,
                           drotg_type f_drotg,drot_type f_drot,int* errcode,
                           char* errstr)

cdef extern from "src/eptools/wrap/eptwrap_choldnrk1.h":
    void eptwrap_choldnrk1(int ain,int aout,fst_matrix* lmat,double* vvec,
                           int nvvec,double* cvec,int ncvec,double* svec,
                           int nsvec,double* wkvec,int nwkvec,int isp,
                           fst_matrix* zmat,double* yvec,int nyvec,int* stat,
                           dcopy_type f_dcopy,dtrsv_type f_dtrsv,
                           ddot_type f_ddot,drotg_type f_drotg,
                           drot_type f_drot,dscal_type f_dscal,
                           daxpy_type f_daxpy,int* errcode,char* errstr)
