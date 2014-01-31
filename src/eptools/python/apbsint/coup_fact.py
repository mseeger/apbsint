"""
coup_fact
=========

HIER: Docstring for coup_fact package. Defines base class Mat as well as
basic subclasses of the B coupling factor hierarchy.

"""

import numpy as np
import scipy.sparse as ssp
import numbers

import apbsint.helpers as helpers

__all__ = ['Mat', 'MatDef', 'MatSparse', 'MatDiag', 'MatEye', 'MatSub',
           'MatContainer', 'MatFactorizedInf']

# Mat class hierarchy

class Mat:
    """
    Mat
    ===

    Base class for B coupling factors. B is m-by-n. If 'transp' is True,
    the n-by-m factor B^T is represented.
    The constructor doubles as copy constructor (shallow copy), mainly
    to implement the transpose method 'T'.
    """
    def __init__(self,m,n=None):
        """
        Constructor doubles as copy constructor. May not be very
        "Pythonic", but we need something like this (shallow copy)
        to implement things like T (transpose).
        NOTE: This is a shallow copy, for temporary objects such as
        returned by T.
        """
        if isinstance(m,Mat):
            # Copy constructor
            self.m = m.m
            self.n = m.n
            self.transp = m.transp
        else:
            if (not isinstance(m,numbers.Integral) or n is None or
                not isinstance(n,numbers.Integral) or m<=0 or n<=0):
                raise ValueError('M, N wrong')
            self.m = m; self.n = n
            self.transp = False

    def shape(self,id=None):
        if self.transp:
            sz = (self.n, self.m)
        else:
            sz = (self.m, self.n)
        if id is None:
            return sz
        else:
            return sz[id]

    def T(self):
        # Make a (shallow) copy of 'self'. Needs copy constructor
        tp = self.__class__
        ret = tp(self)
        ret.transp = not self.transp
        return ret

    # Default implementations (using 'mvm')

    def mvm(self,v,out=None):
        """
        Returns matrix-vector product B*v. If 'out' is given, it must be
        a contiguous vector of the correct size, the result is written into
        its buffer. Otherwise, a new vector is created.
        """
        raise NotImplementedError("Method MVM must be implemented")

    def getcol(self,i,out=None):
        """
        Returns i-th column of B. If 'out' is given, result is written
        into its buffer.
        """
        sz = self.shape(1)
        if not isinstance(i,numbers.Integral) or i<0 or i>=sz:
            raise IndexError('I')
        v = np.zeros(sz)
        v[i] = 1.
        return self.mvm(v,out)

    def mat_btdb(self,v,out=None):
        """
        Returns matrix B^T (diag v) B.
        If 'out' is given, the result is written there directly. 'out' must
        have exactly the right size and must be C contiguous, otherwise an
        exception is thrown.
        """
        m, n = self.shape()
        if not helpers.check_vecsize(v,m):
            raise TypeError('V wrong')
        if out is None:
            out = np.empty((n,n))
        else:
            self._matbtdb_checkout(out)
        bt = self.T()
        tv1 = np.zeros(n)
        tv2 = np.empty(m)
        for i in xrange(n):
            tv1[i] = 1.
            self.mvm(tv1,tv2)
            tv1[i] = 0.
            tv2 *= v
            bt.mvm(tv2,out[i])
        return out

    def diag_bsbt(self,s,out=None):
        """
        Returns vector diag(B S B^T), where S is a symmetric matrix. If
        'out' is given, the result is written into its buffer.
        """
        m, n = self.shape()
        out = self._check_resultvec(out,m)
        if (not isinstance(s,np.ndarray) or s.ndim != 2 or
            s.shape[0] != n or s.shape[1] != n):
            raise TypeError('S wrong')
        out.fill(0.)
        tv1 = np.zeros(n)
        tv2 = np.empty(m)
        tv3 = np.empty(m)
        for i in xrange(n):
            tv1[i] = 1.
            self.mvm(tv1,tv2)
            tv1[i] = 0.
            self.mvm(s[i],tv3)
            tv2 *= tv3
            out += tv2
        return out

    # Internal methods (helpers for implementations in subclasses)

    def _check_resultvec(self,out,sz):
        """
        If 'out' is not None, checks that it is a contiguous vector of
        size 'sz'. Otherwise, a vector of this size is created. In any
        case, the vector is returned.
        """
        if out is None:
            out = np.empty(sz)
        elif not (helpers.check_vecsize(out,sz) and
                  out.flags['C_CONTIGUOUS']):
            raise TypeError('OUT argument wrong (must be contiguous)')
        return out

    def _mvm_checkin(self,v,out):
        """
        Tests input arguments to 'mvm'. Returns buffer vector for result,
        which is 'out' if given, otherwise a new vector.
        """
        m, n = self.shape()
        if not helpers.check_vecsize(v,n):
            raise TypeError('V wrong')
        return self._check_resultvec(out,m)

    def _matbtdb_checkout(self,out):
        """
        'out' must be n-by-n matrix, C contiguous
        """
        if not (isinstance(out,np.ndarray) and out.ndim == 2 and
                out.shape == (self.n,self.n) and
                (out.flags['C_CONTIGUOUS'] or out.flags['F_CONTIGUOUS'])):
            raise TypeError('OUT wrong type or size (internal!)')

    def _matbtdb_setdgout(self,out,v):
        """
        Set 'out' equal to diag(v)
        """
        self._matbtdb_checkout(out)
        oflat = out.ravel()
        if oflat.flags['OWNDATA']:
            raise TypeError('Internal error: Need view here!')
        oflat.fill(0.)
        oflat[0::self.n+1] = v

# We use an internal buffer matrix 'buffmat1', of same size as 'mx'. We are
# then careful to avoid the creation of intermediate m-by-n matrices, using
# 'buffmat1' instead.
class MatDef(Mat):
    """
    MatDef
    ======

    MatDef wraps normal matrix (2D numpy.array) as coupling factor.
    """
    def __init__(self,mx):
        if isinstance(mx,MatDef):
            Mat.__init__(self,mx)
            self.mx = mx.mx
            self.buffmat1 = mx.buffmat1
        else:
            if not isinstance(mx,np.ndarray) or mx.ndim != 2:
                raise TypeError('MX wrong type')
            m, n = mx.shape
            Mat.__init__(self,m,n)
            self.mx = mx
            # Buffer matrix for computations (avoids reallocations)
            self.buffmat1 = np.empty_like(mx)

    def mvm(self,v,out=None):
        out = self._mvm_checkin(v,out)
        if self.transp:
            mx = self.mx.T
        else:
            mx = self.mx
        np.dot(mx,v,out)
        return out

    def getcol(self,i,out=None):
        out = self._check_resultvec(out,self.shape(0))
        if self.transp:
            out[:] = self.mx[i]
        else:
            out[:] = self.mx[:,i]
        return out

    def mat_btdb(self,v,out=None):
        if self.transp:
            raise NotImplementedError("NOT IMPLEMENTED")
        self.buffmat1[:] = self.mx
        self.buffmat1 *= v.reshape(-1,1)
        if out is None:
            return np.dot(self.mx.T,self.buffmat1)
        else:
            self._matbtdb_checkout(out)
            return np.dot(self.mx.T,self.buffmat1,out)

    def diag_bsbt(self,s,out=None):
        out = self._check_resultvec(out,self.shape(0))
        if self.transp:
            raise NotImplementedError("NOT IMPLEMENTED")
        np.dot(self.mx,s,self.buffmat1)
        self.buffmat1 *= self.mx
        return np.sum(self.buffmat1,1,out=out)

class MatSparse(Mat):
    """
    MatSparse
    =========

    MatSparse wraps sparse matrix (scipy.sparse, either csr_matrix or
    csc_matrix) as coupling factor.
    """
    def __init__(self,mx):
        if isinstance(mx,MatSparse):
            Mat.__init__(self,mx)
            self.mx = mx.mx
        else:
            if (not isinstance(mx,ssp.csr_matrix) and
                not isinstance(mx,ssp.csc_matrix)):
                raise TypeError('MX must be scipy.sparse.cs{r|c}_matrix')
            m, n = mx.shape
            Mat.__init__(self,m,n)
            self.mx = mx

    def mvm(self,v,out=None):
        # NOTE: If 'out'==None, one temporary vector could be avoided here
        out = self._mvm_checkin(v,out)
        if self.transp:
            mx = self.mx.T
        else:
            mx = self.mx
        # mx.dot creates temp. vector
        out[:] = mx.dot(v)
        return out

    def getcol(self,i,out=None):
        m = self.shape(0)
        out = self._check_resultvec(out,m)
        if self.transp:
            self.mx[i].toarray(out=np.reshape(out,(1,m)))
        else:
            self.mx[:,i].toarray(out=np.reshape(out,(m,1)))
        return out

    def mat_btdb(self,v,out=None):
        if self.transp:
            raise NotImplementedError("NOT IMPLEMENTED")
        mx = self.mx
        if out is None:
            out = np.empty((self.n,self.n))
        else:
            self._matbtdb_checkout(out)
        return mx.T.dot(ssp.diags(v,0,format=
                                  mx.getformat()).dot(mx)).toarray(out=out)

    def diag_bsbt(self,s,out=None):
        if self.transp:
            raise NotImplementedError("NOT IMPLEMENTED")
        out = self._check_resultvec(out,self.m)
        # This solution is not efficient! Two dense m-by-n
        # matrices are created (not just one, because 'multiply' creates
        # a dense return matrix(!)). If n < m, a for loop over 0:(n-1) is not
        # so bad
        tv = np.ones(self.n)
        mx = self.mx
        out[:] = ssp.csr_matrix(mx.multiply(mx.dot(s))).dot(tv)
        return out

class MatDiag(Mat):
    """
    MatDiag
    =======

    MatDiag represents a diagonal matrix B = (diag dg)
    """
    def __init__(self,dg):
        if isinstance(dg,MatDiag):
            Mat.__init__(self,dg)
            self.dg = dg.dg
            self.sqdg = dg.sqdg
        else:
            if not helpers.check_vecsize(dg) or dg.shape[0] == 0:
                raise TypeError('DG wrong type or size')
            n = dg.shape[0]
            Mat.__init__(self,n,n)
            self.dg = dg
            self.sqdg = dg**2

    def mvm(self,v,out=None):
        out = self._mvm_checkin(v,out)
        out[:] = v
        out *= self.dg
        return out

    def getcol(self,i,out=None):
        out = self._check_resultvec(out,self.m)
        out.fill(0.)
        out[i] = self.dg[i]
        return out

    def mat_btdb(self,v,out=None):
        if out is None:
            return np.diag(v*self.sqdg)
        else:
            self._matbtdb_setdgout(out,v*self.sqdg)
            return out

    def diag_bsbt(self,s,out=None):
        out = self._check_resultvec(out,self.m)
        out[:] = self.sqdg
        out *= np.diag(s)
        return out

class MatEye(Mat):
    """
    MatEye
    ======

    MatEye represents an identity matrix
    """
    def __init__(self,m):
        if isinstance(m,MatEye):
            Mat.__init__(self,m)
        else:
            Mat.__init__(self,m,m)

    def mvm(self,v,out=None):
        out = self._mvm_checkin(v,out)
        out[:] = v
        return out

    def getcol(self,i,out=None):
        out = self._check_resultvec(out,self.m)
        out.fill(0.)
        out[i] = 1.
        return out

    def mat_btdb(self,v,out=None):
        if not helpers.check_vecsize(v,self.m):
            raise TypeError('V wrong')
        if out is None:
            return np.diag(v)
        else:
            self._matbtdb_setdgout(out,v)
            return out

    def diag_bsbt(self,s,out=None):
        if (not isinstance(s,np.ndarray) or s.ndim != 2 or
            s.shape[0] != self.m or s.shape[1] != self.m):
            raise TypeError('S wrong')
        out = self._check_resultvec(out,self.m)
        out[:] = np.diag(s)
        return out

class MatSub(Mat):
    """
    MatSub
    ======

    MatSub represents selection operator I_{sind,.}, where sind is a subindex
    of 0:(n-1)
    """
    def __init__(self,n,sind=None):
        if isinstance(n,MatSub):
            Mat.__init__(self,n)
            self.sind = n.sind
        else:
            if sind is None:
                raise ValueError('SIND required')
            Mat.__init__(self,len(sind),n)
            # Test SIND
            tv = np.empty(n)
            try:
                tv2 = tv[sind]
            except IndexError:
                raise IndexError('SIND not a valid subindex')
            self.sind = sind

    def mvm(self,v,out=None):
        out = self._mvm_checkin(v,out)
        if not self.transp:
            out[:] = v[sind]
        else:
            out.fill(0.)
            out[sind] = v
        return out

    def getcol(self,i,out=None):
        m, n = self.shape()
        if not isinstance(i,numbers.Integral) or i<0 or i>=n:
            raise IndexError('I')
        out = self._check_resultvec(out,m)
        if self.transp:
            out.fill(0.)
            out[self.sind[i]] = 1.
        else:
            out[:] = np.array([float(x==i) for x in self.sind])
        return out

    def mat_btdb(self,v,out=None):
        if self.transp:
            raise NotImplementedError("NOT IMPLEMENTED")
        tv = np.zeros(self.n)
        tv[sind] = v
        if out is None:
            return np.diag(tv)
        else:
            self._matbtdb_setdgout(out,tv)
            return out

    def diag_bsbt(self,s,out=None):
        if self.transp:
            raise NotImplementedError("NOT IMPLEMENTED")
        if (not isinstance(s,np.ndarray) or s.ndim != 2 or
            s.shape[0] != self.n or s.shape[1] != self.n):
            raise TypeError('S wrong')
        out = self._check_resultvec(out,self.m)
        out[:] = np.diag(s)[self.sind]
        return out

class MatContainer(Mat):
    """
    MatContainer
    ============

    Container class, represents np.vstack(B1,B2,...), where Bj are Mat objects
    themselves.
    """
    # chpos_k, chpos_j: Row index i falls into child[chpos_k[i]], relative
    # row position chpos_j[i] there
    def __init__(self,child):
        if isinstance(child,MatContainer):
            Mat.__init__(self,child)
            self.child = child.child
            self.chpos_k = child.chpos_k
            self.chpos_j = child.chpos_j
        else:
            if len(child)<2:
                raise TypeError('CHILD must be list of >=2 Mat objects')
            fst = True
            m = 0
            chpos_k = []
            chpos_j = []
            for k in xrange(len(child)):
                chd = child[k]
                if not isinstance(chd,Mat):
                    raise TypeError('CHILD must be list of >=2 Mat objects')
                cm, cn = chd.shape()
                m += cm
                chpos_k.extend(cm*[k])
                chpos_j.extend(range(cm))
                if fst:
                    n = cn
                    fst = False
                elif cn != n:
                    raise ValueError('CHILD entries must have same number of columns')
            Mat.__init__(self,m,n)
            # We do a shallow copy of CHILD, just to make sure
            self.child = list(child)
            self.chpos_k = chpos_k
            self.chpos_j = chpos_j

    def mvm(self,v,out=None):
        out = self._mvm_checkin(v,out)
        if not self.transp:
            off = 0
            for chd in self.child:
                sz = chd.shape(0)
                chd.mvm(v,out[off:off+sz])
                off += sz
        else:
            off = 0
            out.fill(0.)
            tv = np.empty_like(out)
            for chd in self.child:
                sz = chd.shape(0)
                chd.T().mvm(v[off:off+sz],tv)
                out += tv
                off += sz
        return out

    def getcol(self,i,out=None):
        m, n = self.shape()
        if not isinstance(i,numbers.Integral) or i<0 or i>=n:
            raise IndexError('I')
        out = self._check_resultvec(out,m)
        if self.transp:
            self.child[self.chpos_k[i]].T().getcol(self.chpos_j[i],out)
        else:
            off = 0
            for chd in self.child:
                sz = chd.shape(0)
                chd.getcol(i,out[off:off+sz])
                off += sz
        return out

    def mat_btdb(self,v,out=None):
        if self.transp:
            raise NotImplementedError("NOT IMPLEMENTED")
        if not helpers.check_vecsize(v,self.m):
            raise TypeError('V wrong')
        off = 0
        m, n = self.shape()
        if out is None:
            out = np.empty((n,n))
        else:
            self._matbtdb_checkout(out)
        out.fill(0.)
        # 'buffmat1' is allocated once, then reused
        try:
            self.buffmat1.resize((n,n),refcheck=False)
        except AttributeError:
            self.buffmat1 = np.empty((n,n))
        for chd in self.child:
            sz = chd.shape(0)
            chd.mat_btdb(v[off:off+sz],self.buffmat1)
            out += self.buffmat1
            off += sz
        return out

    def diag_bsbt(self,s,out=None):
        if self.transp:
            raise NotImplementedError("NOT IMPLEMENTED")
        m, n = self.shape()
        if not (isinstance(s,np.ndarray) and s.ndim == 2 and
                s.shape[0] == n and s.shape[1] == n):
            raise TypeError('S wrong')
        out = self._check_resultvec(out,m)
        off = 0
        for chd in self.child:
            sz = chd.shape(0)
            chd.diag_bsbt(s,out[off:off+sz])
            off += sz
        return out

class MatFactorizedInf(MatSparse):
    """
    MatFactorizedInf
    ================

    Coupling factor for EP in factorized mode. The 'Mat' methods are
    inherited from 'MatSparse', but not really required.
    'mx' must be of type scipy.sparse.csr_matrix here.
    New services:
    - Internal representation of B: 'rowind', 'colind', 'bvals'. See
      comments in C++ class 'FactorizedEPRepresentation' for details
    - Sparse matrix B**2 in 'b2fact' (required for variance computations)
    """
    def __init__(self,mx):
        if isinstance(mx,MatFactorizedInf):
            MatSparse.__init__(self,mx)
            self.rowind = mx.rowind
            self.colind = mx.colind
            self.bvals = mx.bvals
            self.b2fact = mx.b2fact
        else:
            if not isinstance(mx,ssp.csr_matrix):
                raise TypeError('MX must be scipy.sparse.csr_matrix')
            mx.sort_indices()
            MatSparse.__init__(self,mx)
            m, n = mx.shape
            # The ssp.csr_matrix format is pretty much what we need for
            # 'rowind' and 'bvals'
            self.bvals = mx.data.copy()
            self.rowind = np.empty(mx.nnz+m+1,dtype=np.int32)
            self.rowind[:m+1] = mx.indptr
            self.rowind[m+1:] = mx.indices
            # 'colind': V_i are obtained by doing the same on the transpose.
            # For J_i, we fill a matrix of the same sparsity pattern as 'mx'
            # with 0:(nnz-1), then read out the values of the transpose. In
            # fact, we fill in 1:nnz and subtract 1 later, as otherwise the 0
            # does not count as data entry
            tmpm = ssp.csr_matrix((np.arange(1,mx.nnz+1),mx.indices,
                                   mx.indptr),shape=(m,n),dtype=np.int32)
            tmpm = tmpm.T.tocsr()  # Transpose
            tmpm.sort_indices()
            self.colind = np.empty(2*mx.nnz+n+1,dtype=np.int32)
            self.colind[:n+1] = 2*tmpm.indptr+(n+1)
            off = 0; off2 = n+1
            for i in xrange(1,n+1):
                sz = tmpm.indptr[i]-off
                self.colind[off2:off2+sz] = tmpm.indices[off:off+sz]
                self.colind[off2+sz:off2+2*sz] = tmpm.data[off:off+sz]-1
                off += sz
                off2 += 2*sz
            # B**2 factor into 'b2fact'
            self.b2fact = ssp.csr_matrix((mx.data**2,mx.indices,mx.indptr),
                                         shape=mx.shape)

    def nnz(self):
        return self.mx.getnnz()

    def get_mat(self):
        return self.mx

# Testcode (really basic)

if __name__ == "__main__":
    # __all__ has all Mat subclasses, 1st element is 'Mat'
    for cl in __all__[1:]:
        if cl == 'MatDef':
            m, n = np.random.randint(5,21,2)
            amat = np.random.randn(m,n)
            afct = MatDef(amat)
        elif cl == 'MatSparse':
            m, n = np.random.randint(1000,3001,2)
            spmat = ssp.rand(m,n,0.01).tocsr()
            amat = spmat.toarray()
            afct = MatSparse(spmat)
        elif cl == 'MatDiag':
            m = np.random.randint(5,21); n = m
            dg = np.random.randn(n)
            amat = np.diag(dg)
            afct = MatDiag(dg)
        elif cl == 'MatEye':
            m = np.random.randint(5,21); n = m
            amat = np.eye(n)
            afct = MatEye(n)
        elif cl == 'MatSub':
            n = np.random.randint(5,21)
            m = np.random.randint(2,n)
            sind = np.random.choice(range(n),m,False)
            amat = np.eye(n)[sind]
            afct = MatSub(n,sind)
        else:
            n = np.random.randint(5,21)
            m1 = np.random.randint(2,n)
            sind = np.random.choice(range(n),m1,False)
            a1mat = np.eye(n)[sind]
            a1fct = MatSub(n,sind)
            m = m1
            a2mat = np.eye(n)
            a2fct = MatEye(n)
            m += n
            m3 = np.random.randint(5,21)
            a3mat = np.random.randn(m3,n)
            a3fct = MatDef(a3mat)
            m += m3
            dg = np.random.randn(n)
            a4mat = np.diag(dg)
            a4fct = MatDiag(dg)
            m += n
            m5 = np.random.randint(20,101)
            spmat = ssp.rand(m5,n,0.05).tocsr()
            a5mat = spmat.toarray()
            a5fct = MatSparse(spmat)
            m += m5
            amat = np.vstack((a1mat, a2mat, a3mat, a4mat, a5mat))
            afct = MatContainer([a1fct, a2fct, a3fct, a4fct, a5fct])
        if afct.shape() != (m,n) or afct.T().shape() != (n,m):
            raise IndexError('Mat.shape')
        x = np.random.randn(n)
        y1 = np.dot(amat,x)
        y2 = afct.mvm(x)
        print '%s.mvm: %f' % (cl,helpers.maxreldiff(y1,y2))
        #print y1
        #print y2
        x = np.random.randn(m)
        y1 = np.dot(amat.T,x)
        y2 = afct.T().mvm(x)
        print '%s.T().mvm: %f' % (cl,helpers.maxreldiff(y1,y2))
        i = np.random.randint(0,n)
        y1 = amat[:,i]
        y2 = afct.getcol(i)
        print '%s.getcol: %f' % (cl,helpers.maxreldiff(y1,y2))
        i = np.random.randint(0,m)
        y1 = amat[i]
        y2 = afct.T().getcol(i)
        print '%s.T().getcol: %f' % (cl,helpers.maxreldiff(y1,y2))
        if cl == 'MatSparse':
            print 'm=%d, n=%d' % (m,n)
        dg = np.random.randn(m)
        y1 = np.dot(amat.T,np.dot(np.diag(dg),amat))
        if cl == 'MatSparse':
            print 'Start'
        y2 = afct.mat_btdb(dg)
        print '%s.T().mat_btdb: %f' % (cl,helpers.maxreldiff(y1,y2))
        tmp = np.random.randn(n,n)
        tmp += tmp.T
        y1 = np.diag(np.dot(amat,np.dot(tmp,amat.T)))
        if cl == 'MatSparse':
            print 'Start'
        y2 = afct.diag_bsbt(tmp)
        print '%s.T().diag_bsbt: %f' % (cl,helpers.maxreldiff(y1,y2))
