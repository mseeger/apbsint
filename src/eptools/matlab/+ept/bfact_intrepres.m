function irp = bfact_intrepres(bmat)
%BFACT_INTREPRES Internal representation for sparse B matrix
%  IRP = EPT.BFACT_INTREPRES(BMAT)
%  Creates internal representation IRP from sparse matrix BMAT.
%  This is used for MEX transfer.

[m,n] = size(bmat);
if m==0 || n==0 || ~issparse(bmat)
  error('BMAT wrong');
end

% BVALS and ROWIND
bT = bmat'; % Much simpler to traverse
irp.bvals = nonzeros(bT);
irp.rowind = int32(zeros(nnz(bmat)+m+1,1));
% ATTENTION: Indexes must be 0-floor
off = 0;
for j = 1:m
  irp.rowind(j) = int32(off);
  ri = find(bT(:,j)); sz=numel(ri);
  if sz==0
    error(sprintf('Row %d of BMAT is zero',j));
  end
  off2 = off+m+1;
  irp.rowind(off2+1:off2+sz) = int32(ri(:)-1);
  off = off+sz;
end
if off~=nnz(bmat)
  error('Internal error!');
end
irp.rowind(m+1) = off;

% COLIND. For each column i, we need J_i, the position of nonzero
% elements in the flat array BVALS. This is done by building a
% matrix with the same sparsity pattern as BMAT, but filling the
% nonzeros with [0 1 2 ...] in row major ordering
[bi,bj] = find(bT);
tmpT = sparse(bi,bj,(1:numel(bi))',n,m);
% Transpose does reordering (because NONZEROS reads out in
% column-major ordering)
jivec = nonzeros(tmpT')-1;
irp.colind = int32(zeros(2*nnz(bmat)+n+1,1));
off = n+1; offj = 0;
for i = 1:n
  irp.colind(i) = int32(off);
  ri = find(bmat(:,i)); sz=numel(ri);
  if sz>0
    irp.colind(off+1:off+sz) = int32(ri(:)-1);
    irp.colind(off+sz+1:off+2*sz) = int32(jivec(offj+1:offj+sz));
  else
    % DEBUG
    fprintf(1,'Warning: Column %d of BMAT is zero\n',i);
  end
  off = off+2*sz; offj = offj+sz;
end
if off~=2*nnz(bmat)+n+1 || offj~=nnz(bmat)
  error('Internal error!');
end
irp.colind(n+1) = off;
