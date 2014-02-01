function bv = mvm(B,v)
%MatContainer.mvm: Matrix-vector product
%  bv = mvm(B,v)  computes MVM bv = B*v.

[m,n] = size(B);
if numel(v)~=n
  error('V has wrong size');
end
bv = zeros(m,1);
if ~B.transp
  off = 0;
  for k = 1:numel(B.child)
    sz = size(B.child{k},1);
    bv(off+1:off+sz) = mvm(B.child{k},v);
    off = off+sz;
  end
else
  off = 0;
  for k = 1:numel(B.child)
    sz = size(B.child{k},1);
    bv = bv + mvm(ctranspose(B.child{k}),v(off+1:off+sz));
    off = off+sz;
  end
end
