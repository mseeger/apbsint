function S = matBtDgB(B,v)

if B.transp
  % This sucks: Not required anyway
  error('NOT IMPLEMENTED!');
end
[m,n] = size(B);
if numel(v)~=m
  error('V has wrong size');
end
S = zeros(n,n);
off = 0;
for k = 1:numel(B.child)
  sz = size(B.child{k},1);
  S = S + matBtDgB(B.child{k},v(off+1:off+sz));
  off = off+sz;
end
