function v = diagBSymBt(B,S)

if B.transp
  % This sucks: Not required anyway
  error('NOT IMPLEMENTED!');
end
[m,n] = size(B);
if size(S,1)~=n || size(S,2)~=n
  error('M has wrong size');
end
v = zeros(m,1);
off = 0;
for k = 1:numel(B.child)
  sz = size(B.child{k},1);
  v(off+1:off+sz) = diagBSymBt(B.child{k},S);
  off = off+sz;
end
