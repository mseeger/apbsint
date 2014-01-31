function v = getCol(B,i)

[m,n] = size(B);
if i<=0 || i>n
  error('I out of range');
end
if ~B.transp
  v = zeros(m,1);
  off = 0;
  for k = 1:numel(B.child)
    sz = size(B.child{k},1);
    v(off+1:off+sz) = getCol(B.child{k},i);
    off = off+sz;
  end
else
  [k,j] = B.getRelPos(i);
  v = getCol(ctranspose(B.child{k}),j);
end
