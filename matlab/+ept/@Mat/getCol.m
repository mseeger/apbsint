function v = getCol(B,i)
%Mat.getCol: Return column of B factor
%  v = getCol(B,i)  returns i-th column of B

% Default implementation by MVM
sz = size(B,2);
if i<=0 || i>sz
  error('I out of range');
end
dv = zeros(sz,1); dv(i)=1;
v = mvm(B,dv);
