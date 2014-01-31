function S = matBtDgB(B,v)
%Mat.matBtDgB: Compute B'*diag(v)*B
%  S = matBtDgB(B,V)  computes S = B'*diag(V)*B
%
%  NOTE: Components of V may be negative.
%  NOTE: This method need not be implemented for B.TRANSP==true,
%  this service is not required.

% Default implementation by MVM
[m,n] = size(B);
if numel(v)~=m
  error('V has wrong size');
end
S = zeros(n,n);
for i=1:n
  S(:,i) = mvm(ctranspose(B),v(:).*getCol(B,i));
end
