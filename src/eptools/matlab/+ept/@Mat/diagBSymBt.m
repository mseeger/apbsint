function v = diagBSymBt(B,S)
%Mat.diagBSymBt: Compute diag(B*S*B') for symmetric S
%  v = diagBSymBt(B,S)  computes v = diag(B*S*B') for a symmetric
%  matrix S.
%
%  NOTE: This method need not be implemented for B.TRANSP==true,
%  this service is not required.

% Default implementation by MVM (does not use fact that S is
% symmetric)
[m,n] = size(B);
if size(S,1)~=n || size(S,2)~=n
  error('M has wrong size');
end
v = zeros(m,1);
for i=1:n
  v = v + mvm(B,S(:,i)).*getCol(B,i);
end
