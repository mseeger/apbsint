function d = reldiff(a,b,eps)
%RELDIFF Computes relative difference
%  D = RELDIFF(A,B,{EPS=1e-8}). A, B must be matrices of the same
%  total number of entries. D is a column vector.
%  NOTE: If A or B is a row vector, it is converted into a column
%  vector first.

if nargin<3
  eps = 1e-8;
end
d = abs(a(:)-b(:))./max(max(abs([a(:) b(:)]),[],2),eps);
