function ret = test_fact_maxpi(model,repr,flag)
%TEST_FACT_MAXPI Test code for selective damping data structure
%  RET = EPT.TEST_FACT_MAXPI(MODEL,REPR,{FLAG=1})
%  Used to test the max-pi data structure maintained in
%  REPR.SD_XXX. We test whether the maximum values are correct,
%  return RET=1 if not, RET=0 if things are OK.
%  Details: C++ 'FactEPMaximumPiValues', design technical report.

if nargin<3
  flag = 1;
end
[m,n] = size(model.matB);
if ~isfield(repr,'sd_numk') || ~isfield(repr,'sd_numvalid') || ...
      ~isfield(repr,'sd_topind') || ~isfield(repr,'sd_topval')
  error('Need SD_XXX fields in REPR');
end
if isfield(repr,'sd_subind');
  if ~repr.sd_subexcl
    subind = repr.sd_subind+1;
  else
    subind = setdiff((1:m)',repr.sd_subind+1);
  end
else
  subind = (1:m)';
end
numk = repr.sd_numk;
% Test whether maximum values are correct
tmat = reshape(repr.sd_topval,numk+1,n);
mxv1 = tmat(1,:);
tmat = reshape(repr.sd_topind,numk+1,n);
mxi1 = tmat(1,:)+1;
[ri,ci,vl] = find(model.matB'); % Structure of B'
tmpS = sparse(ri,ci,repr.epPi,n,m)';
[mxv2,mxi2] = max(tmpS(subind,:),[],1);
mxv2 = full(mxv2); mxi2 = subind(mxi2)';
if max(abs(mxv1-mxv2))>1e-12 || any(mxi1~=mxi2)
  ret = 1;
  if flag>1
    % DEBUG:
    fprintf(1,'TEST_FACT_MAXPI: Differences:\ntopval:\n');
    [mxv1; mxv2]
    fprintf(1,'topind:\n');
    [mxi1; mxi2]
    pause;
  end
else
  ret = 0;
end
