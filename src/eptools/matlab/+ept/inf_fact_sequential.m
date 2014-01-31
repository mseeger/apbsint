function [model,repr,delta,skip,sd_damp] = ...
    inf_fact_sequential(model,repr,updind,damp,piminthres,deb_fid)
%INF_FACT_SEQUENTIAL EP inference (factorized backbone)
%  [MODEL,REPR,{DELTA},{SKIP},{SD_DAMP}] =
%    EPT.INF_FACT_SEQUENTIAL(MODEL,REPR,UPDIND,{DAMP=0},
%                            {PIMINTHRES=1e-8})
%  EP inference, factorized backbone distribution: EP updates are
%  run sequentially on potentials in UPDIND. DAMP is damping factor
%  (0: No damping). PIMINTHRES: See EPTOOLS_FACT_SEQUPDATES.
%
%  An update can be skipped for a number of reasons. The skip
%  status for each update is returned in SKIP (0 means no skip; see
%  EPTOOLS_FACT_SEQUPDATES).
%  DELTA contains relative change in marginal moments (max over
%  nonzeros) for each update, or 0 for skipped updates.
%  SD_DAMP can be used only if the selective damping mechanism is
%  active (fields SD_XXX in REPR). In this case, it returns the
%  effective damping factor used for each update (entries reliable
%  only for entries which are not skipped).

if nargin<3
  error('Not enough input arguments');
end
[m,n] = size(model.matB);
imode = 'Factorized';
if ~isempty(find(~ismember(updind,(1:m)')))
  error('UPDIND has invalid entries');
end
if nargin<5
  piminthres = 1e-8;
  if nargin<4
    damp=0;
  end
end
% Check, and compute internal representations
model = ept.check_model_repres(model,repr,imode);
irp = model.potManInt;
irb = model.matBInt;

% Loop over potentials to update
% Note that EPTOOLS_FACT_SEQUPDATES overwrites the buffers of some
% input arguments (undocumented and risky, but we cannot allocate
% lots of return argument memory for each of these
% calls). Therefore, we work on copies inside the loop, and force a
% real copy by using (:).
% ATTENTION: This could be expensive for large models. If so, try
% what happens if no copies are drawn.
epBeta = repr.epBeta(:);
epPi = repr.epPi(:);
mBeta = repr.mBeta(:);
mPi = repr.mPi(:);
if isfield(repr,'sd_numk')
  if ~isfield(repr,'sd_numvalid') || ~isfield(repr,'sd_topind') || ...
	~isfield(repr,'sd_topval')
    error('REPR must have SD_XXX fields');
  end
  sd_numvalid = repr.sd_numvalid(:);
  sd_topind = repr.sd_topind(:);
  sd_topval = repr.sd_topval(:);
end
if isa(updind,'int32')
  updind0 = updind-1; % 0-floor
else
  updind0 = int32(updind-1);
end
% Everything is done by EPTOOLS_FACT_SEQUPDATES
%fprintf(1,'Calling MEX\n'); % DEBUG
args = {n,m,updind0,irp.potids,irp.numpot,irp.parvec,irp.parshrd, ...
	irb.rowind,irb.colind,irb.bvals,epPi,epBeta,mPi,mBeta, ...
	piminthres,damp};
if isfield(repr,'sd_numk')
  args(end+1:end+3) = {sd_numvalid,sd_topind,sd_topval};
  if isfield(repr,'sd_subind')
    args{end+1} = repr.sd_subind;
    if isfield(repr,'sd_subexcl')
      args{end+1} = repr.sd_subexcl;
    end
  end
end
if nargin>5
  % DEBUG: Print input (and I/O) args
  ept.debug_printvec(int32([n m]),deb_fid);
  for k=3:length(args)
    ept.debug_printvec(args{k},deb_fid);
  end
end
if nargout>2
  if isfield(repr,'sd_numk')
    % DEBUG
    %fprintf(1,['n=%d,K=%d,n(K+1)=%d,numvalid=%d,topind=%d,' ...
    %       'topval=%d\n'],n,repr.sd_numk,n*(repr.sd_numk+1), ...
    %    numel(sd_numvalid),numel(sd_topind),numel(sd_topval));
    [skip,delta,sd_damp] = eptools_fact_sequpdates(args{:});
  else
    [skip,delta] = eptools_fact_sequpdates(args{:});
  end
else
  eptools_fact_sequpdates(args{:});
end
if nargin>5
  % DEBUG: Print return and I/O args
  ept.debug_printvec(epPi,deb_fid);
  ept.debug_printvec(epBeta,deb_fid);
  ept.debug_printvec(mPi,deb_fid);
  ept.debug_printvec(mBeta,deb_fid);
  if isfield(repr,'sd_numk')
    ept.debug_printvec(sd_numvalid,deb_fid);
    ept.debug_printvec(sd_topind,deb_fid);
    ept.debug_printvec(sd_topval,deb_fid);
  end
  if nargout>2
    ept.debug_printvec(skip,deb_fid);
    ept.debug_printvec(delta,deb_fid);
    if isfield(repr,'sd_numk')
      ept.debug_printvec(sd_damp,deb_fid);
    end
  end
end
%fprintf(1,'Done\n'); % DEBUG
% Write back
repr.epBeta = epBeta(:);
repr.epPi = epPi(:);
repr.mBeta = mBeta(:);
repr.mPi = mPi(:);
if isfield(repr,'sd_numk')
  repr.sd_numvalid = sd_numvalid(:);
  repr.sd_topind = sd_topind(:);
  repr.sd_topval = sd_topval(:);
end
