function repr = refresh_repres(model,repr,mode,st_cov,deb_fid)
%REFRESH_REPRES Recompute representation from EP parameters
%  REPR = EPT.REFRESH_REPRES(MODEL,REPR,MODE,{ST_COV=1})
%  Recomputes representation from EP parameters.
%  If MODE == 'CoupParallel', the posterior covariance matrix is
%  stored in REPR.qCov. This is not done if ST_COV is false.

global debug_bmat; % DEBUG!
do_debug = 0;

if nargin<4
  st_cov = 1;
end
[m,n] = size(model.matB);
st_cov = st_cov & strcmp(mode,'CoupParallel');
if ~strcmp(mode,'Factorized')
  if isempty(debug_bmat)
    if ~st_cov
      amat = matBtDgB(model.matB,repr.epPi);
      amat = 0.5*(amat+amat');
      repr.lFact = chol(amat,'lower');
      %deb_amat = amat; % DEBUG!
    else
      repr.qCov = matBtDgB(model.matB,repr.epPi);
      repr.qCov = 0.5*(repr.qCov+repr.qCov');
      repr.lFact = chol(repr.qCov,'lower');
      %deb_amat = repr.qCov; % DEBUG!
    end
    repr.cVec = repr.lFact\mvm(model.matB',repr.epBeta);
    %deb_lfact = repr.lFact; % DEBUG
    %deb_cvec = repr.cVec; % DEBUG
  else
    % DEBUG!
    amat = full(debug_bmat'*(repmat(repr.epPi(:),1,n).*debug_bmat));
    amat = (amat+amat')/2;
    repr.lFact = chol(amat,'lower');
    repr.cVec = repr.lFact\full(debug_bmat'*repr.epBeta);
    % END DEBUG
  end
  if ~strcmp(mode,'CoupSequential')
    % Compute marginal moments
    if isempty(debug_bmat)
      repr.mMeans = mvm(model.matB,repr.lFact'\repr.cVec);
    else
      repr.mMeans = full(debug_bmat*(repr.lFact'\repr.cVec));
    end
    %deb_mmeans = repr.mMeans; % DEBUG
    % Compute inverse from Cholesky factor
    % TODO: This may be a bottleneck! If so, check whether LAPACK
    % DPOTRI is faster
    if ~do_debug
      linv = repr.lFact\eye(n);
      if ~st_cov
	amat = linv'*linv;
	amat = 0.5*(amat+amat');
	%deb_inva = amat; % DEBUG
	if isempty(debug_bmat)
	  repr.mVars = diagBSymBt(model.matB,amat);
	else
	  repr.mVars = full(sum(debug_bmat.*(debug_bmat*amat),2));
	end
      else
	repr.qCov = linv'*linv; % Posterior covariance
	repr.qCov = 0.5*(repr.qCov+repr.qCov');
	%deb_inva = repr.qCov; % DEBUG
	if isempty(debug_bmat)
	  repr.mVars = diagBSymBt(model.matB,repr.qCov);
	else
	  repr.mVars = full(sum(debug_bmat.*(debug_bmat*repr.qCov),2));
	end
      end
      % DEBUG:
      %deb_mvars = repr.mVars;
      %save('deb_vars','deb_amat','deb_lfact','deb_cvec', ...
      %   'deb_mmeans','deb_inva','deb_mvars');
      %fprintf(1,'DEBUG(refresh_repres): Stored stuff\n');
      %pause;
      % END DEBUG
    else
      % DEBUG: Code from glm-ie
      ei = zeros(n,1);
      z = zeros(m,1);
      for i=1:n
	ei(i) = 1;
	w = (repr.lFact')\ei;
	if isempty(debug_bmat)
	  v = mvm(model.matB,w);
	else
	  v = debug_bmat*w;
	end
	z = z+(v.^2);
	ei(i) = 0;
      end
      repr.mVars = z;
      % END DEBUG
    end
  end
else
  % Factorized mode
  % NOTE: EPTOOLS_FACT_COMPMARGINALS overwrites input argument
  % (undocumented). We use copies here to avoid issues with Matlab
  mBeta = zeros(n,1); mPi = zeros(n,1);
  if isfield(model,'matBInt')
    irb = model.matBInt;
  else
    irb = ept.bfact_intrepres(model.matB);
  end
  eptools_fact_compmarginals(n,m,irb.rowind,irb.colind,irb.bvals, ...
			     repr.epPi,repr.epBeta,mPi,mBeta);
  repr.mBeta = mBeta(:); repr.mPi = mPi(:); % Write-back
  if nargin>4
    ept.debug_printvec(int32([n m]),deb_fid);
    ept.debug_printvec(irb.rowind,deb_fid);
    ept.debug_printvec(irb.colind,deb_fid);
    ept.debug_printvec(irb.bvals,deb_fid);
    ept.debug_printvec(repr.epPi,deb_fid);
    ept.debug_printvec(repr.epBeta,deb_fid);
    ept.debug_printvec(mPi,deb_fid);
    ept.debug_printvec(mBeta,deb_fid);
  end
end
