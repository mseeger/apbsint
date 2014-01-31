function [model,repr,delta,skip] = ...
    inf_coup_sequential(model,repr,updind,imode,damp,caveps, ...
			skipeps,deb_fid)
%INF_COUP_SEQUENTIAL EP inference (coupled model, sequential updates)
%  [MODEL,REPR,{DELTA},{SKIP}] =
%    EPT.INF_COUP_SEQUENTIAL(MODEL,REPR,UPDIND,IMODE,{DAMP=0},
%                            {CAVEPS=1e-5},{SKIPEPS=1e-8})
%  EP inference, sequential mode: Run EP updates on potentials in
%  UPDIND. IMODE is 'CoupSequential' or 'CoupSeqMargsUp2Date'. For
%  the latter, REPR contains marginal moments which are kept
%  up-2-date. DAMP is damping factor (0: No damping).
%
%  Update on k is skipped if cavVar(k)/mVar(k) > 1/CAVEPS, if
%  local EP update fails for numerical reasons, or if the absolute
%  change |Delta pi_k| is smaller than SKIPEPS.
%  NOTE: The last condition is not counted as a skip (only if it
%  comes about due to selective damping).
%  We apply selective damping to ensure that
%    1 + (Delta pi_k) mVar(k) >= CAVEPS
%  DELTA, SKIP are vectors same size as UPDIND. DELTA contains
%  relative change in marginal (only at potential k), or 0 for
%  skipped updates. SKIP is vector for skipped updates, entries:
%  - 0: Not skipped
%  - 1: Skipped due to cavity marginal or local EP failure
%  - 2: Skipped due to small change, selective damping (would not
%       be skipped without that)
%  - 3: Skipped due to Cholesky up/downdate error
%
%  TODO:
%  - Selective damping/skipping in order to keep cavity marginals
%    well-defined (in mode 'CoupSeqMargsUp2Date')

if nargin<4
  error('Not enough input arguments');
end
doU2D = strcmp(imode,'CoupSeqMargsUp2Date');
if ~doU2D && ~strcmp(imode,'CoupSequential')
  error('IMODE wrong');
end
[m,n] = size(model.matB);
if ~isempty(find(~ismember(updind,(1:m)')))
  error('UPDIND has invalid entries');
end
if nargin<7
  skipeps = 1e-8;
  if nargin<6
    caveps = 1e-5;
    if nargin<5
      damp=0;
    end
  end
end
% Check, and compute internal representation
model = ept.check_model_repres(model,repr,imode);
irp = model.potManInt;

% Buffers for Cholesky up/downdates
cup_c = zeros(n,1); cup_s = zeros(n,1);
cup_wk = zeros(n,1); cup_z = zeros(1,n);
% Loop over potentials to update
nupd = length(updind);
skip = zeros(nupd,1);
delta = zeros(nupd,1);
for k = 1:nupd
  j = updind(k);
  epPi = repr.epPi(j); epBeta = repr.epBeta(j);
  doSkip = 0;
  % Compute marginal moments. This is done also in "up-2-date" mode,
  % since VVEC is required anyway. We refresh the moments in this
  % case.
  bvec = getCol(model.matB',j);
  vvec = repr.lFact\bvec;
  mu = vvec'*repr.cVec;
  rho = vvec'*vvec;
  if doU2D
    % Refresh
    if ept.reldiff(repr.mMeans(j),mu)>1e-3
      fprintf(1,['UUPS(inf_coup_sequential): j=%d, ' ...
		 'mMeans=%f, mu=%f\n'],j,repr.mMeans(j),mu);
    end
    if ept.reldiff(repr.mVars(j),rho)>1e-3
      fprintf(1,['UUPS(inf_coup_sequential): j=%d, ' ...
		 'mVars=%f, rho=%f\n'],j,repr.mVars(j),rho);
    end
    repr.mMeans(j) = mu;
    repr.mVars(j) = rho;
  end
  tscal = 1-epPi*rho;
  if tscal >= caveps
    % Cavity marginal
    crho = rho/tscal;
    cmu = (mu-epBeta*rho)/tscal;
    % Local EP update (NOTE: 1-floor -> 0-floor)
    [rstat,alpha,nu] = ...
	eptools_epupdate_single(irp.potids,irp.numpot,irp.parvec, ...
				irp.parshrd,j-1,cmu,crho);
    if nargin>7
      % DEBUG: Print input and return args
      ept.debug_printvec(irp.potids,deb_fid);
      ept.debug_printvec(irp.numpot,deb_fid);
      ept.debug_printvec(irp.parvec,deb_fid);
      ept.debug_printvec(irp.parshrd,deb_fid);
      ept.debug_printvec(int32(j-1),deb_fid);
      ept.debug_printvec(cmu,deb_fid);
      ept.debug_printvec(crho,deb_fid);
      ept.debug_printvec(int32(rstat),deb_fid);
      ept.debug_printvec(alpha,deb_fid);
      ept.debug_printvec(nu,deb_fid); 
    end
    if ~rstat
      % Skip (local EP update fails)
      doSkip = 1;
    else
      tscal = 1-nu*crho;
      if tscal>=1e-7
	newPi = nu/tscal;
	newBeta = (cmu*nu+alpha)/tscal;
      else
	% Skip (local EP update fails)
	doSkip = 1;
      end
    end
  else
    % Skip (cavity marginal invalid)
    doSkip = 1;
  end
  if ~doSkip
    % Damping
    dflPi = newPi-epPi; % Full update
    dflBeta = newBeta-epBeta;
    delPi = (1-damp)*dflPi;
    delBeta = (1-damp)*dflBeta;
    delPi2 = delPi; delBeta2 = delBeta; % For 'doSkip' decision
    % Selective damping
    if delPi*rho+1<caveps
      delPi = (caveps-1)/rho;
      delBeta = (delPi/dflPi)*delBeta;
      % DEBUG
      fprintf(1,'Sel.damp: j=%d, damp=%f\n',j,1-delPi/dflPi);
    end
    newPi = epPi+delPi; newBeta = epBeta+delBeta;
    if abs(delPi)>=skipeps
      % Update representation
      if doU2D
	% w = B L^-T v, required below
	wvec = mvm(model.matB,repr.lFact'\vvec);
      end
      if delPi>0
	% Cholesky update
	tscal = sqrt(delPi);
	bvec = tscal*bvec;
	yscal = delBeta/tscal;
	% NOTE: (:) should force copy here: MEX function simply
        % overwrites CUP_Z!
	cup_z = (repr.cVec(:))';
	% L matrix: SCODE is 'LN', where L is for lower triangular
	rstat = eptools_choluprk1({repr.lFact,[1 1 n n],'LN'},bvec, ...
				  cup_c,cup_s,cup_wk,cup_z,yscal);
	if ~rstat
	  repr.cVec = cup_z'; % New c vector
	else
	  fprintf(1,['UUPS: Cholesky update failed. Representation' ...
		     ' may be corrupt!']);
	  doSkip = 3;
	end
      else
	% Cholesky downdate
	% NOTE: delPi*rho+1 >= caveps by selective damping, so the
        % downdate should work out
	tscal = sqrt(-delPi);
	vvec = tscal*vvec; % L^-1 b
	yscal = -delBeta/tscal;
	% NOTE: (:) should force copy here: MEX function simply
        % overwrites CUP_Z!
	cup_z = (repr.cVec(:))';
	% L matrix: SCODE is 'LN', where L is for lower triangular
	rstat = eptools_choldnrk1({repr.lFact,[1 1 n n],'LN'},vvec, ...
				  cup_c,cup_s,cup_wk,1,cup_z, ...
				  yscal);
	if ~rstat
	  repr.cVec = cup_z'; % New c vector
	else
	  fprintf(1,['UUPS: Cholesky downdate failed. Representation' ...
		     ' may be corrupt!']);
	  doSkip = 3;
	end
      end
      if doU2D && ~doSkip
	% Update marginal moments
	tscal=1/(delPi*rho+1);
	repr.mMeans = repr.mMeans + ((delBeta-delPi*mu)*tscal)* ...
	    wvec;
	repr.mVars = repr.mVars - (delPi*tscal)*(wvec.^2);
      end
    else
      % Skip (small change)
      if abs(delPi2)>=skipeps
	doSkip = 2; % Small change due to selective damping
      else
	doSkip = 4; % Not counted as a skip
      end
    end
    if doSkip<4
      skip(k) = doSkip;
    else
      skip(k) = 0; % 'doSkip'==4 not counted as skip
    end
    if ~doSkip
      repr.epPi(j) = newPi;
      repr.epBeta(j) = newBeta;
      % DELTA depends on new marginal moments
      hrho = crho*(1-nu*crho);
      hmu = cmu+alpha*crho;
      delta(k) = max(ept.reldiff([hmu sqrt(hrho)],[mu sqrt(rho)]));
    end
  end
end
