function [model,repr,delta,nskip] = ...
    inf_coup_parallel(model,repr,damp,caveps,deb_fid)
%INF_COUP_PARALLEL EP inference (coupled model, parallel updates)
%  [MODEL,REPR,{DELTA},{NSKIP}] =
%    EPT.INF_COUP_PARALLEL(MODEL,REPR,{DAMP=0},{CAVEPS=1e-5})
%  EP inference, coupled mode: Run one round of parallel updating.
%  DAMP is damping factor (0: No damping).
%  Update on k is skipped if cavVar(k)/mVar(k) > 1/CAVEPS, or if
%  local EP update fails for numerical reasons. NSKIP returns
%  number of skipped updates. Potentials of type 'Gaussian' are not
%  updated on, but they are not counted as skipped either.
%  DELTA returns measure of relative change of marginal moments
%  (means, stddev's), can be used for deciding EP convergence.
%
%  TODO:
%  - Selective damping/skipping in order to keep cavity marginals
%    well-defined.
%  - Recovery from error when recomputing the representation

global debug_pydeb debug_pydeb_iter;

do_debug = 0; % DEBUG

if nargin<4
  caveps = 1e-5;
  if nargin<3
    damp=0;
  end
end
infMode = 'CoupParallel';
% Check, and compute internal representation
model = ept.check_model_repres(model,repr,infMode);
irp = model.potManInt;
updind = int32(model.potInd-1); % Non-Gaussian potentials (0-floor)
[m,n] = size(model.matB);
% Local EP updates
% We update only on potentials in MODEL.potInd (excludes
% Gaussians), but those not in MODEL.potInd are not counted as
% skipped
indok = model.potInd; mm = length(indok);
if ~do_debug
  cmu = repr.mMeans(indok); crho = repr.mVars(indok);
  tvec = 1-repr.epPi(indok).*crho;
  indok2 = find(tvec>=caveps);
  % DEBUG:
  if numel(indok2)<mm
    indnok = setdiff((1:mm)',indok2);
    fprintf(1,'Skip %d due to undef. cavity\n',numel(indnok));
    for i = indnok'
      j = indok(i);
      fprintf(1,'%d: rho=%f,pi=%f,den=%f,h=%f\n',j,crho(i), ...
	      repr.epPi(j),tvec(i),cmu(i));
    end
  end
  % END DEBUG
  tvec = 1./tvec(indok2);
  indok = indok(indok2);
  cmu(indok2) = (cmu(indok2)-repr.epBeta(indok).*crho(indok2)).*tvec;
  crho(indok2) = crho(indok2).*tvec;
  if debug_pydeb
    deb_cmu = cmu; deb_crho = crho;
  end
else
  % DEBUG: Compute cavity margs as in glm-ie
  cpi = 1./repr.mVars(indok)-repr.epPi(indok);
  cbeta = repr.mMeans(indok)./repr.mVars(indok)-repr.epBeta(indok);
  crho = 1./cpi;
  cmu = cbeta./cpi;
  % END DEBUG
end
if ~do_debug
  [rstat,alpha,nu] = eptools_epupdate_parallel(irp.potids,irp.numpot, ...
					       irp.parvec, ...
					       irp.parshrd,cmu,crho, ...
					       updind);
  if nargin>4
    % DEBUG: Store input/return
    ept.debug_printvec(irp.potids,deb_fid);
    ept.debug_printvec(irp.numpot,deb_fid);
    ept.debug_printvec(irp.parvec,deb_fid);
    ept.debug_printvec(irp.parshrd,deb_fid);
    ept.debug_printvec(cmu,deb_fid);
    ept.debug_printvec(crho,deb_fid);
    ept.debug_printvec(rstat,deb_fid);
    ept.debug_printvec(alpha,deb_fid);
    ept.debug_printvec(nu,deb_fid);
  end
  if debug_pydeb
    deb_alpha = alpha; deb_nu = nu;
  end
else
  % DEBUG: Local EP updates using GPML code
  % ATTENTION, specific to EPTEST_BINCLASS_COUP with Gaussian
  % prior!
  if numel(updind)~=model.potMan{2}.size
    error('DEBUG ERROR');
  end
  [logz,alpha,mnu] = likErf([],model.potMan{2}.pars{1},cmu,crho, ...
			    'infEP');
  nu = -mnu;
  rstat = true(size(alpha));
  % END DEBUG
end
if ~do_debug
  % DEBUG:
  indnok = find(~rstat);
  if ~isempty(indnok)
    fprintf(1,'Skip %d due to EP update error\n',numel(indnok));
    for i = indnok'
      j = indok(i);
      fprintf(1,'%d: cmu=%f,crho=%f\n',j,cmu(i),crho(i));
    end
  end
  % END DEBUG
  indnok = setdiff((1:mm)',indok2); % Complement
  rstat(indnok) = 0; % Filter out undef. cavity positions
  indok2 = find(rstat); % Positions to update on
  newPi = repr.epPi; newBeta = repr.epBeta;
  nu = nu(indok2); alpha = alpha(indok2);
  tvec = 1-nu.*crho(indok2);
  % Just to be sure. These should not fire for positions where EP
  % updates went fine
  indok3 = find(tvec>=1e-7);
  if length(indok3)<length(indok2)
    fprintf(1,'UUPS(inf_coup_parallel): On %d\n',length(indok2)- ...
	    length(indok3));
    tvec = tvec(indok3);
    nu = nu(indok3); alpha = alpha(indok3);
    indok2 = indok2(indok3);
  end
  indok = model.potInd(indok2);
  tvec = 1./tvec;
  newPi(indok) = nu.*tvec;
  newBeta(indok) = (cmu(indok2).*nu+alpha).*tvec;
else
  % DEBUG: Compute new EP pars as in glm-ie
  newPi = repr.epPi; newBeta = repr.epBeta;
  newPi(indok) = nu./(1-nu./cpi);
  newBeta(indok) = (alpha+cmu.*nu)./(1-nu./cpi);
  %END DEBUG
end
% Damping
if damp>0
  newPi(indok) = (1-damp)*newPi(indok) + damp*repr.epPi(indok);
  newBeta(indok) = (1-damp)*newBeta(indok) + damp* ...
      repr.epBeta(indok);
end
nskip = mm-length(indok); % Number of skips
% Recompute representation
indok = model.potInd;
if nargout>2
  oldMargs = [repr.mMeans(indok); sqrt(repr.mVars(indok))];
end
repr.epPi = newPi; repr.epBeta = newBeta;
repr = ept.refresh_repres(model,repr,infMode);
if nargout>2
  delta = max(ept.reldiff(oldMargs,[repr.mMeans(indok); ...
		    sqrt(repr.mVars(indok))]));
end
if debug_pydeb
  cmu = deb_cmu; crho = deb_crho;
  alpha = deb_alpha; nu = deb_nu;
  ep_pi = repr.epPi; ep_beta = repr.epBeta;
  lfact = repr.lFact; cvec = repr.cVec;
  mmeans = repr.mMeans; mvars = repr.mVars;
  save(sprintf('eptb9_vars_it%d.mat',debug_pydeb_iter),'cmu', ...
       'crho','alpha','nu','ep_pi','ep_beta','lfact','cvec', ...
       'mmeans','mvars','-v7');
end
