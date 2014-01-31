function [varargout] = predict_ep(model,repr,imode,ptype,pmodel)
%PREDICT_EP Compute posterior predictions
%  [...] = EPT.PREDICT_EP(MODEL,REPR,IMODE,PTYPE,PMODEL)
%  MODEL, REPR represent model and posterior, IMODE is the
%  inference mode. Compute predictive expectations for test set
%  represented by PMODEL, depending on PTYPE. PMODEL.matB is the
%  coupling factor B_*. For predictive moments, PMODEL.potMan must
%  be given as well (potential manager).
%
%  PTYPE == 1: Gaussian moments.
%  [...] == [H_Q,{RHO_Q}]
%  Means and variances (optional) of Gaussian marginals.
%
%  PTYPE == 2: Predictive moments.
%  [...] == [LOGZ_P,{H_P},{RHO_P}]
%  log normalization constants, means, variances of predictive
%  marginals ("tilted" distributions). If EP computations fail,
%  entries are set to 0, H_Q, RHO_Q resp.
%
%  PTYPE == 3: Gaussian and predictive moments.
%  [...] == [H_Q,RHO_Q,LOGZ_P,H_P,RHO_P]

if nargin<5
  error('Need 5 input arguments');
end
switch imode
 case 'CoupParallel'
  mod = 1;
 case 'CoupSequential'
  mod = 2;
 case 'CoupSeqMargsUp2Date'
  mod = 3;
 case 'Factorized'
  mod = 4;
 otherwise
  error('Unknown IMODE');
end
isCoup = (mod<4);
if ptype~=1 && ptype~=2 && ptype~=3
  error('PTYPE wrong');
end
[pm,n] = size(pmodel.matB);
if n~=size(model.matB,2)
  error('MODEL, PMODEL not compatible');
end
if nargout==0 || (ptype==1 && nargout>2) || ...
	   (ptype==2 && nargout>3) || (ptype==3 && nargout~=5)
  error('Wrong number of return arguments');
end

% Gaussian moments
if isCoup
  % Coupled posterior
  h_q = mvm(pmodel.matB,(repr.lFact')\repr.cVec);
  if ptype~=2
    varargout{1} = h_q;
  end
  if nargout>1 || ptype~=1
    % Gaussian variances: Posterior covariance matrix may already
    % be in REPR
    if mod==1 && isfield(repr,'qCov')
      rho_q = diagBSymBt(pmodel.matB,repr.qCov);
    else
      linv = repr.lFact\eye(size(repr.lFact,1));
      qCov = linv'*linv;
      rho_q = diagBSymBt(pmodel.matB,qCov);
    end
    if ptype~=2
      varargout{2} = rho_q;
    end
  end
else
  % Factorized posterior
  h_q = full(pmodel.matB*(repr.mBeta./repr.mPi));
  rho_q = full((pmodel.matB.^2)*(1./repr.mPi));
  if ptype~=2
    varargout{1} = h_q;
    varargout{2} = rho_q;
  end
end

% Predictive moments
if ptype>1
  % Need potential manager
  if ~isfield(pmodel,'potMan') || ept.potman_size(pmodel.potMan)~= ...
	pm
    error('PMODEL: Potential manager missing or wrong size');
  end
  irp = ept.potman_intrepres(pmodel.potMan);
  [rstat,alpha,nu,logz_p] = ...
      eptools_epupdate_parallel(irp.potids,irp.numpot,irp.parvec, ...
				irp.parshrd,h_q,rho_q);
  indok = find(rstat);
  tvec = 1-nu(indok).*rho_q(indok);
  indok2 = find(tvec>=1e-9);
  if length(indok2)<pm
    indok = indok(indok2);
    indnok = setdiff((1:pm)',indok);
    logz_p(indnok) = 0;
    h_p = h_q; rho_p = rho_q;
    tvec2 = rho_q(indok);
    h_p(indok) = h_q(indok) + alpha(indok).*tvec2;
    rho_p(indok) = tvec2.*tvec(indok2);
  else
    h_p = h_q + alpha.*rho_q;
    rho_p = rho_q.*tvec;
  end
  if ptype==2
    varargout{1} = logz_p;
    if nargout>1
      varargout{2} = h_p;
      if nargout>2
	varargout{3} = rho_p;
      end
    end
  else
    % PTYPE==3: H_Q, RHO_Q already assigned
    varargout(3:5) = {logz_p,h_p,rho_p};
  end
end
