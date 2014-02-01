function [model,repr,res,res_det] = inf_ep(model,repr,imode,opts)
%INF_EP Gateway function for EP inference
%  [MODEL,REPR,{RES},{RES_DET}] = EPT.INF_EP(MODEL,REPR,IMODE,OPTS)
%  Runs expectation propagation (EP) inference. MODEL describes the
%  model, REPR is the representation of the posterior. Inference
%  operates on REPR, depending on the mode IMODE and options in
%  OPTS. The MODEL entries remain constant, with the exception of
%  internal representations (which is why MODEL is returned). Basic
%  return arguments are fields in RES, more details can be
%  returned in RES_DET. OPTS, RES, RES_DET depend on IMODE.
%  General OPTS fields:
%  - VERBOSE: Verbosity level. 0: No output
%  - TEST_MODEL: Optional. Model for test set. If given, the test
%    set accuracy and log likelihood is computed and output after
%    each sweep.
%    NOTE: We require that the test set target vector (values
%    -1,+1) is given by the concatenation of
%      opts.test_model.potMan{*}.pars{1}, * = 1,2,...
%    This is the case for binary classification potentials.
%
%  Coupled posterior, parallel updating (IMODE == 'CoupParallel'):
%  Calls internal EPT.INF_COUP_PARALLEL.
%  OPTS fields:
%  - MAXIT: Maximum number of sweeps (over non-Gaussian potentials)
%  - DELTAEPS: Converged if DELTA < DELTAEPS, return of
%    EPT.INF_COUP_PARALLEL
%  - DAMP: Damping constant. Def.: 0
%  - CAVEPS: See EPT.INF_COUP_PARALLEL. Def.: 1e-5
%  RES fields:
%  - RSTAT: Return status
%    - 0: Converged to accuracy DELTAEPS
%    - 1: Done MAXIT sweeps
%  - NIT:   Number of sweeps done
%  - DELTA: Return from final EPT.INF_COUP_PARALLEL call
%  - NSKIP: Total number of skipped updates over all sweeps
%  RES_DET fields:
%  - DELTA: Return EPT.INF_COUP_PARALLEL for every sweep
%  - NSKIP: Number of skips for every sweep
%
%  Coupled posterior, sequential updating (IMODE ==
%  'CoupSequential' or 'CoupSeqMargsUp2Date'):
%  Calls internal EPT.INF_COUP_SEQUENTIAL.
%  For IMODE == 'CoupSeqMargsUp2Date', all marginals are kept
%  up-2-date at all times. Right now, this is not really used, so a
%  waste of time.
%  TODO: Implement optimized scheduling, by forward scoring of EP
%  updates.
%  We run sweeps over all potentials of type != 'Gaussian'. For
%  each sweep, we pick a random partition of these potentials.
%  OPTS fields:
%  - MAXIT: Maximum number of sweeps (over non-Gaussian
%    potentials)
%  - DELTAEPS: Converged if max of DELTA < DELTAEPS, return of
%    EPT.INF_COUP_SEQUENTIAL
%  - DAMP: Damping constant. Def.: 0
%  - CAVEPS: See EPT.INF_COUP_SEQUENTIAL. Def.: 1e-5
%  - SKIPEPS: See EPT.INF_COUP_SEQUENTIAL. Def.: 1e-8
%  - REFRESH: If true, the representation is recomputed from
%    scratch after each sweep. Def.: true
%  - UPD_1STSWEEP: Cell vector of potential type names. First sweep
%    updates on potentials of listed types only. Restriction
%    applies on top of others. Optional
%  RES fields:
%  - RSTAT: Return status
%    - 0: Converged to accuracy DELTAEPS
%    - 1: Done MAXIT sweeps
%  - NIT:   Number of sweeps done
%  - DELTA: Return from final EPT.INF_COUP_SEQUENTIAL call (max
%    over potentials)
%  - NSKIP: Histogram of SKIP entries returned by
%    EPT.INF_COUP_SEQUENTIAL calls, summed over SKIP vectors and
%    sweeps. 4 entries for values 0:3. NSKIP(1) is total number of
%    non-skipped updates, ...
%  RES_DET fields:
%  - DELTA: Return EPT.INF_COUP_SEQUENTIAL for every sweep (max
%    over potentials)
%  - NSKIP: Matrix, one row per sweep. Row is histogram, see
%    RES.NSKIP. RES.NSKIP is sum of rows here
%
%  Factorized posterior (IMODE == 'Factorized'):
%  Calls internal EPT.INF_FACT_SEQUENTIAL.
%  We run sweeps over all potentials. For each sweep, we pick a
%  random partition of these potentials.
%  OPTS fields:
%  - MAXIT: Maximum number of sweeps (over all potentials)
%  - DELTAEPS: Converged if max of DELTA < DELTAEPS, return of
%    EPT.INF_FACT_SEQUENTIAL
%  - DAMP: Damping constant. Def.: 0
%  - PIMINTHRES: See EPT.INF_FACT_SEQUENTIAL. Def.: 1e-8
%  - REFRESH: If true, the marginals are recomputed from the EP
%    message parameters after each sweep. Def.: true
%  - SKIP_GAUSS: If true, we do not update on potentials of type
%    'Gaussian'. Def.: false
%  - UPD_1STSWEEP: See "coupled, sequential". Optional
%  RES fields:
%  - RSTAT: Return status
%    - 0: Converged to accuracy DELTAEPS
%    - 1: Done MAXIT sweeps
%  - NIT:   Number of sweeps done
%  - DELTA: Return from final EPT.INF_FACT_SEQUENTIAL call (max
%    over potentials)
%  - NSKIP: Histogram of SKIP entries returned by
%    EPT.INF_FACT_SEQUENTIAL calls, summed over SKIP vectors and
%    sweeps. 5 entries for values 0:4. NSKIP(1) is total number of
%    non-skipped updates, ...
%  - NSDAMP: Only if selective damping is active (fields SD_XXX in
%    REPR). Total number of non-skipped updates for which selective
%    damping was applied.
%  RES_DET fields:
%  - DELTA: Return EPT.INF_FACT_SEQUENTIAL for every sweep (max
%    over potentials)
%  - NSKIP: Matrix, one row per sweep. Row is histogram, see
%    RES.NSKIP. RES.NSKIP is sum of rows here
%  - NSDAMP: Only if selective damping is active (fields SD_XXX in
%    REPR). For each sweep: Number of non-skipped updates for which
%    selective damping was applied.
%  TODO:
%  - Damping constants depend on potential type
%  - Allow for a skip list (right now, just OPTS.SKIP_GAUSS)
%  - Allow for scheduling heuristics other than random

%global deb_ept_pi deb_ept_beta deb_ept_h deb_ept_rho; % DEBUG!

%deb_ept_pi = []; deb_ept_beta = []; deb_ept_h = []; % DEBUG!
%deb_ept_rho = [];

global debug_pydeb debug_pydeb_iter;

[m,n] = size(model.matB);
if nargin~=4
  error('Need 4 input arguments');
end
if nargout<2
  error('Need at least 2 return arguments');
end
% Check common OPTS fields
if ~isstruct(opts)
  error('OPTS must be structure');
end
if ~isfield(opts,'verbose')
  opts.verbose = 0;
end
if ~isfield(opts,'maxit') || ~isnumeric(opts.maxit) || ...
      opts.maxit~=ceil(opts.maxit) || opts.maxit<1
  error('OPTS: maxit missing or wrong');
end
if ~isfield(opts,'damp')
  opts.damp = 0;
elseif ~isnumeric(opts.damp) || opts.damp<0 || opts.damp>=1
  error('OPTS: damp wrong');
end
if isfield(opts,'test_model')
  [m2,n2] = size(opts.test_model.matB);
  if n2~=n
    error('OPTS.TEST_MODEL: Wrong number of variables');
  end
  test_targs=[];
  for k=1:numel(opts.test_model.potMan)
    vv = opts.test_model.potMan{k}.pars{1};
    test_targs=[test_targs; vv(:)];
  end
end
if ~strcmp(imode,'Factorized')
  if ~isfield(opts,'deltaeps') || ~isnumeric(opts.deltaeps) || ...
	opts.deltaeps<=0
    error('OPTS: deltaeps missing or wrong');
  end
  if ~isfield(opts,'caveps')
    opts.caveps = 1e-5;
  elseif ~isnumeric(opts.caveps) || opts.caveps<=0
    error('OPTS: caveps wrong');
  end
end

% Gateway: Case distinction
switch(imode)
 case 'CoupParallel'
  % Coupled mode, parallel updating
  res.rstat = 1;
  res.nskip = 0;
  if nargout>3
    res_det.delta = [];
    res_det.nskip = [];
  end
  for nit = 1:opts.maxit
    if debug_pydeb
      debug_pydeb_iter = nit;
    end
    [model,repr,res.delta,nskip] = ...
	ept.inf_coup_parallel(model,repr,opts.damp,opts.caveps);
    res.nskip = res.nskip+nskip;
    if nargout>3
      res_det.delta = [res_det.delta; res.delta];
      res_det.nskip = [res_det.nskip; nskip];
    end
    if opts.verbose>0
      fprintf(1,'It. %d: delta=%f, nskip=%d\n',nit,res.delta, ...
	      nskip);
    end
    if isfield(opts,'test_model')
      [h_q,rho_q,logz,h_p,rho_p] = ...
	  ept.predict_ep(model,repr,imode,3,opts.test_model);
      nte = numel(test_targs);
      acc = 100*sum(test_targs==sign(h_q))/nte;
      loglh = sum(logz)/nte;
      fprintf(1,['Test set: Accuracy: %4.2f%%,  ' ...
		 'Log likelihood: %.6f\n'],acc,loglh);
    end
    % DEBUG:
    %indok = model.potInd;
    %deb_ept_pi = [deb_ept_pi repr.epPi(indok)];
    %deb_ept_beta = [deb_ept_beta repr.epBeta(indok)];
    %deb_ept_h = [deb_ept_h repr.mMeans(indok)];
    %deb_ept_rho = [deb_ept_rho repr.mVars(indok)];
    % END DEBUG
    if res.delta<opts.deltaeps
      res.rstat = 0;
      break; % Converged: Leave sweeps loop
    end
  end
  res.nit = nit;
 case {'CoupSequential','CoupSeqMargsUp2Date'}
  % Coupled mode, sequential updating
  if ~isfield(opts,'skipeps')
    opts.skipeps = 1e-8;
  elseif ~isnumeric(opts.skipeps) || opts.skipeps<=0
    error('OPTS: skipeps wrong');
  end
  if ~isfield(opts,'refresh')
    opts.refresh = 1;
  end
  if isfield(opts,'upd_1stsweep')
    if ~iscellstr(opts.upd_1stsweep) || isempty(opts.upd_1stsweep)
      error('OPTS.UPD_1STSWEEP wrong');
    end
    ind_swp1 = ept.potman_filterpots(model.potMan, ...
				     opts.upd_1stsweep);
    if isempty(ind_swp1)
      error('OPTS.UPD_1STSWEEP: No updates in first sweep?');
    end
    do_1stsweep = 1;
  else
    do_1stsweep = 0;
  end
  % MODEL.potInd: Non-Gaussian potentials
  model = ept.check_model_repres(model,repr,imode);
  res.rstat = 1;
  res.nskip = zeros(1,4);
  if nargout>3
    res_det.delta = [];
    res_det.nskip = [];
  end
  for nit = 1:opts.maxit
    rp = randperm(numel(model.potInd));
    updind = model.potInd(rp);
    if do_1stsweep && nit==1
      updind = intersect(updind,ind_swp1);
    end
    [model,repr,delta,skip] = ...
	ept.inf_coup_sequential(model,repr,updind,imode, ...
				opts.damp,opts.caveps, ...
				opts.skipeps);
    res.delta = max(delta);
    nskip = hist(skip,0:3); nskip = nskip(:)';
    res.nskip = res.nskip + nskip;
    if nargout>3
      res_det.delta = [res_det.delta; res.delta];
      res_det.nskip = [res_det.nskip; nskip];
    end
    if opts.refresh
      repr = ept.refresh_repres(model,repr,imode);
    end
    if opts.verbose>0
      fprintf(1,'It. %d: delta=%f, nnskip=%d\n',nit,res.delta, ...
	      sum(nskip(2:end)));
      fprintf(1,'   nskip=[');
      fprintf('%d ',nskip); fprintf(1,']\n');
    end
    if isfield(opts,'test_model')
      [h_q,rho_q,logz,h_p,rho_p] = ...
	  ept.predict_ep(model,repr,imode,3,opts.test_model);
      nte = numel(test_targs);
      acc = 100*sum(test_targs==sign(h_q))/nte;
      loglh = sum(logz)/nte;
      fprintf(1,['Test set: Accuracy: %4.2f%%,  ' ...
		 'Log likelihood: %.6f\n'],acc,loglh);
    end
    if res.delta<opts.deltaeps
      res.rstat = 0;
      break; % Converged: Leave sweeps loop
    end
  end
  res.nit = nit;
 case 'Factorized'
  % Factorized mode, sequential updating
  if ~isfield(opts,'piminthres')
    opts.piminthres = 1e-8;
  elseif ~isnumeric(opts.piminthres) || opts.piminthres<=0
    error('OPTS: piminthres wrong');
  end
  if ~isfield(opts,'refresh')
    opts.refresh = 1;
  end
  if ~isfield(opts,'skip_gauss')
    opts.skip_gauss = 0;
  end
  if isfield(opts,'upd_1stsweep')
    if ~iscellstr(opts.upd_1stsweep) || isempty(opts.upd_1stsweep)
      error('OPTS.UPD_1STSWEEP wrong');
    end
    ind_swp1 = int32(ept.potman_filterpots(model.potMan, ...
					   opts.upd_1stsweep));
    if isempty(ind_swp1)
      error('OPTS.UPD_1STSWEEP: No updates in first sweep?');
    end
    do_1stsweep = 1;
    % DEBUG:
    %fprintf(1,'DEBUG: Index in 1st sweep:\n');
    %ind_swp1'
    % END DEBUG
  else
    do_1stsweep = 0;
  end
  do_seldamp = isfield(repr,'sd_numk'); % Selective damping?
  do_deb_testmaxpi = do_seldamp && isfield(opts,'debug_testmaxpi') ...
      && opts.debug_testmaxpi;
  do_testmodel = isfield(opts,'test_model');
  do_deb_plotdiff = isfield(opts,'debug_plotabsdiff') && ...
      opts.debug_plotabsdiff;
  % Check MODEL, compute internal representations
  % MODEL.potInd: Non-Gaussian potentials
  model = ept.check_model_repres(model,repr,imode);
  res.rstat = 1;
  res.nskip = zeros(1,5); % BAD: Could be more states in future
  if do_seldamp
    res.nsdamp = 0;
  end
  if nargout>3
    res_det.delta = [];
    res_det.nskip = [];
    if do_seldamp
      res_det.nsdamp = [];
    end
  end
  for nit = 1:opts.maxit
    % DEBUG:
    if isfield(opts,'debug_plotabsdiff') && opts.debug_plotabsdiff
      old_mean = repr.mBeta./repr.mPi;
      old_std = 1./sqrt(repr.mPi);
    end
    %if debug_pydeb
    %  debug_pydeb_iter = nit;
    %end
    if ~opts.skip_gauss
      updind = int32(randperm(m)); % Permutation of all potentials
    else
      rp = randperm(numel(model.potInd));
      updind = int32(model.potInd(rp));
    end
    if do_1stsweep && nit==1
      updind = int32(intersect(updind,ind_swp1));
    end
    if ~do_seldamp
      [model,repr,delta,skip] = ...
	  ept.inf_fact_sequential(model,repr,updind,opts.damp, ...
				  opts.piminthres);
    else
      [model,repr,delta,skip,sd_damp] = ...
	  ept.inf_fact_sequential(model,repr,updind,opts.damp, ...
				  opts.piminthres);
      % Among non-skipped updates, count those with SD_DAMP larger
      % than OPTS.DAMP (this is where selective damping was used)
      nsdamp = sum(sd_damp(find(skip==0))>opts.damp);
      res.nsdamp = res.nsdamp+nsdamp;
      if nargout>3
	res_det.nsdamp = [res_det.nsdamp; nsdamp];
      end
    end
    res.delta = max(delta);
    nskip = hist(skip,0:4); nskip = nskip(:)';
    res.nskip = res.nskip + nskip;
    if nargout>3
      res_det.delta = [res_det.delta; res.delta];
      res_det.nskip = [res_det.nskip; nskip];
    end
    % DEBUG:
    if do_deb_testmaxpi
      if test_fact_maxpi(model,repr,opts.debug_testmaxpi)
	fprintf(1,['DEBUG(inf_ep), it=%d: Max-pi data structure is' ...
		   ' wrong!'],nit);
      end
    end
    if opts.refresh
      repr = ept.refresh_repres(model,repr,imode);
    end
    if opts.verbose>0
      fprintf(1,'It. %d: delta=%f, nnskip=%d\n',nit,res.delta, ...
	      sum(nskip(2:end)));
      fprintf(1,'   nskip=[');
      fprintf(1,'%d ',nskip); fprintf(1,']');
      if do_seldamp
	fprintf(1,', nsdamp=%d\n',nsdamp);
      else
	fprintf(1,'\n');
      end
    end
    if do_testmodel
      [h_q,rho_q,logz,h_p,rho_p] = ...
	  ept.predict_ep(model,repr,imode,3,opts.test_model);
      nte = numel(test_targs);
      acc = 100*sum(test_targs==sign(h_q))/nte;
      loglh = sum(logz)/nte;
      fprintf(1,['Test set: Accuracy: %4.2f%%,  ' ...
		 'Log likelihood: %.6f\n'],acc,loglh);
    end
    % DEBUG:
    if debug_pydeb
      updind = updind-1;
      ep_pi = repr.epPi; ep_beta = repr.epBeta;
      marg_pi = repr.mPi; marg_beta = repr.mBeta;
      save(sprintf('eptb13_vars_it%d.mat',nit),'ep_pi', ...
	   'ep_beta','marg_pi','marg_beta','updind','-v7');
    end
    % DEBUG:
    if do_deb_plotdiff
      df_mean = repr.mBeta./repr.mPi - old_mean;
      df_std = 1./sqrt(repr.mPi) - old_std;
      [tmp,ind] = sort(abs(df_mean),1,'descend');
      figure(1);
      plot(1:25,df_mean(ind(1:25)),'o');
      title('Means: Absolute differences');
      lab = [];
      for i=1:25
	lab{i} = num2str(ind(i));
      end
      set(gca,'XTick',1:25);
      set(gca,'XTickLabel',lab);
      [tmp,ind] = sort(abs(df_std),1,'descend');
      figure(2);
      plot(1:25,df_std(ind(1:25)),'o');
      title('Standard devs: Absolute differences');
      lab = [];
      for i=1:25
	lab{i} = num2str(ind(i));
      end
      set(gca,'XTick',1:25);
      set(gca,'XTickLabel',lab);
      pause;
    end
    if res.delta<opts.deltaeps
      res.rstat = 0;
      break; % Converged: Leave sweeps loop
    end
  end
  res.nit = nit;
 otherwise
  error('Unknown IMODE value');
end
