% EPTOOLS Matlab Interface
% Example: Binary classification, coupled or factorized posterior
% mode.
% Can be configured with Gaussian or Laplace prior.
%
% Same data as glm-ie_v1.5/doc/classify_inference.m:
% UCI Adult (a9a): Dimension 123, first 30000 for testing, last 2561
% for training.

clear all, close all

%global debug_bmat; % DEBUG!
%global deb_ept_pi deb_ept_beta deb_ept_h deb_ept_rho; % DEBUG!
global debug_pydeb debug_pydeb_iter;
debug_pydeb = 0;
%debug_pydeb = 1;
debug_pydeb_iter = 0;

% Load and split into train and test set
load('/home/seeger/lhotse/src/eptools/matlab/test/binclass/classify.mat');
nte = 30000;      % say how many test examples we wich to keep
Bte = B(1:nte,:); B = B(nte+1:end,:); cte = c(1:nte);
c = c(nte+1:end); % split
% Remove variables not touched by any observations
ind = find(sum(B~=0,1));
B = B(:,ind); Bte = Bte(:,ind);
[q,n] = size(B);
% DEBUG:
%debug_bmat = [eye(n); B];

% Nearest neighbour classification
dp = sum(bsxfun(@minus, Bte, full(mean(B(c==+1,:)))).^2,2);  % dist. to pos. ctr
dm = sum(bsxfun(@minus, Bte, full(mean(B(c==-1,:)))).^2,2);  % dist. to neg. ctr
cc = 2*double(dp<dm)-1;
fprintf('nearest neighbor accuracy=%1.2f%%\n',100*sum(cte==cc)/ ...
	numel(cte))

% Setup
%imode = 'CoupParallel';   % Coupled, parallel updating
imode = 'CoupSequential'; % Coupled, sequential updating
%imode = 'Factorized';     % Factorized backbone
doLaplace = 1;            % Laplace prior
tau_lapl = 2/5;
%doLaplace = 0;            % Gaussian prior
%ssq_gauss = 25/4;
isFact = strcmp(imode,'Factorized');
if isFact
  seldampK = 5;           % Selective damping: K
end

% Coupling factors
if ~isFact
  model.matB = ept.MatContainer({ept.MatEye(n),ept.MatDef(B)});
  pmodel.matB = ept.MatDef(Bte);
else
  % Sparse matrices
  model.matB = [eye(n); B];
  pmodel.matB = Bte;
end
m = size(model.matB,1);
% Potential managers
if ~doLaplace
  model.potMan{1}.name = 'Gaussian';
  model.potMan{1}.size = n;
  model.potMan{1}.pars{1}=0;
  model.potMan{1}.pars{2}=ssq_gauss;
else
  model.potMan{1}.name = 'Laplace';
  model.potMan{1}.size = n;
  model.potMan{1}.pars{1}=0;
  model.potMan{1}.pars{2}=tau_lapl;
end
model.potMan{2}.name = 'Probit';
model.potMan{2}.size = q;
model.potMan{2}.pars{1}=c;
model.potMan{2}.pars{2}=0;
pmodel.potMan{1}.name = 'Probit';
pmodel.potMan{1}.size = nte;
pmodel.potMan{1}.pars{1}=cte;
pmodel.potMan{1}.pars{2}=0;
fprintf(1,'Prior:      %s\nLikelihood: %s\n', ...
	model.potMan{1}.name,model.potMan{2}.name);
% Initialization:
if ~doLaplace
  % NOTE: In factorized mode, OPTS.CAV_VAR is not needed, since all
  % Gaussian potentials have support 1
  [model,repr] = ept.init_ep(model,imode,'ADF');
else
  % Initialization as in glm-ie
  if ~isFact
    repr.epPi = ones(m,1);
    repr.epBeta = zeros(m,1);
    repr.epBeta(n+1:end) = 0.5*c;
  else
    % Factorized: Do not update on prior potentials in 1st sweep
    repr.epPi = zeros(nnz(model.matB),1);
    repr.epBeta = zeros(nnz(model.matB),1);
    repr.epPi(1:n) = 1;
  end
  repr = ept.refresh_repres(model,repr,imode);
end
% DEBUG:
%[repr.epPi(1:126) repr.epBeta(1:126)]
%pause;
if debug_pydeb
  if strcmp(imode,'CoupParallel')
    ep_pi = repr.epPi; ep_beta = repr.epBeta;
    lfact = repr.lFact; cvec = repr.cVec;
    mmeans = repr.mMeans; mvars = repr.mVars;
    save('eptb9_vars_it0.mat','ep_pi','ep_beta','lfact','cvec', ...
	 'mmeans','mvars','-v7');
  elseif strcmp(imode,'Factorized')
    ep_pi = repr.epPi; ep_beta = repr.epBeta;
    marg_pi = repr.mPi; marg_beta = repr.mBeta;
    save('eptb13_vars_it0.mat','ep_pi','ep_beta','marg_pi', ...
	 'marg_beta','-v7');
  end
end
% Predictions (before we start)
[h_q,rho_q,logz,h_p,rho_p] = ept.predict_ep(model,repr,imode,3, ...
					    pmodel);
acc = 100*sum(cte==sign(h_q))/numel(cte);
loglh = sum(logz)/numel(cte);
fprintf(1,['Initial predictions:\nAccuracy: %4.2f%%\n' ...
	   'Test set log likelihood: %.6f\n'],acc,loglh);

% EP inference
fprintf(1,'\nRunning EP inference (mode: %s)...\n',imode);
switch imode
 case 'CoupParallel'
  opts.maxit = 20;       % Max. number sweeps
  opts.deltaeps = 1e-4;  % Convergence threshold
  opts.damp = 0;         % No damping
  opts.caveps = 1e-5;
  opts.verbose = 1;
  %opts.test_model = pmodel;  % Test set predictions after each sweep
  t_start = tic;
  [model,repr,res,res_det] = ept.inf_ep(model,repr,imode,opts);
  t_inf = toc(t_start);
  fprintf(1,'\nTime for inference: %.6f\n',t_inf);
  if res.rstat==1
    fprintf(1,'\nDone MAXIT iterations.\n');
  else
    fprintf(1,'\nConverged to DELTAEPS accuracy.\n');
  end
  fprintf(1,['Number sweeps:   %d\n' ...
	     'Final delta:     %f\n' ...
	     'Number of skips: %d\n'],res.nit,res.delta,res.nskip);
 case 'CoupSequential'
  opts.maxit = 15;       % Max. number sweeps
  opts.deltaeps = 1e-4;  % Convergence threshold
  opts.damp = 0;         % No damping
  opts.caveps = 1e-5;
  opts.skipeps = 1e-7;
  opts.refresh = 1;      % Recompute repres. after each sweep
  opts.verbose = 1;
  %opts.test_model = pmodel;  % Test set predictions after each sweep
  t_start = tic;
  [model,repr,res,res_det] = ept.inf_ep(model,repr,imode,opts);
  t_inf = toc(t_start);
  fprintf(1,'\nTime for inference: %.6f\n',t_inf);
  if res.rstat==1
    fprintf(1,'\nDone MAXIT iterations.\n');
  else
    fprintf(1,'\nConverged to DELTAEPS accuracy.\n');
  end
  fprintf(1,['Number sweeps:   %d\n' ...
	     'Final delta:     %f\n' ...
	     'Skip histogram:  ['],res.nit,res.delta);
  fprintf(1,'%d ',res.nskip); fprintf(1,']\n');
 case 'Factorized'
  opts.maxit = 100;      % Max. number sweeps
  opts.deltaeps = 1e-4;  % Convergence threshold
  opts.damp = 0;         % No damping
  opts.piminthres = 1e-8;
  opts.refresh = 1;      % Recompute repres. after each sweep
  opts.verbose = 1;
  if ~doLaplace
    opts.skip_gauss = 1; % Skip Gaussian prior potentials
  else
    opts.upd_1stsweep = {'Probit'}; % Skip Laplace pots. 1st sweep
  end
  if exist('seldampK')
    % Selective damping mechanism
    repr.sd_numk = seldampK;
    model.matBInt = ept.bfact_intrepres(model.matB);
    irb = model.matBInt;
    if ~doLaplace
      % Exclude Gaussian prior potentials
      repr.sd_subind = ...
	  int32(ept.potman_filterpots(model.potMan, ...
				      {'Gaussian'})-1);
      repr.sd_subexcl = 1;
      [repr.sd_numvalid,repr.sd_topind,repr.sd_topval] = ...
	  eptools_fact_compmaxpi(n,m,irb.rowind,irb.colind, ...
				 irb.bvals,repr.epPi,repr.epBeta, ...
				 seldampK,repr.sd_subind, ...
				 repr.sd_subexcl);
    else
      [repr.sd_numvalid,repr.sd_topind,repr.sd_topval] = ...
	  eptools_fact_compmaxpi(n,m,irb.rowind,irb.colind, ...
				 irb.bvals,repr.epPi,repr.epBeta, ...
				 seldampK);
    end
    % DEBUG:
    %opts.debug_testmaxpi = 2; % Test after each sweep, full output
  end
  % DEBUG:
  %opts.debug_plotabsdiff = 1; % Plots after each sweep
  %opts.test_model = pmodel;  % Test set predictions after each sweep
  t_start = tic;
  [model,repr,res,res_det] = ept.inf_ep(model,repr,imode,opts);
  t_inf = toc(t_start);
  fprintf(1,'\nTime for inference: %.6f\n',t_inf);
  if res.rstat==1
    fprintf(1,'\nDone MAXIT iterations.\n');
  else
    fprintf(1,'\nConverged to DELTAEPS accuracy.\n');
  end
  fprintf(1,['Number sweeps:   %d\n' ...
	     'Final delta:     %f\n' ...
	     'Skip histogram:  ['],res.nit,res.delta);
  fprintf(1,'%d ',res.nskip); fprintf(1,']\n');
  if exist('seldampK')
    fprintf(1,'Select. damps:   %d\n',res.nsdamp);
  end
 otherwise
  error('Unknown IMODE');
end

% Predictions
[h_q,rho_q,logz,h_p,rho_p] = ept.predict_ep(model,repr,imode,3, ...
					    pmodel);
acc = 100*sum(cte==sign(h_q))/numel(cte);
loglh = sum(logz)/numel(cte);
fprintf(1,['\nFinal predictions:\nAccuracy: %4.2f%%\n' ...
	   'Test set log likelihood: %.6f\n'],acc,loglh);
