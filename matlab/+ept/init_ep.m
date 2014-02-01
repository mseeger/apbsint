function [model,repr] = init_ep(model,imode,init_proc,opts)
%INIT_EP Initialisation of representation for EP inference
%  [MODEL,REPR] = EPT.INIT_EP(MODEL,IMODE,INIT_PROC,{OPTS})
%  Initialization of EP parameters in representation REPR, using a
%  initialization procedure selected by INIT_PROC (may depend on
%  mode IMODE as well). Some initialization procedures require
%  additional options in OPTS.
%
%  Coupled posterior (parallel or sequential updating):
%  IMODE == 'CoupXXX'. INIT_PROC selects:
%  - 'ADF': EP parameters for non-Gaussian potentials are set to 0,
%    those for Gaussian potentials fixed to their constant
%    value. Works only if the Gaussian part alone is normalizable
%    and sufficiently well-conditioned
%
%  Factorized posterior:
%  IMODE == 'Factorized'. INIT_PROC selects:
%  - 'ADF': EP parameters for non-Gaussian potentials are set to 0,
%    those for Gaussian potentials fixed to their constant
%    value. For potentials j with V_j = {i} and b_ji = 1, the EP
%    parameters remain fixed there. In general, this is just a
%    heuristic, which also depends on OPTS.CAV_VAR (see design
%    technical report for details).
%
%  TODO: Add further initialization procedures for common use cases!

[m,n] = size(model.matB);
if nargin<3
  error('Need 3 input arguments');
end
if nargout<2
  error('Need 2 return arguments');
end

% Gateway: Case distinction
switch(imode)
 case {'CoupParallel','CoupSequential','CoupSeqMargsUp2Date'}
  % Coupled mode
  switch(init_proc)
   case 'ADF'
    % Loop over Gaussian potentials. The EP parameters are pi =
    % 1/ssq, beta = y/ssq.
    pman = model.potMan;
    nb = numel(pman);
    if ~iscell(pman)
      ppm{1} = pman;
    else
      ppm = pman(:);
    end
    pidGauss = eptools_getpotid('Gaussian');
    repr.epPi = []; repr.epBeta = [];
    for k=1:nb
      pid = eptools_getpotid(ppm{k}.name);
      sz = ppm{k}.size;
      if pid == -1
	error(sprintf('Block %d: Unknown potential name',k));
      end
      if pid == pidGauss
	y = ppm{k}.pars{1};
	if numel(y)==1
	  y = repmat(y,sz,1);
	else
	  y = y(:);
	end
	ssq = ppm{k}.pars{2};
	if numel(ssq)==1
	  ssq = repmat(ssq,sz,1);
	else
	  ssq = ssq(:);
	end
	repr.epPi = [repr.epPi; 1./ssq];
	repr.epBeta = [repr.epBeta; y./ssq];
      else
	tv = zeros(sz,1);
	repr.epPi = [repr.epPi; tv];
	repr.epBeta = [repr.epBeta; tv];
      end
    end
  otherwise
    error('Unknown INIT_PROC value');
  end
 case 'Factorized'
  % Factorized mode
  switch(init_proc)
   case 'ADF'
    % Loop over Gaussian potentials. If potential j is
    % N(s|y_j,ssq_j), then the EP parameters are
    %   pi_ji = b_ji^2 / ( (|V_j|-1) cv + ssq_j)
    %   beta_ji = b_ji y_j / ( (|V_j|-1) cv + ssq_j)
    % Here, cv==OPTS.CAV_VAR. In the standard case |V_j|=1 and
    % b_ji=1 for V_j = {i}, this becomes 1/ssq_j and y_j/ssq_j
    % as in the coupled mode
    if nargin<4 || ~isfield(opts,'cav_var')
      opts.cav_var = 1;
    end
    pman = model.potMan;
    nb = numel(pman);
    if ~iscell(pman)
      ppm{1} = pman;
    else
      ppm = pman(:);
    end
    pidGauss = eptools_getpotid('Gaussian');
    repr.epPi = zeros(nnz(model.matB),1);
    repr.epBeta = zeros(nnz(model.matB),1);
    off = 0; % Offset into epPi, epBeta
    offj = 0; % Offset into potential set
    % bT is B^T. It is used to navigate and fill epPi, epBeta
    bT = model.matB';
    for k=1:nb
      pid = eptools_getpotid(ppm{k}.name);
      sz = ppm{k}.size;
      if pid == -1
	error(sprintf('Block %d: Unknown potential name',k));
      end
      if pid ~= pidGauss
	% Non-Gaussian: Just skip
	off = off + nnz(bT(:,offj+1:offj+sz));
      else
	%fprintf(1,'Gauss: sz=%d\n',sz); % DEBUG
	y = ppm{k}.pars{1};
	if numel(y)==1
	  y = repmat(y,sz,1);
	else
	  y = y(:);
	end
	ssq = ppm{k}.pars{2};
	%fprintf(1,'ssq=%f\n',ssq); % DEBUG
	if numel(ssq)==1
	  ssq = repmat(ssq,sz,1);
	else
	  ssq = ssq(:);
	end
	% Slice of B^T corr. to potential
	[ri,ci,vl] = find(bT(:,offj+1:offj+sz));
	sz2 = numel(ri);
	%fprintf(1,'sz2=%d\n',sz2); % DEBUG
	vjsz = hist(ci,1:sz)';
	% DEBUG
	%fprintf(1,'Should be empty for us now:\n');
	%find(vjsz~=1)
	%fprintf(1,'This as well:\n');
	%find(vl~=1)
	% END DEBUG
	tvec = 1./((vjsz-1)*opts.cav_var+ssq);
	repr.epPi(off+1:off+sz2) = (vl.^2).*tvec(ci);
	tvec = tvec.*y;
	repr.epBeta(off+1:off+sz2) = vl.*tvec(ci);
	off = off+sz2;
      end
      offj = offj+sz;
    end
   otherwise
    error('Unknown INIT_PROC value');
  end
 otherwise
  error('Unknown IMODE value');
end
% Compute representation from scratch
repr = ept.refresh_repres(model,repr,imode);
