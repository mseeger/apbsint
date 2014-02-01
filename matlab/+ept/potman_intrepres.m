function [irp,upind] = potman_intrepres(pman,deb_fid)
%POTMAN_INTREPRES Internal representation for potential manager
%  [IRP,{UPIND}] = EPT.POTMAN_INTREPRES(PMAN)
%  Creates internal representation IRP from potential manager
%  object PMAN (IRP can be passed to MEX functions).
%  Optionally, UPIND returns the index of potentials which are not
%  Gaussian (type 'Gaussian').
%
%  Use 'help ept.datatype_potentialmanager' for data type of
%  PMAN.

nb = numel(pman);
if ~iscell(pman)
  if nb~=1
    error('Wrong input PMAN');
  end
  ppm{1} = pman;
else
  ppm = pman(:);
end
% Loop 1: Everything except PARVEC, determine size
irp.potids = zeros(nb,1,'int32');
irp.numpot = zeros(nb,1,'int32');
parshrd = []; pvsz = 0;
if nargout>1
  upind = [];
  pidGauss = eptools_getpotid('Gaussian');
  %fprintf(1,'pidGauss=%d\n',pidGauss); % DEBUG!
  off = 0;
end
for k=1:nb
  pid = eptools_getpotid(ppm{k}.name);
  %fprintf(1,'k=%d, pid=%d\n',k,pid); % DEBUG!
  if pid == -1
    error(sprintf('Block %d: Unknown potential name',k));
  end
  irp.potids(k) = int32(pid);
  numk = ppm{k}.size;
  if numk~=floor(ppm{k}.size) || numk<=0
    error(sprintf('Block %d: size field invalid',k));
  end
  irp.numpot(k) = int32(numk);
  nump = numel(ppm{k}.pars);
  if nump>0
    if ~iscell(ppm{k}.pars)
      error(sprintf('Block %d: pars field must be cell array',k));
    end
    for p=1:nump
      sz = numel(ppm{k}.pars{p});
      if sz==1
	parshrd = [parshrd; 1];
      else
	if sz~=numk
	  error(sprintf('Block %d, pars{%d}: Wrong size',k,p));
	end
	parshrd = [parshrd; 0];
      end
      pvsz = pvsz+sz;
    end
  end
  if nargout>1
    if pid~=pidGauss
      upind = [upind; (off+1:off+numk)'];
    end
    off = off+numk;
  end
end
irp.parshrd = int32(parshrd);
% Loop 2: Assemble PARVEC
irp.parvec = zeros(pvsz,1);
off=0;
for k=1:nb
  nump = numel(ppm{k}.pars);
  for p=1:nump
    sz = numel(ppm{k}.pars{p});
    irp.parvec(off+1:off+sz) = ppm{k}.pars{p};
    off = off+sz;
  end
end
% Test whether all parameter values are valid
msg = eptools_potmanager_isvalid(irp.potids,irp.numpot,irp.parvec, ...
				 irp.parshrd);
if length(msg)~=0
  error(msg);
end
if nargin>1
  % DEBUG: Print input args
  ept.debug_printvec(irp.potids,deb_fid);
  ept.debug_printvec(irp.numpot,deb_fid);
  ept.debug_printvec(irp.parvec,deb_fid);
  ept.debug_printvec(irp.parshrd,deb_fid);
end
