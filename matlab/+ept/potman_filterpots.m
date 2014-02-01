function ind = potman_filterpots(pman,typlst,excl)
%POTMAN_FILTERPOTS Subindex for certain potential types
%  IND = EPT.POTMAN_FILTERPOTS(PMAN,TYPLST,{EXCL=0})
%  PMAN is a potential manager, TYPLST a cell array of potential
%  typenames. Returns position index of potentials whose type is in
%  TYPLST. If EXCL is true: ... whose type is not in TYPLST.

if nargin<2
  error('Must have 2 input arguments');
end
if ~iscellstr(typlst)
  error('TYPLST wrong');
end
if nargin<3
  excl = 0;
end
nb = numel(pman);
if ~iscell(pman)
  if nb~=1
    error('Wrong input PMAN');
  end
  ppm{1} = pman;
else
  ppm = pman(:);
end
% Loop over potential types
ind = [];
off = 0;
for k=1:nb
  numk = ppm{k}.size;
  if xor(excl,any(cellfun(@(x) strcmp(x,ppm{k}.name),typlst)))
    ind = [ind; (off+1:off+numk)'];
  end
  off = off+numk;
end
