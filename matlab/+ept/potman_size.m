function sz = potman_size(pman)
%POTMAN_SIZE Number of potentials of manager
%  SZ = EPT.POTMAN_SIZE(PMAN)
%  Returns size (number of potentials) of manager.

if ~iscell(pman)
  sz = pman.size;
else
  sz = sum(cellfun(@(x) x.size,pman));
end
