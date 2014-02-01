function model = check_model_repres(model,repr,mode,recomp)
%CHECK_MODEL_REPRES Check model and representation
%  [{MODEL}] = EPT.CHECK_MODEL_REPRES(MODEL,REPR,MODE,{RECOMP=0})
%  Checks whether model MODEL and representation REPR are
%  consistent. These structs depend on the mode MODE.
%
%  If there is a return argument, internal and auxiliary
%  representations are computed for MODEL, unless these are already
%  present and RECOMP is false.

if nargin<4
  recomp = 0;
end
switch mode
 case 'CoupParallel'
  mod = 1;
 case 'CoupSequential'
  mod = 2;
 case 'CoupSeqMargsUp2Date'
  mod = 3;
 case 'Factorized'
  mod = 4;
 otherwise
  error('Unknown MODE');
end
isCoup = (mod<4);
if ~isstruct(model) || ~isfield(model,'potMan') || ...
      ~isfield(model,'matB') || ...
      (isCoup && ~isa(model.matB,'ept.Mat')) || ...
      (~isCoup && ~issparse(model.matB))
  error('MODEL has wrong fields');
end
[m,n] = size(model.matB);
if ept.potman_size(model.potMan)~=m
  error('MODEL: potMan, matB have incomplete sizes');
end
if isCoup
  sz = m;
else
  sz = nnz(model.matB);
end
if numel(repr.epPi)~=sz
  error('REPR: epPi');
end
if numel(repr.epBeta)~=sz
  error('REPR: epBeta');
end
if isCoup
  if size(repr.lFact,1)~=n || size(repr.lFact,2)~=n
    error('REPR: lFact');
  end
  if numel(repr.cVec)~=n
    error('REPR: cVec');
  end
  if mod~=2
    if numel(repr.mMeans)~=m
      error('REPR: mMeans');
    end
    if numel(repr.mVars)~=m
      error('REPR: mVars');
    end
  end
else
  if numel(repr.mBeta)~=n
    error('REPR: mBeta');
  end
  if numel(repr.mPi)~=n
    error('REPR: mPi');
  end
end

if nargout>0
  if recomp || ~isfield(model,'potManInt') || ...
	~isfield(model,'potInd')
    [model.potManInt,model.potInd] = ...
	ept.potman_intrepres(model.potMan);
  end
  if ~isCoup && (recomp || ~isfield(model,'matBInt'))
    model.matBInt = ept.bfact_intrepres(model.matB);
  end
end
