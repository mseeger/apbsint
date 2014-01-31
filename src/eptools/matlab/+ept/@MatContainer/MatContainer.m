%MatContainer: Container factor [B1; B2; ...]

classdef MatContainer < ept.Mat
  properties (SetAccess = private, GetAccess = protected)
    child;  % Cell of child ept.Mat objects
    chPos;
  end

  methods
    function B = MatContainer(child)
      if ~iscell(child)
	error('CHILD must be cell array');
      end
      fst = 1; m = 0;
      for k = 1:numel(child)
	obj = child{k};
	if ~isa(obj,'ept.Mat')
	  error('CHILD must be array of Mat objects');
	end
	m = m+size(obj,1);
	if fst
	  n = size(obj,2); fst=0;
	elseif size(obj,2)~=n
	  error(['All CHILD entries must have same number of' ...
		 ' columns']);
	end
      end
      B = B@ept.Mat(m,n); % Superclass constructor
      B.child = child;
      % chPos: See getRelPos
      B.chPos = zeros(m,2);
      off = 0;
      for k = 1:numel(child)
	sz = size(child{k},1);
	B.chPos(off+1:off+sz,:) = [repmat(k,sz,1) (1:sz)'];
	off = off+sz;
      end
    end
  end % methods

  methods (Access = private)
    % Map row index i -> [k,j]: Row j in child{k}
    function [k,j] = getRelPos(B,i)
      k = B.chPos(i,1); j = B.chPos(i,2);
    end
  end
end % classdef
