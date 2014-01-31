%MatDef: Wrapper for normal matrix (dense or sparse)

classdef MatDef < ept.Mat
  properties (SetAccess = private, GetAccess = protected)
    mx;  % Explicit matrix
  end

  methods
    function B = MatDef(mx)
      if ndims(mx)>2
	error('MX must be matrix');
      end
      [m,n] = size(mx);
      B = B@ept.Mat(m,n); % Superclass constructor
      B.mx = mx;
    end

    function bv = mvm(B,v)
      if ~B.transp
	bv = full(B.mx*v);
      else
	bv = full(B.mx'*v);
      end
    end

    function v = getCol(B,i)
      if ~B.transp
	v = full(B.mx(:,i));
      else
	v = full(B.mx(i,:)');
      end
    end

    function S = matBtDgB(B,v)
      if B.transp
	% This sucks: Not required anyway
	error('NOT IMPLEMENTED!');
      end
      if ~issparse(B.mx)
	S = B.mx'*bsxfun(@times,v(:),B.mx);
      else
	% S = full(B.mx'*(diag(v)*B.mx)); % Real slow!
	S = full(B.mx'*(repmat(v(:),1,size(B.mx,2)).*B.mx));
      end
    end

    function v = diagBSymBt(B,S)
      if B.transp
	% This sucks: Not required anyway
	error('NOT IMPLEMENTED!');
      end
      % Compute B*S intermediate. This may be bad if B is sparse
      % and m>>n, since B*S is dense
      v = full(sum((B.mx*S).*(B.mx),2));
    end
  end % methods
end % classdef
