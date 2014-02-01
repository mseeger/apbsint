%Mat: Base class for B coupling factors

classdef Mat
  properties (SetAccess = private, GetAccess = public)
    m = 1;       % Number of rows
    n = 1;       % Number of columns
    transp = 0;  % Transpose (n-by-m)?
  end % properties

  methods
    function B = Mat(m,n)
      if m~=ceil(m) || m<=0 || n~=ceil(n) || n<=0
	error('Wrong size M, N');
      end
      B.m = m; B.n = n;
      B.transp = 0;
    end

    function [m,n] = size(B,id)
      sz = [B.m B.n];
      if B.transp, sz = fliplr(sz); end
      if nargin==1
	if nargout==2
	  m = sz(1); n = sz(2);
	else
	  m = sz;
	end
      elseif nargin==2
	if nargout>1, error('Too many output arguments'); end
	if id(1)>numel(sz)
	  m = 1;
	else
	  m = sz(id);
	end
      else
	error('Wrong number of input arguments')
      end
    end

    function B = ctranspose(B)
      B.transp = xor(B.transp,1);
    end
    
    % Further public methods:
    % getCol
    % matBtDgB
    % diagBMxBt
  end % methods

  % Abstract methods
  methods (Abstract = true)
    % MVM: Compute bv = B*v     [transp==false]
    %              bv = B'*v    [transp==true]
    bv = mvm(B,v);
  end
end % classdef
