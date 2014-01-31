%MatDiag: Diagonal matrix

classdef MatDiag < ept.Mat
  properties (SetAccess = private, GetAccess = protected)
    dg;  % Diagonal
  end

  methods
    function B = MatDiag(dg)
      m = numel(dg);
      if m==0
	error('DG must not be empty');
      end
      B = B@ept.Mat(m,m); % Superclass constructor
      B.dg = dg(:);
    end

    function bv = mvm(B,v)
      bv = B.dg.*v;
    end

    function v = getCol(B,i)
      if i~=ceil(i) || i<1 || i>B.m
	error('I wrong');
      end
      v = zeros(B.m,1); v(i)=B.dg(i);
    end

    function S = matBtDgB(B,v)
      S = diag(v(:).*(B.dg.^2));
    end

    function v = diagBSymBt(B,S)
      if size(S,1)~=B.m || size(S,2)~=B.m
	error('S has wrong size');
      end
      v = diag(S).*(B.dg.^2);
    end
  end
end % classdef
