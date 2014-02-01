%MatEye: Identity matrix

classdef MatEye < ept.Mat
  methods
    function B = MatEye(m)
      B = B@ept.Mat(m,m); % Superclass constructor
    end

    function bv = mvm(B,v)
      if numel(v)~=B.m
	error('V has wrong size');
      end
      bv = v;
    end

    function v = getCol(B,i)
      if i~=ceil(i) || i<1 || i>B.m
	error('I wrong');
      end
      v = zeros(B.m,1); v(i)=1;
    end

    function S = matBtDgB(B,v)
      if numel(v)~=B.m
	error('V has wrong size');
      end
      S = diag(v);
    end

    function v = diagBSymBt(B,S)
      if size(S,1)~=B.m || size(S,2)~=B.m
	error('S has wrong size');
      end
      v = diag(S);
    end
  end % methods
end % classdef
