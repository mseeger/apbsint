%MatSub: Selection operator

classdef MatSub < ept.Mat
  properties (SetAccess = private, GetAccess = protected)
    sind;  % Selection index
  end

  methods
    function B = MatSub(n,sind)
      if n~=ceil(n) || ~isempty(find(sind~=ceil(sind))) || ...
	    min(sind)<1 || max(sind)>n
	error('Wrong arguments');
      end
      B = B@ept.Mat(numel(sind),n); % Superclass constructor
      B.sind = sind;
    end

    function bv = mvm(B,v)
      if ~B.transp
	bv = v(B.sind);
      else
	bv = zeros(B.n,1);
	bv(B.sind) = v;
      end
    end

    function v = getCol(B,i)
      if i~=ceil(i) || i<1 || i>size(B,2)
	error('I wrong');
      end
      if ~B.transp
	v = double(B.sind==i);
      else
	v = zeros(B.n,1);
	v(B.sind(i)) = 1;
      end
    end

    function S = matBtDgB(B,v)
      if B.transp
	% This sucks: Not required anyway
	error('NOT IMPLEMENTED!');
      end
      tv = zeros(B.n,1);
      tv(B.sind) = v;
      S = zeros(B.n);
      S(1:B.n+1:end) = tv;
    end

    function v = diagBSymBt(B,S)
      if B.transp
	% This sucks: Not required anyway
	error('NOT IMPLEMENTED!');
      end
      if size(S,1)~=B.n || size(S,2)~=B.n
	error('S has wrong size');
      end
      tv = diag(S);
      v = tv(B.sind);
    end
  end % methods
end % classdef
