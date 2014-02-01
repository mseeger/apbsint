% MS: Really INCOMPLETE implementation of probit potential, just
% for debugging
%
%   See also POTFUNCTIONS.M.

function P = potMSErf(s,type,z)

%fprintf(1,'PotMSErf: nargin=%d\n',nargin);
if nargin==1
  % Have no clue what to do here! Just copied from 'potLogistic'
  q = numel(s); P = zeros(q,4);                                % allocate memory
  P(:,1) = min(0,s(:)) - log(1+exp(-abs(s(:))));          % safe since abs(s)>=0
  p = exp(P(:,1));
  P(:,2) = 1-p;           % 1st derivative of log potential, exp(-s)/(1+exp(-s))
  P(:,3) = p.*(p-1);   % 2nd derivative of log potential, -exp(-s)/(1+exp(-s))^2
  P(:,4) = ones(size(s))/2;
else
  %fprintf(1,'PotMSErf: type=%s\n',type);
  if strcmp(type,'VB')
    % Otherwise it does not work!
    P = potLogistic(s);
  elseif strcmp(type,'EP')
    q = numel(s); P = zeros(q,3); % allocate memory
    % Use GPML 'likErf' to implement this here
    [lZ,dlZ,d2lZ] = likErf([],ones(q,1),s(:),z(:),'infEP');
    P(:,1) = lZ; P(:,2) = dlZ; P(:,3) = d2lZ;
  else
    error('Unknown type')
  end
end
