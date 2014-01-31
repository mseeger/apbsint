% Modification M. Seeger:
% This is to compare against our EPTEST_BINCLASS_COUP.M
% - Gaussian prior (also Laplace; optional)
% - Parallel updating EP, with a single IL iteration
% - Output same delta
% - Uses logit likelihood, while we use probit
%
% Inference example for logistic regression using glm-ie demonstrating the 
% sparsity occuring in estimation but similar classification performance in
% inference where probabilities can readily be evaluated.
%
% We use the well known a9a dataset from the UCI machine learning repository
% with 32561 instances of dimension 123.
% B data matrix
% c binary labels
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 30
clear all, close all

global deb_glm_pi deb_glm_beta deb_glm_h deb_glm_rho; % DEBUG!
%global deb_ept_pi deb_ept_beta deb_ept_h deb_ept_rho; % DEBUG!

%load deb_ept % DEBUG

doLaplace = 0; % Laplace prior? Otherwise: Gaussian

potS = @potLaplace; % Laplace weight prior
%potLh = @potLogistic;
potLh = @potMSErf;

% load and split into train and test set
load classify                               % load a part of the UCI a9a dataset
nte = 30000;                        % say how many test examples we wich to keep
Bte = B(1:nte,:); B = B(nte+1:end,:); cte = c(1:nte); c = c(nte+1:end);  % split
% Remove variables not touched by any observations
ind = find(sum(B~=0,1));
B = B(:,ind); Bte = Bte(:,ind);
[q,n] = size(B);

% nearest neighbour classification
dp = sum(bsxfun(@minus, Bte, full(mean(B(c==+1,:)))).^2,2);  % dist. to pos. ctr
dm = sum(bsxfun(@minus, Bte, full(mean(B(c==-1,:)))).^2,2);  % dist. to neg. ctr
cc = 2*double(dp<dm)-1;
fprintf('nearest neighbor accuracy=%1.2f%%\n',100*sum(cte==cc)/numel(cte))

% classification
%fprintf('Do regression using a ')
%str = func2str(potS);
%fprintf('%s',str(8:end))
if doLaplace
  % Laplace prior
  m = 0; X = zeros(0,n); y = zeros(0,1); s2 = 1; % no Gaussian part
  %tau = [(2/5)*ones(n,1); c.*ones(q,1)];
  tau = [2*ones(n,1); 5*c.*ones(q,1)];
  B = [eye(n); B]; t = 0;
  pot = @(s,varargin) potCat(s, varargin{:}, {potS,potLh}, ...
			     {1:n,n+(1:q)} );
else
  % Gaussian prior
  m = n; X = eye(n); y = zeros(n,1); s2 = 25/4; % Gaussian part
  tau = c(:); t = 0;
  pot = @(s,varargin) potCat(s, varargin{:},{potLh},{1:q});
end

%fprintf(' weight prior.\n')

%% Inference
opts.innerType = 'EP';
opts.innerOutput = 1;
opts.outerOutput = 1;
opts.outerNiter  = 20;        % Number of outer loop iterations
opts.outerMethod = 'full';    % Exact variances
opts.innerExact = 1;          % Exact means
%opts.innerMVM   =  50;        % CG steps
%opts.innerVBpls = 'plsTN';    % PLS algorithm, use LBFGS if compiled
[uinf,ga,b,z,zu,nlZ] = ms_dli(X,y,s2,B,t,pot,tau,opts);
sdinf = sqrt(zu);

%% Estimation
%opt.nMVM = 150; opt.output = 1;
%u0 = zeros(n,1);
%uest = feval(opts.innerVBpls,u0,X,y,B,t,opt,s2,'penVB',pot,tau/sqrt(s2));

fprintf('accuracies on test set\n')
acc_inference  = 100*sum(cte==sign(Bte*uinf))/numel(cte)
%acc_estimation = 100*sum(cte==sign(Bte*uest))/numel(cte)

% ci-f1.png
%sz = [500,300]; figure('Position',[50,50,sz]), set(gca,'FontSize',16)
%plot([-1,1],[-1,1],'r'), hold on
%for i=1:length(uinf)
%  plot(uinf(i)+sdinf*[-1,1],uest(i)*[1,1],'c')
%end
%plot(uinf,uest,'bo')
%plot(uinf,uest,'k.')
%axis([-1,1,-1,1])
%xlabel('u_{VB/EP}'), ylabel('u_{MAP}')
%title('sparsity pattern: inference vs. estimation')
