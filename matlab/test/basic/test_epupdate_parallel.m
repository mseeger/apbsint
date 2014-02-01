% Test script for EPTOOLS_EPUPDATE_PARALLEL:
% Construct potential manager with Laplace, Probit and Gaussian
% potentials. We test against GPML code.
%
% See EPT.TEST_EPUPDATE_SINGLE for details.

num_it = 500;

% Create potential manager and cavity marginals. This is done in
% the same way as in EPT.TEST_EPUPDATE_SINGLE
clear('pman');
cmu=zeros(3*num_it,1); crho=zeros(3*num_it,1);
pman{1}.name = 'Gaussian';
pman{1}.size = num_it;
pman{1}.pars{1}=randn(num_it,1);
pman{1}.pars{2}=exp(0.5*randn(num_it,1));
pman{2}.name = 'Probit';
pman{2}.size = num_it;
pman{2}.pars{1}=2*(rand(num_it,1)>0.5)-1;
pman{2}.pars{2}=[0];
cmu(1:2*num_it) = 5*randn(2*num_it,1);
crho(1:2*num_it) = exp(1.5*randn(2*num_it,1)-1);
pman{3}.name = 'Laplace';
pman{3}.size = num_it;
pman{3}.pars{1}=[0];
pman{3}.pars{2}=rand(num_it,1)+0.05;
cmu(2*num_it+1:end) = 2*randn(num_it,1);
crho(2*num_it+1:end) = exp(randn(num_it,1));
% Create internal representation
irp = ept.potman_intrepres(pman);

% Run EPTOOLS_EPUPDATE_PARALLEL
[rstat,alpha1,nu1,logz1] = ...
    eptools_epupdate_parallel(irp.potids,irp.numpot,irp.parvec, ...
			      irp.parshrd,cmu,crho);
fprintf(1,'Number of numerical failures: %d\n',sum(~rstat));

% Comparison (GPML)
alpha2 = zeros(3*num_it,1); nu2 = zeros(3*num_it,1);
logz2 = zeros(3*num_it,1);
off = 0; rng = off+1:off+num_it;
y = pman{1}.pars{1}; ssq = pman{1}.pars{2};
nu2(rng) = 1./(crho(rng)+ssq);
tmp = y-cmu(rng);
alpha2(rng) = nu2(rng).*tmp;
logz2(rng) = -0.5*(nu2(rng).*(tmp.^2)+log(2*pi)-log(nu2(rng)));
off = num_it; rng = off+1:off+num_it;
y = pman{2}.pars{1};
[logz2(rng),alpha2(rng),mnu2] = likErf([],y,cmu(rng),crho(rng), ...
				       'infEP');
nu2(rng) = -mnu2;
off = 2*num_it; rng = off+1:off+num_it;
hyp = 0.5*log(2)-log(pman{3}.pars{2});
y = zeros(num_it,1);
[logz2(rng),alpha2(rng),mnu2] = likLaplace(hyp,y,cmu(rng), ...
					   crho(rng),'infEP');
nu2(rng) = -mnu2;

% Print largest errors
fstr{1} = ['%s: df=%f (1: %f; 2: %f), cmu=%f,crho=%f, ' ...
	   'y=%f,ssq=%f\n'];
fstr{2} = ['%s: df=%f (1: %f; 2: %f), cmu=%f,crho=%f, ' ...
	   'y=%f,      %d\n'];
fstr{3} = ['%s: df=%f (1: %f; 2: %f), cmu=%f,crho=%f, ' ...
	   'tau=%f,      %d\n'];
off=0;
for typ = 1:3
  rng = off+1:off+num_it;
  fprintf(1,'\n%s: Largest errors:\n',pman{typ}.name);
  v1 = nu1(rng); v2 = nu2(rng);
  err = ept.reldiff(v1,v2);
  [mx,i] = max(err);
  if typ~=3
    p1 = pman{typ}.pars{1}(i);
  else
    p1 = pman{3}.pars{2}(i);
  end
  if typ==1
    p2 = pman{typ}.pars{2}(i);
  else
    p2 = 0;
  end
  fprintf(1,fstr{typ},'nu',mx,v1(i),v2(i),cmu(off+i),crho(off+i), ...
	  p1,p2);
  v1 = alpha1(rng); v2 = alpha2(rng);
  err = ept.reldiff(v1,v2);
  [mx,i] = max(err);
  if typ~=3
    p1 = pman{typ}.pars{1}(i);
  else
    p1 = pman{3}.pars{2}(i);
  end
  if typ==1
    p2 = pman{typ}.pars{2}(i);
  else
    p2 = 0;
  end
  fprintf(1,fstr{typ},'alpha',mx,v1(i),v2(i),cmu(off+i),crho(off+i), ...
	  p1,p2);
  v1 = logz1(rng); v2 = logz2(rng);
  err = ept.reldiff(v1,v2);
  [mx,i] = max(err);
  if typ~=3
    p1 = pman{typ}.pars{1}(i);
  else
    p1 = pman{3}.pars{2}(i);
  end
  if typ==1
    p2 = pman{typ}.pars{2}(i);
  else
    p2 = 0;
  end
  fprintf(1,fstr{typ},'logz',mx,v1(i),v2(i),cmu(off+i),crho(off+i), ...
	  p1,p2);
  if typ<3
    pause;
  end
  off = off+num_it;
end
