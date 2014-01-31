% Test script for EPTOOLS_EPUPDATE_SINGLE:
% For now, we test Laplace and Probit against GPML code,checking
% whether they approximately agree.
%
% GPML EP code:
%   [logz,alpha,mnu] = lik(hyp,y,cmu,crho,'infEP')
% Here, mnu is -nu, hyp are hyperparameters, y the target, cmu,
% crho cavity mean, variance.
%
% (1) Probit:
% GPML: likErf, hyp=[], y in {+1,-1}
%
% (2) Laplace:
% GPML: likLaplace, hyp=[log(b*sqrt(2))]
%   likLaplace(t) = exp(-|t-y|/b)/(2*b)
% We use tau=1/b, so hyp=[log(sqrt(2)/tau)]

num_it = 500;

nam=[]; pid=[]; rstr=[];
nam{1}='Probit';
nam{2}='Laplace';
pid{1} = eptools_getpotid(nam{1});
pid{2} = eptools_getpotid(nam{2});
rstr{1} = ['%d: y=%f, cmu=%f,crho=%f.\n' ...
	   '  nu:    df=%f (1: %f; 2: %f)\n' ...
	   '  alpha: df=%f (1: %f; 2: %f)\n' ...
	   '  logz:  df=%f (1: %f; 2: %f)\n'];
rstr{2} = ['%d: y=%f, tau=%f, cmu=%f,crho=%f.\n' ...
	   '  nu:    df=%f (1: %f; 2: %f)\n' ...
	   '  alpha: df=%f (1: %f; 2: %f)\n' ...
	   '  logz:  df=%f (1: %f; 2: %f)\n'];
for typ = 1:2
  mxnu    = [0];
  mxalpha = [0];
  mxlogz  = [0];
  fprintf(1,'\nTests for %d:\n',nam{typ});
  for it = 1:num_it
    if typ==1
      y = 2*(rand>0.5)-1;
      tau = 0; % Dummy!
      ppars = [y 0];
      hyp = [];
      lkfun = @likErf;
      cmu = 5*randn;
      crho = exp(1.5*randn-1);
    else
      y = 0; % Does not matter anyway
      %tau = exp(1-1.5*randn);
      tau = 0.05+rand;
      ppars = [y tau];
      hyp = [0.5*log(2)-log(tau)];
      lkfun = @likLaplace;
      cmu = 2*randn;
      crho = exp(randn);
    end
    [rstat,alpha1,nu1,logz1] = eptools_epupdate_single(pid{typ}, ...
						  ppars,cmu,crho);
    [logz2,alpha2,mnu2] = feval(lkfun,hyp,y,cmu,crho,'infEP');
    nu2=-mnu2;
    dnu=ept.reldiff(nu1,nu2);
    if dnu>mxnu(1)
      mxnu = [dnu nu1 nu2 it y cmu crho tau];
    end
    dalpha=ept.reldiff(alpha1,alpha2);
    if dalpha>mxalpha(1)
      mxalpha = [dalpha alpha1 alpha2 it y cmu crho tau];
    end
    dlogz=ept.reldiff(logz1,logz2);
    if dlogz>mxlogz(1)
      mxlogz = [dlogz logz1 logz2 it y cmu crho tau];
    end
    if typ==1
      fprintf(1,rstr{1},it,y,cmu,crho,dnu,nu1,nu2,dalpha,alpha1, ...
	      alpha2,dlogz,logz1,logz2);
    else
      fprintf(1,rstr{2},it,y,tau,cmu,crho,dnu,nu1,nu2,dalpha,alpha1, ...
	      alpha2,dlogz,logz1,logz2);
    end
  end
  fprintf(1,'\n%s: Largest errors:\n',nam{typ});
  mx=mxnu;
  fprintf(1,'nu:    df=%f (1: %f; 2: %f), cmu=%f,crho=%f,y=%f,tau=%f\n', ...
	  mx(1),mx(2),mx(3),mx(6),mx(7),mx(5),mx(8));
  mx=mxalpha;
  fprintf(1,'alpha: df=%f (1: %f; 2: %f), cmu=%f,crho=%f,y=%f,tau=%f\n', ...
	  mx(1),mx(2),mx(3),mx(6),mx(7),mx(5),mx(8));
  mx=mxlogz;
  fprintf(1,'logz:  df=%f (1: %f; 2: %f), cmu=%f,crho=%f,y=%f,tau=%f\n', ...
	  mx(1),mx(2),mx(3),mx(6),mx(7),mx(5),mx(8));
  if typ==1
    pause;
  end
end
