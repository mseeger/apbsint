%TEST_CHOLUPDNRK1 Test script for Cholesky update/downdate
%  Small example. Used to compare Python against Matlab implementation

n = 10; r = 2;
print_python = 1;

lfact = tril(randn(n));
lfact(1:n+1:end) = 0.1+rand(n,1);
amat = lfact*lfact';
zmat = randn(r,n);
vec = randn(n,1);
yvec = randn(r,1);
cvec = zeros(n,1); svec = zeros(n,1);
workv = zeros(max(n,r),1);

% CHOLUPRK1
lfact1 = lfact(:,:); % Real copy!
zmat1 = zmat(:,:);
stat = eptools_choluprk1({lfact1,[1 1 n n],'LN'},vec,cvec,svec, ...
			 workv,zmat1,yvec);
if stat~=0
  error('EPTOOLS_CHOLUPRK1: Numerical error!');
end
lfact1b = chol(amat+vec*vec','lower');
zmat1b = (lfact1b\(lfact*zmat'+vec*yvec'))';
fprintf('CHOLUPRK1: rdf(L)=%f, rdf(Z)=%f\n', ...
	max(ept.reldiff(lfact1,lfact1b)),max(ept.reldiff(zmat1, ...
						  zmat1b)));

% CHOLDNRK1
lfact2 = lfact1(:,:);
zmat2 = zmat1(:,:);
pvec = lfact1\vec;
stat = eptools_choldnrk1({lfact2,[1 1 n n],'LN'},pvec,cvec,svec, ...
			 workv,1,zmat2,yvec);
fprintf('CHOLDNRK1: rdf(L)=%f, rdf(Z)=%f\n', ...
	max(ept.reldiff(lfact2,lfact)),max(ept.reldiff(zmat2, ...
						  zmat)));

if print_python
  % Outputs for Python comparison
  fprintf('lfact = np.array(');
  ept.print_python_mat(lfact); fprintf(',order=''F'')\n');
  fprintf('zmat = np.array(');
  ept.print_python_mat(zmat); fprintf(',order=''F'')\n');
  fprintf('vec = np.array(');
  ept.print_python_mat(vec); fprintf(')\n');
  fprintf('yvec = np.array(');
  ept.print_python_mat(yvec); fprintf(')\n');
  fprintf('pvec = np.array(');
  ept.print_python_mat(pvec); fprintf(')\n');
  fprintf('lfact1 = np.array(');
  ept.print_python_mat(lfact1); fprintf(',order=''F'')\n');
  fprintf('zmat1 = np.array(');
  ept.print_python_mat(zmat1); fprintf(',order=''F'')\n');
  fprintf('lfact2 = np.array(');
  ept.print_python_mat(lfact2); fprintf(',order=''F'')\n');
  fprintf('zmat2 = np.array(');
  ept.print_python_mat(zmat2); fprintf(',order=''F'')\n');
end
