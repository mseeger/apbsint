function print_python_mat(a)
%PRINT_PYTHON_MAT Prints matrix for insertion into Python script

[m,n] = size(a);
if m==1 || n==1
  a = a(:)';
  [m,n] = size(a);
  isv = 1;
else
  isv = 0;
  fprintf('[');
end
for j=1:m
  fprintf('[');
  fprintf('%.25x, ',a(j,1:end-1));
  if j<m
    fprintf('%.25x], \\\n ',a(j,end));
  else
    fprintf('%.25x]',a(j,end));
  end
end
if ~isv
  fprintf(']');
end
