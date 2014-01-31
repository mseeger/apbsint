function debug_printvec(vec,fid)
%DEBUG_PRINTVEC Print vector into file (text format)

if isa(vec,'int32')
  estr = '%d';
else
  estr = '%f';
end
fprintf(fid,[estr ' '],vec(1:end-1));
fprintf(fid,[estr '\n'],vec(end));
