% EPTOOLS Matlab Interface
% Startup script. Mainly sets PATH.
%
% Adapted frmo glm-ie_v1.5 startup.m (Hannes Nickisch)

disp(['ApBsInT startup script...']);

me = mfilename; % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me)); % where am I?
addpath(mydir(1:end-1))
addpath([mydir,'bin']); % MEX libraries
addpath([mydir,'mex']); % Help for MEX functions
addpath([mydir,'test'])
addpath([mydir,'test/basic'])
addpath([mydir,'test/binclass'])

clear me mydir
