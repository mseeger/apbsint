%EPTOOLS: Data type for potential manager
%  A potential manager object PMAN is either a struct or a cell
%  array of structs. In latter case, each entry is called a
%  block. In the former case, there is just one block. The PMs for
%  blocks are appended in order.
%
%  Each block has a potential data type, number of potentials, and
%  potential parameters. If PMAN is a block, PMAN.NAME is the type
%  name (f.ex. 'Gaussian', 'Laplace'), see
%  EPTOOLS_GETPOTID. PMAN.SIZE is the number of
%  potentials. PMAN.PARS is a cell array of vectors, one entry for
%  each potential attribute (number and restrictions depend on
%  potential type). An entry can either be a vector of size
%  PMAN.SIZE or a scalar (shared value).
