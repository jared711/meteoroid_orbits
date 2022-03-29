function [x_dim, A] = dim(x_nondim, DU, VU)
%DIM Summary of this function goes here
% 
% [OUTPUTARGS] = DIM(INPUTARGS)
% 
% Inputs: 
% 
% Outputs: 
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/02/23 13:49:01 	Revision: 0.1 $

if ~iscolumn(x_nondim);    x_nondim = x_nondim'; end

A = [eye(3)*DU,  zeros(3);
      zeros(3), eye(3)*VU];
x_dim = A*x_nondim;

end
