function [x] = make_column_vecs(x)
%MAKE_COLUMN_VECS Summary of this function goes here
% 
% [x] = MAKE_COLUMN_VECS(x)
% 
% Inputs:   x [] (6XN OR NX6) matrix of vectors  
% 
% Outputs:  x [] (6XN) matrix of vectors
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/03/28 09:06:15 	Revision: 0.1 $

[m,n] = size(x);
if m == 6 || m == 3
    if n == 6 || n == 3
        warning('Dimensions of x_dim are ambiguous, leaving as is');
    end
elseif n == 6 || n == 3
    x = x';
else
    warning('Dimensions of x_dim are likely incorrect, leaving as is')
end

end
