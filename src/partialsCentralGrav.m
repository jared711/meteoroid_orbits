function [dadr, dadv] = partialsCentralGrav(r, gm)
%PARTIALSCENTRALGRAV Summary of this function goes here
% 
% [OUTPUTARGS] = PARTIALSCENTRALGRAV(INPUTARGS)
% 
% Inputs: 
% 
% Outputs: 
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/02/22 08:25:03 	Revision: 0.1 $

dadv = zeros(3);
dadr = gm/(norm(r)^5)*...
       [3*r(1)^2 - norm(r)^2, 3*r(1)*r(2), 3*r(1)*r(3);
        3*r(1)*r(2), 3*r(2)^2 - norm(r)^2, 3*r(2)*r(3);
        3*r(1)*r(3), 3*r(2)*r(3), 3*r(3)^2 - norm(r)^2];
  

end
