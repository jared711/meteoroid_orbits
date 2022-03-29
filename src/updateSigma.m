function [Sigma2] = updateSigma(Sigma1, R)
%UPDATESIGMA Summary of this function goes here
% 
% [Sigma2] = UPDATESIGMA(Sigma1, R)
% 
% Inputs:   Sigma1 [] (6x6) initial covariance matrix
%           R [] (6x6) Jacobian Matrix
% 
% Outputs:  Sigma2 [] (6x6) updated covariance matrix
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/02/21 10:20:20 	Revision: 0.1 $

Sigma2 = R*Sigma1*R';

end
