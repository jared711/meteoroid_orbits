function [D] = oe_similar(oeA, oeB)
%OE_SIMILAR Similarity function between two orbits as described by orbital
%elements. 
% Southworth, R. B., & Hawkins, G. S. (1963). Statistics of meteor streams. Smithsonian Contributions to Astrophysics, 261–285.
% 
% [D] = OE_SIMILAR(oeA, oeB)
% 
% Inputs:   as [2x1] (AU)
% 
% Outputs: 
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/01/21 17:37:31 	Revision: 0.1 $


aA = oeA(1);    aB = oeB(1);
eA = oeA(2);    eB = oeB(2);
iA = oeA(3);    iB = oeB(3);
oA = oeA(4);    oB = oeB(4);
uA = oeA(5);    uB = oeB(5);
qA = aA*(1-eA); qB = aB*(1-eB);

iA = deg2rad(iA);   iB = deg2rad(iB);
oA = deg2rad(oA);   oB = deg2rad(oB);
uA = deg2rad(uA);   uB = deg2rad(uB);


D = (eB - eA)^2 + (qB - qA)^2 + (2*sin((iB-iA)/2))^2 +...
    sin(iA)*sin(iB)*(2*sin((oB-oA)/2))^2 +...
    ((eA+eB)/2*2*sin((oB+uB-oA-uA)/2))^2;

end
