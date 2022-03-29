function [a_std, e_std, i_std, Omega_std, omega_std, ME_std] = oe_err(x, GM, Sigma)
%OE_ERR Montenbruck pg 238
% 
% [OUTPUTARGS] = OE_ERR(INPUTARGS)
% 
% Inputs: 
% 
% Outputs: 
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/02/21 10:07:45 	Revision: 0.1 $
if ~iscolumn(x);    x = x'; end

rVec = x(1:3);
r = norm(rVec);
vVec = x(4:6);
v = norm(vVec);

hVec = cross(rVec, vVec);
nVec = cross([0; 0; 1], hVec); % 3/28/22 I think I can get rid of this now, but I need to verify what changes first
eVec = (1/GM)*((v^2 - GM/r)*rVec - dot(rVec, vVec)*vVec); % verified from D'Amico notes 3/28/22
qVec = cross(hVec,eVec);

P = unit(eVec);
Q = unit(qVec);
W = unit(hVec);

[a,e,i,Omega,omega,nu] = cart2oe(x,GM);
E = nu2E(nu,e);
n = sqrt(GM/a^3); % mean motion
xhat = a*(cosd(E) - e);
xhatdot = -sqrt(GM*a)*sind(E)/r;
yhat = a*sqrt(1-e^2)*sind(E);
yhatdot = sqrt(GM*a)*sqrt(1-e^2)*cosd(E)/r;

% a, e, M
dxhatda = xhat/a;
dyhatda = yhat/a;
dxhatdotda = -xhatdot/(2*a);
dyhatdotda = -yhatdot/(2*a);
drda = dxhatda*P + dyhatda*Q;
dvda = dxhatdotda*P + dyhatdotda*Q;

dxhatde = -a - yhat^2/(r*(1-e^2));
dyhatde = xhat*yhat/(r*(1-e^2));
dxhatdotde = xhatdot*(a/r)^2*(2*xhat/a + e/(1-e^2)*(yhat/a)^2);
dyhatdotde = n/sqrt(1-e^2)*(a/r)^2*(xhat^2/r - yhat^2/(a*(1-e^2)));
drde = dxhatde*P + dyhatde*Q;
dvde = dxhatdotde*P + dyhatdotde*Q;

dxhatdM = xhatdot/n;
dyhatdM = yhatdot/n;
dxhatdotdM = -n*(a/r)^3*xhat;
dyhatdotdM = -n*(a/r)^3*yhat;
drdM = dxhatdM*P + dyhatdM*Q;
dvdM = dxhatdotdM*P + dyhatdotdM*Q;

% i, Omega, omega
dPdi = sind(omega)*W;
dQdi = cosd(omega)*W;
drdi = xhat*dPdi + yhat*dQdi;
dvdi = xhatdot*dPdi + yhatdot*dQdi;

dPdOmega = [-P(2);P(1);0];
dQdOmega = [-Q(2);Q(1);0];
drdOmega = xhat*dPdOmega + yhat*dQdOmega;
dvdOmega = xhatdot*dPdOmega + yhatdot*dQdOmega;

dPdomega = Q;
dQdomega = -P;
drdomega = xhat*dPdomega + yhat*dQdomega;
dvdomega = xhatdot*dPdomega + yhatdot*dQdomega;

% I don't know why this isn't equivalent 3/3/22
% dxdalpha = [drda, drde, drdi, drdOmega, drdomega, drdM;
%             dvda, dvde, dvdi, dvdOmega, dvdomega, dvdM];
% dalphadx = pinv(dxdalpha);
        
PP = zeros(6);
PP(1,6) = -2/(n*a);
PP(6,1) = -PP(1,6);
PP(2,5) = sqrt(1-e^2)/(n*a^2*e);
PP(5,2) = -PP(2,5);
PP(2,6) = -(1-e^2)/(n*a^2*e);
PP(6,2) = -PP(2,6);
PP(3,4) = 1/(n*a^2*sqrt(1-e^2)*sind(i));
PP(4,3) = -PP(3,4);
PP(3,5) = -cotd(i)/(n*a^2*sqrt(1-e^2));
PP(5,3) = -PP(3,5);


dalphadx = PP*[[dvda, dvde, dvdi, dvdOmega, dvdomega, dvdM]', -[drda, drde, drdi, drdOmega, drdomega drdM]']; %(6x6)


J_ME = [GM/r^3*rVec', vVec'];
ME_std = sqrt(J_ME*Sigma*J_ME');

Sigma = dalphadx*Sigma*dalphadx';

a_std = sqrt(Sigma(1,1));
e_std = sqrt(Sigma(2,2));
i_std = rad2deg(sqrt(Sigma(3,3)));
Omega_std = rad2deg(sqrt(Sigma(4,4)));
omega_std = rad2deg(sqrt(Sigma(5,5)));


end
