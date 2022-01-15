function lmcosi_out = GravityContinuation(lmcosi_in, R_in, R_out, type)
% INPUTS:
% lmcosi_in - Input spherical harmonic coefficients of gravity
% type      - Type of gravity anomaly (for both input and output):
%             1 - Potential (m^2/s^2)
%             2 - Geoid (meters)
%             3 - Gravity disturbance (m/s^2)
%             4 - Free air gravity (m/s^2)
% R_in      - the radius used for lmcosi_in
% R_out     - The radius where you want new gravity coefficients
%
% OUTPUT:
% lmcosi_out
%
% Created by rhdame@gmail.com, 9/10/19

degrees = lmcosi_in(:,1);
lmcosi_out = lmcosi_in;
ncols = size(lmcosi_in,2);

if type==1 % Potential
    factor = (R_in/R_out) .^ (degrees+1);
    lmcosi_out(:,3:ncols) = lmcosi_in(:,3:ncols).*repmat(factor,1,ncols-2);
elseif type==3 || type==4 % Gravity acceleration
    factor = (R_in/R_out) .^ (degrees+2);
    lmcosi_out(:,3:ncols) = lmcosi_in(:,3:ncols).*repmat(factor,1,ncols-2);

end


end