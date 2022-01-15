function [lmcosi_out] = clm2grav(lmcosi_in, type, GM, R)
% INPUTS:
% lmcosi_in - Dimentionless SH coefficients
% type      - Type of grav anomaly:
%              1 - Potential (m^2/s^2)
%              2 - Geoid (meters)
%              3 - Gravity Disturbance (m/s^2)
%              4 - Free air (m/s^2)
% GM        - gravitational parameter
% R         - Radius
%
% OUTPUT:
% lmcosi_out - Gravity coefficients
%Created by rudger_dame1@baylor.edu, 9/5/19

degrees = lmcosi_in(:,1);
lmcosi_out = lmcosi_in;
ncols = size(lmcosi_in, 2);


if type == 1
    % Calculate potential
    lmcosi_out(:,3:ncols) = GM/R *lmcosi_in(:,3:ncols);
    
elseif type == 2
    % Geoid
    lmcosi_out(:,3:ncols) = R *lmcosi_in(:,3:ncols);
    
elseif type == 3
    % Gravity Disturbance
    factor1 = GM/R^2 * (degrees + 1);   
    lmcosi_out(:, 3:ncols) = repmat(factor1,1,4) .* lmcosi_in(:, 3:ncols);
    
elseif type == 4
    %Free air
    factor2 = GM/R^2 * (degrees - 1);   
    lmcosi_out(:, 3:ncols) = repmat(factor2,1,4) .* lmcosi_in(:, 3:ncols);
    
end 

end