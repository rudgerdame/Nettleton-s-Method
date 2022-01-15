function glmcosi = Topo2Grav(Hlmcosi,drho,RH,type,nmax)
% INPUTS:
% Hlmcosi - coefficients for the interface relief (topography, Moho, etc.)
% drho    - the density contrast across that interface
% RH      - the mean radius of that interface
% type    - type of gravity (1=potential, 2=geoid, 3=gravity disturbance, 
%           4=free air gravity)
% nmax    - Maximum power
%
% OUTPUT:
% glmcosi - spherical harmonic coefficients of the expected gravity
%
% Created by P_James@baylor.edu, 9/17/19

G = 6.67408 * 10^-11;
degrees = Hlmcosi(:,1);
lmax = max(degrees);
dres = 180/lmax;

glmcosi = Hlmcosi; glmcosi(:,3:4)=0;

% For n=1 (the same as a mass sheet)
if type==1 % Gravitational potential
    for i=2:length(degrees)
        l=degrees(i);
        glmcosi(i,3) = 4*pi*G*RH /(2*l+1) * drho * Hlmcosi(i,3);
        glmcosi(i,4) = 4*pi*G*RH /(2*l+1) * drho * Hlmcosi(i,4);        
    end    
elseif type==4 % Free-air gravity
    for i=2:length(degrees)
        l=degrees(i);
        glmcosi(i,3) = 4*pi*G * (l-1)/(2*l+1) * drho * Hlmcosi(i,3);
        glmcosi(i,4) = 4*pi*G * (l-1)/(2*l+1) * drho * Hlmcosi(i,4);        
    end
end

% For n>1

if nmax>1
    for n=2:nmax 
        if type==1 % Gravitational potential
            Hmap = plm2xyz(Hlmcosi,dres);
            nHrho = Hmap.^n*drho;
            nHlmcosi = xyz2plm(nHrho,lmax);
            for i=2:length(degrees)
                l=degrees(i);
                prod=1/(l+3);
                for j=1:n
                    prod = prod * (l+4-j)/j;
                end
                glmcosi(i,3) = glmcosi(i,3) + 4*pi*G*RH^2 /(2*l+1) * ...
                    nHlmcosi(i,3)/ (RH^n) * prod;
                glmcosi(i,4) = glmcosi(i,4) + 4*pi*G*RH^2 /(2*l+1) * ...
                    nHlmcosi(i,4)/ (RH^n) * prod;                
            end
        elseif type==4 % Free air gravity
            Hmap = plm2xyz(Hlmcosi,dres);
            nHrho = Hmap.^n*drho;
            nHlmcosi = xyz2plm(nHrho,lmax);
            for i=2:length(degrees)
                l=degrees(i);
                prod=1/(l+3);
                for j=1:n
                    prod = prod * (l+4-j)/j;
                end
                glmcosi(i,3) = glmcosi(i,3) + 4*pi*G*RH *(l-1)/(2*l+1) * ...
                    nHlmcosi(i,3)/ (RH^n) * prod;
                glmcosi(i,4) = glmcosi(i,4) + 4*pi*G*RH *(l-1)/(2*l+1) * ...
                    nHlmcosi(i,4)/ (RH^n) * prod;                
            end
        end
    end
end

end



