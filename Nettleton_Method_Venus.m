% Nettleton Method on Venus
%% Path setup
clear all; close all;
% set up path to function folders
current_folder = pwd; func = append(current_folder,'/functions'); 
path(func,path)

%% Study Area
% Study Area: Input Venus longitude and latitude of an location that you wish to calculate
%       the bulk density of the crust.
Longitude = 80;
Latitude = 25;

%% Topography of study area
RH = 6051.880596000000E+03;
Hlmcosi = load('VenusData/shtjv360.a02');
lmax = 180;
Hlmcosi = Hlmcosi(1:addmup(lmax),:);
Hlmcosi(:,3:4)=Hlmcosi(:,3:4)*RH;

% 10x10 degree study area
lon_max = Longitude + 5;   
lon_min = Longitude - 5;
lat_max = Latitude + 5;
lat_min = Latitude - 5;

c11cmn = [lon_min lat_max lon_max lat_min];
dres = .1;
Hmap=plm2xyz(Hlmcosi,dres,c11cmn);

% Plot Topography of Study Area
figure
imagesc(lon_min:dres:lon_max,lat_max:-dres:lat_min,(Hmap+RH)*10^-3)
set(gca,'YDir','normal','FontSize',20)
axis square
c = colorbar;
c.Label.String = 'km';
set(c,'FontSize',20)
title('Venus SH Topography','FontSize',30)
MeanTopo = mean(mean(Hmap));

%% Finite Amplitude Gravity from Topography
drho = 1;
type = 4;
nmax = 4;

% Range of Spherical Harmonic Degrees
lmax = 80;
lmin=40;

gftlmcosi = Topo2Grav(Hlmcosi, drho, RH, type, nmax);
gftlmcosi = gftlmcosi(1:addmup(lmax),:);
gftlmcosi(1:addmup(lmin-1),3:4) = 0;

%% Load Magellan gravity
filename = 'VenusData/shgj180u.a01';
[glmcosi, GM, R] = read_TAB(filename);
R = .6051000000000000E+04 * 10^3;           %From shgj180u.a01 file since read_TAB function 
                                            %doesn't read this file correctly.                               
GM = .3248585920790000E+06 * 10^9;          %from shgj180u.a01 file ''
glmcosi = glmcosi(1:addmup(lmax),:);
glmcosi(1:addmup(lmin-1),3:4) = 0;
type = 4;
glmcosi_FA = clm2grav(glmcosi,type,GM,R);

%% Upward Continue Gravity to mean Topography
type = 4;

glmcosi = GravityContinuation(glmcosi_FA,R,RH+MeanTopo,type);
gftlmcosi = GravityContinuation(gftlmcosi,RH,RH+MeanTopo,type);

%% Gravity from the Moho
rho_over_drho = 2500/500;
drho_over_rho = 1/rho_over_drho;
D = 16e3;     %Haastse-baad Tessera
RW = RH-D;

Wlmcosi = Hlmcosi;
Wlmcosi(:,3:4) = -rho_over_drho * Hlmcosi(:,3:4);
Wlmcosi = Wlmcosi(1:addmup(lmax),:);

nmax = 4;

GFWlmcosi = Topo2Grav(Wlmcosi,drho_over_rho,RW,type,nmax);
GFWlmcosi = GFWlmcosi(1:addmup(lmax),:);
GFWlmcosi(1:addmup(lmin-1),3:4) = 0;

FAFWlmcosi = GravityContinuation(GFWlmcosi,RW,RH+MeanTopo,type);
[gfwmap,lon,lat,Plm] = plm2xyz(FAFWlmcosi,dres,c11cmn);

%% Plot of Free Air Gravity, Topography and MOHO 
[gmap,lon,lat,Plm] = plm2xyz(glmcosi,dres,c11cmn);
[gftmap,lon,lat] = plm2xyz(gftlmcosi,dres,c11cmn);
figure
imagesc(lon_min:dres:lon_max,lat_max:-dres:lat_min,gmap*10^5)
set(gca,'YDir','normal','FontSize',20)
axis square
c = colorbar;
c.Label.String = 'mGal';
set(c,'FontSize',20)
title('Free Air Gravity','FontSize',30)
figure
imagesc(lon_min:dres:lon_max,lat_max:-dres:lat_min,gftmap*10^5)
set(gca,'YDir','normal','FontSize',20)
axis square
c = colorbar;
c.Label.String = 'mGal';
set(c,'FontSize',20)
title('Gravity from Topography (Assuming p = 1 kg/m^2)','FontSize',30)
figure
imagesc(lon_min:dres:lon_max,lat_max:-dres:lat_min,gfwmap*10^5)
set(gca,'YDir','normal','FontSize',20)
axis square
c = colorbar;
c.Label.String = 'mGal';
set(c,'FontSize',20)
title('Gravity of the Moho','FontSize',30)

%% Density Calculation
Gvec = mat2vec(gmap);
Tvec = mat2vec(gftmap+gfwmap);
A =  [Tvec ones(size(Tvec))];
b = Gvec;
x = A\b;
BulkDensity = x(1);

str1 = 'The bulk density of the area is: ';  str2 = ' kg/m^2'; 
str = append(str1,num2str(BulkDensity),str2);
disp(' ')
disp(str)