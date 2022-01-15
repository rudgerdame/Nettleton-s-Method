function [lmcosi,GM,R]=read_TAB(filename)
% function [lmcosi,GM,R]=READ_TAB(filename)
% Reads a file of spherical harmonic coefficients in the ".tab" format

fid=fopen(filename);
if fid==-1
  error('File not found or permission denied');
end

% disp('Make sure to check that the first line starts with one space')
line=fgetl(fid);
Rstr=line(1:19);
Rexp=line(21:23);
GMstr=line(26:43); 
GMexp=line(45:47);

R=str2double(Rstr)*10^str2double(Rexp)*10^3; % In meters!
GM=str2double(GMstr)*10^str2double(GMexp)*10^9; % m^3/s^2

fclose(fid);
lmcosi = csvread(filename, 1, 0);
try lmcosi2 = dlmread(filename,'\t', 1, 0);
catch
    lmcosi2 = 0;
end

% Deal with tab delimited values...
if size(lmcosi2,2)>size(lmcosi,2)
    lmcosi=lmcosi2;
end

% [header,lmcosi]=hdrload(filename);
% Rstr=header(2:19); Rexp=header(22:23);
% GMstr=header(26:43); GMexp=header(46:47);
% 
% R=str2double(Rstr)*10^str2double(Rexp)*10^3; % In meters!
% GM=str2double(GMstr)*10^str2double(GMexp)*10^9; % m^3/s^2

if lmcosi(1,1)==1
    lmcosi=[zeros(1,size(lmcosi,2)); lmcosi];
end



end