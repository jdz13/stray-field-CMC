function [Bx, By, Bz,range] = Mumax_data_tool()
%MUMAX_DATA_TOOL Summary of this function goes here
%   Detailed explanation goes here

[file,path] = uigetfile('*.*'); % interactively gets the filename
filename = [path,file]; % put it all back together

S = load(filename); % loads data from filename.

real = S.OOMMFData; % extract the actual B data out

rep = size(real); % find out how big it is, to locate the data in struct

usize = rep(2:4); % this is the x,y,z size


% extract the data. You have to resize as it thinks its an n-n-n-1 matrix. 
Bx = real(1,:,:,:,1);
Bx = reshape(Bx, usize); 

By = real(1,:,:,:,2); 
By = reshape(By, usize);

Bz = real(1,:,:,:,3); 
Bz = reshape(Bz, usize);


range = S.GridSize'.*usize; % calculate the range - equal to n(x,y,z) by...
% gridsize

fclose('all'); % close the file, as not to piss matlab off. 

end

%linex = linspace(-range(1)/2, range(1)/2,usize(1));
%liney = linspace(-range(2)/2, range(2)/2,usize(2));
%linez = linspace(-range(3)/2, range(3)/2,usize(3));
%[X,Y,Z] = meshgrid(linex, liney,linez);
%figure(1)
%quiver3(X,Y,Z, Bx,By,Bz)
