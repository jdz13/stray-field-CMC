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

cell_size = S.GridSize;
world_range = cell_size.*(size(Bz));
grid_size = size(Bz);
BMumax = sqrt(Bx.^2 + By.^2 + Bz.^2);

displace = [0,0,0];
px = linspace(-world_range(1)+cell_size(1)/2,world_range(1)-cell_size(1)/2,grid_size(1))+displace(1);
py = linspace(-world_range(2)+cell_size(2)/2,world_range(2)-cell_size(2)/2,grid_size(2))+displace(2);
pZ = linspace(-world_range(3)+cell_size(3)/2,world_range(3)-cell_size(3)/2,grid_size(3))+displace(3);
[X,Y,Z] = meshgrid(px,py,pZ);
[max] = MiddleVariableLine(BMumax);


figure(140)

subplot(2,2,1)
plot(pZ,max)
 xlabel 'Distance from magnet centre (m)'; ylabel 'Field (T)'
 title 'Field - distance plot for my PM - equivalent square'
subplot(2,2,3)
semilogy(pZ,max)
 xlabel 'Distance from magnet centre (m)'; ylabel 'Field (T)'
 title 'Log-scale'
subplot(2,2,[2 4])
slice(X,Y,Z,BMumax,0,0,0)

title '3D display of field'; colorbar
xlabel 'X (m)'; ylabel 'Y (m)'; zlabel 'Z(m)'


