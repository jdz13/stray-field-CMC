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
world_range = cell_size'.*(size(Bz));
grid_size = size(Bz);
BMumax = sqrt(Bx.^2 + By.^2 + Bz.^2);

pm_world_length = cell_size'.*(size(Bz))./2;

displace = [0,0,0];
px = linspace(-pm_world_length(1)+cell_size(1)/2,pm_world_length(1)-cell_size(1)/2,grid_size(1))+displace(1);
py = linspace(-pm_world_length(2)+cell_size(2)/2,pm_world_length(2)-cell_size(2)/2,grid_size(2))+displace(2);
pZ = linspace(-pm_world_length(3)+cell_size(3)/2,pm_world_length(3)-cell_size(3)/2,grid_size(3))+displace(3);
[X,Y,Z] = meshgrid(px,py,pZ);
[max] = MiddleVariableLine(BMumax);

%%
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
slice(X,Y,Z,Bz,0,0,0)

title '3D display of field'; colorbar
xlabel 'X (m)'; ylabel 'Y (m)'; zlabel 'Z(m)'

%%

figure(141)

subplot(2,2,1)
plot(pZ,max)
 xlabel 'Distance from magnet centre (m)'; ylabel 'Field (T)'
 title 'Field - distance plot for my PM - equivalent square'
vline(0.01); vline(-0.01); vline(0.003); vline(-0.003);
subplot(2,2,3)
semilogy(pZ,max)
 xlabel 'Distance from magnet centre (m)'; ylabel 'Field (T)'
 title 'Log-scale'
vline(0.01); vline(-0.01); vline(0.003); vline(-0.003);


subplot(2,2,4)
imagesc(pZ,py,reshape(BMumax(51,:,:),[100,150]))
vline(0.01); vline(-0.01); vline(0.003); vline(-0.003);
hline(0.01); hline(-0.01); hline(0.003); hline(-0.003);
xlabel 'Z(m)'; ylabel 'Y(m)'; title 'YZ planar view'
axis equal

subplot(2,2,2)
imagesc(pZ,px,reshape(BMumax(:,51,:),[100,150]))
vline(0.01); vline(-0.01); vline(0.003); vline(-0.003);
hline(0.01); hline(-0.01); hline(0.003); hline(-0.003);
xlabel 'Z(m)'; ylabel 'X(m)'; title 'XZ planar view'
axis equal