tic
data_before = BMumax;
theta = 23*(pi/180);


%for theta = linspace(0,2*pi,37)

% A is your arbitrary affine transformation matrix. You define your 3x3 
% rotation matrix as normal, and put it in the (1:3,1:3) of the 4x4 
% matrix.
A = ([1 0 0 0;
    0 cos(theta)  -sin(theta) 0;
      0 sin(theta) cos(theta) 0;
      0 0 0 1]);

% Define your current grid, where "data_before" is the input data. We 
% make it symmetric about 0 to define the origin.
xn = [1:size(data_before,1)]; 
xn = xn-mean(xn);
yn = [1:size(data_before,2)];
yn = yn-mean(yn);
zn = [1:size(data_before,3)];
zn = zn-mean(zn);

[Xn,Yn,Zn] = meshgrid(xn,yn,zn);

% Here you create your transformation data from the rotation matrix:
tform = affine3d(A);

% And define the new set of coordinates, Xn etc
[XN,YN,ZN] = transformPointsForward(tform,Xn,Yn,Zn);

% and now interpolate the data :)
V = data_before;

Vn = interp3(Xn,Yn,Zn,V,XN,YN,ZN,'nearest',0); % This is total field.
%display('Data is interpolated')


% 
Bxpri = interp3(Xn,Yn,Zn,Bx,XN,YN,ZN,'nearest',0);
Bypri = interp3(Xn,Yn,Zn,By,XN,YN,ZN,'nearest',0);
Bzpri = interp3(Xn,Yn,Zn,Bz,XN,YN,ZN,'nearest',0);

%toc
%%

difhu = V - Vn;
difBX = Bxpri - Bx;
difBY = Bypri - By;
difBZ = Bzpri - Bz;

figure(7)
subplot(3,3,1); imagesc(squeeze(Bx(:,50,:))); title 'Original - X'
subplot(3,3,2); imagesc(squeeze(Bxpri(:,50,:))); title 'Transformed'
subplot(3,3,3); imagesc(squeeze(difBX(:,50,:))); title 'Difference'

subplot(3,3,4); imagesc(squeeze(By(:,50,:))); title 'Original - Y'
subplot(3,3,5); imagesc(squeeze(Bypri(:,50,:))); title 'Transformed'
subplot(3,3,6); imagesc(squeeze(difBY(:,50,:))); title 'Difference'

subplot(3,3,7); imagesc(squeeze(Bz(:,50,:))); title 'Original - Z'
subplot(3,3,8); imagesc(squeeze(Bzpri(:,50,:))); title 'Transformed'
subplot(3,3,9); imagesc(squeeze(difBZ(:,50,:))); title 'Difference'

%%

backX = Bxpri.*cos(-theta) - Bzpri.*sin(-theta);
backY = Bypri;
backZ = Bxpri.*sin(-theta) + Bzpri.*cos(-theta);

figure(8)
subplot(1,3,1); imagesc(squeeze(backX(:,50,:))); title 'Actual X comp'
subplot(1,3,2); imagesc(squeeze(backY(:,50,:))); title 'Actual Y comp'
subplot(1,3,3); imagesc(squeeze(backZ(:,50,:))); title 'Actual Z comp'

%end 
toc
