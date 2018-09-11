%% The old version of this code did not take into account the fact that this needed to be for a dipole. This takes this into account. 

% SI units unless stated. 

% JDZ 04/09/2018.

%clc 
clear

tic
% Start looking at the spatial parameters. 

% Think about how this might want to be done, could be changed for a GUI.

mu0 = 4* pi * 10^-7; % [H/m] SI.

grid_size = [50,50,50];

world_range = [10^-3, 10^-3, 10^-3];

extra = 0; %  10^-10; % [m]

space.Xline = linspace(-world_range(1),world_range(1),grid_size(1))+extra;
space.Yline = linspace(-world_range(2),world_range(2),grid_size(2))+extra;
space.Zline = linspace(-world_range(3),world_range(3),grid_size(3))+extra;

[space.X,space.Y,space.Z] = meshgrid(space.Xline, space.Yline, space.Zline);

space.radialN = sqrt(space.X.^2 + space.Y.^2 + space.Z.^2); %[m^3]

space.volume = prod(world_range./grid_size); % [m^3]


% Now making the magnetic vector.

Mag = zeros(grid_size(1),grid_size(2), grid_size(3));

unimagV = [ 0 0 1 ];

% Add in the dirac delta function somewhere

Mag(grid_size(1)/2:grid_size(1)/2+1,grid_size(2)/2:grid_size(2)/2+1,grid_size(3)/2)=1;
Mag(grid_size(1)/2:grid_size(1)/2+1,grid_size(2)/2:grid_size(2)/2+1,grid_size(3)/2+1)=-1;

Msat = 10^6; % [Am^2]

scaling = Msat*space.volume; % To be factored into the Dirac delta functions 

Mag = Mag.*scaling;

[space(1).find,space(2).find,space(3).find]=ind2sub(size(Mag), find(Mag));

space(1).totX = zeros(size(Mag));   space(1).totY = zeros(size(Mag));
space(1).totZ = zeros(size(Mag));   space(1).totB = zeros(size(Mag));

space(1).totXffeq = zeros(size(Mag));   space(1).totYffeq = zeros(size(Mag));
space(1).totZffeq = zeros(size(Mag));   space(1).totBffeq = zeros(size(Mag));


space(1).totXcrl = zeros(size(Mag));   space(1).totYcrl = zeros(size(Mag));
space(1).totZcrl = zeros(size(Mag));   space(1).totBcrl = zeros(size(Mag));   

for nP = 1:length(find(Mag))
    space(nP).par = [space(1).find(nP) space(2).find(nP) space(3).find(nP)];
    space(nP).mom = Mag(space(nP).par(1),space(nP).par(2),space(nP).par(3));
    space(nP).momV = space(nP).mom.*unimagV;
    space(nP).Rx = space(1).X + space(1).X(space(nP).par(1),space(nP).par(2),space(nP).par(3));
    space(nP).Ry = space(1).Y + space(1).Y(space(nP).par(1),space(nP).par(2),space(nP).par(3));
    space(nP).Rz = space(1).Z + space(1).Z(space(nP).par(1),space(nP).par(2),space(nP).par(3));
    space(nP).modR = sqrt(space(nP).Rx.^2 + space(nP).Ry.^2 + space(nP).Rz.^2);
    
    
    % finding the field
    
    for nX = 1:length(space(1).Xline)
        for nY = 1:length(space(1).Yline)
            for nZ = 1:length(space(1).Zline)
                
        rvec = [space(nP).Rx(nX,nY,nZ) space(nP).Ry(nX,nY,nZ) space(nP).Rz(nX,nY,nZ)];
        B = mu0/4/pi*(3*((dot(space(nP).momV,rvec).*rvec/space(nP).modR(nX,nY,nZ))-(space(nP).momV/space(nP).modR(nX,nY,nZ).^3)));
        [space(nP).Bx(nX,nY,nZ), space(nP).By(nX,nY,nZ), space(nP).Bz(nX,nY,nZ)] = field_comps(B(1), B(2), B(3));  
        
        Bffeq = mu0/4/pi*space(nP).mom.*rvec./(space(nP).modR(nX,nY,nZ)^3);
        [space(nP).Bxffeq(nX,nY,nZ), space(nP).Byffeq(nX,nY,nZ), space(nP).Bzffeq(nX,nY,nZ)] = field_comps(Bffeq(1), Bffeq(2), Bffeq(3)); 
        
        A = mu0/4/pi.*cross(space(nP).momV, rvec)./(space(nP).modR(nX,nY,nZ)^3);
        [space(nP).Ax(nX,nY,nZ), space(nP).Ay(nX,nY,nZ), space(nP).Az(nX,nY,nZ)] = field_comps(A(1), A(2), A(3));
        
            end 
        end
    end
    
    space(nP).modB = sqrt(space(nP).Bx.^2 + space(nP).By.^2 + space(nP).Bz.^2);
    space(1).totX = space(1).totX + space(nP).Bx;
    space(1).totY = space(1).totY + space(nP).By;
    space(1).totZ = space(1).totZ + space(nP).Bz;
    space(1).totB = space(1).totB + space(nP).modB;
    
    space(nP).modBffeq = sqrt(space(nP).Bxffeq.^2 + space(nP).Byffeq.^2 + space(nP).Bzffeq.^2);
    space(1).totXffeq = space(1).totXffeq + space(nP).Bxffeq;
    space(1).totYffeq = space(1).totYffeq + space(nP).Byffeq;
    space(1).totZffeq = space(1).totZffeq + space(nP).Bzffeq;
    space(1).totBffeq = space(1).totBffeq + space(nP).modBffeq;
    
    [space(nP).Ax(nX,nY,nZ), space(nP).Ay(nX,nY,nZ), space(nP).Az(nX,nY,nZ)] = field_comps(A(1), A(2), A(3));
    [space(nP).curlx,space(nP).curly,space(nP).curlz,space(nP).cav] = curl(space(1).X,space(1).Y,space(1).Z,space(nP).Ax, space(nP).Ay, space(nP).Az);
    space(1).totXcrl = space(1).totXcrl + space(nP).curlx;
    space(1).totYcrl = space(1).totYcrl + space(nP).curly;
    space(1).totZcrl = space(1).totZcrl + space(nP).curlz;
    space(nP).modBcrl = sqrt(space(nP).curlx.^2 + space(nP).curly.^2 + space(nP).curlz.^2); 
    space(1).totBcrl = space(1).totBcrl + space(nP).modBcrl;
    
end 

% Just to debug the system
debug.odd = space(1).Bz + space(2).Bz + space(3).Bz + space(4).Bz;
debug.even = space(5).Bz + space(6).Bz + space(7).Bz + space(8).Bz;
debug.btot = debug.even + debug.odd;

debug.oddffeq = space(1).Bzffeq + space(2).Bzffeq + space(3).Bzffeq + space(4).Bzffeq;
debug.evenffeq = space(5).Bzffeq + space(6).Bzffeq + space(7).Bzffeq + space(8).Bzffeq;
debug.bffeqtot = debug.evenffeq + debug.oddffeq;



filter = [10^-4, 10^-4, 10^-4];
[space(1).totXfil, space(1).totYfil, space(1).totZfil] = quiver_filter(filter, space(1).totX, space(1).totX, space(1).totX);

figure(1)
subplot(1,2,1)
slice(space(1).X,space(1).Y,space(1).Z, debug.odd, 0,0,0)
colorbar 
caxis([-0.001,0])
title 'Positive charges'
subplot(1,2,2)
slice(space(1).X,space(1).Y,space(1).Z, debug.even, 0,0,0)
colorbar 
caxis([0,0.001])
title 'Negative charges'

figure(2)
slice(space(1).X,space(1).Y,space(1).Z, space(1).totZ, 0,0,0)
caxis([-0.0001,0.0001])
colorbar
%hold on 
%quiver3(space(1).X,space(1).Y,space(1).Z,space(1).totXfil,space(1).totYfil,space(1).totZfil);


toc