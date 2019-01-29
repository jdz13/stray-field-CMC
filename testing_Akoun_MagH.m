%clc 
clear

tic
% Start looking at the spatial parameters. 

% Think about how this might want to be done, could be changed for a GUI.

mu0 = 4* pi * 10^-7; % [H/m] SI.

grid_size = [50,50,50];

world_range = [10^-3, 10^-3, 10^-3];
cell_size = 2.*world_range./grid_size;

extra = 0; %  10^-10; % [m]

space.Xline = linspace(-world_range(1),world_range(1),grid_size(1))+extra;
space.Yline = linspace(-world_range(2),world_range(2),grid_size(2))+extra;
space.Zline = linspace(-world_range(3),world_range(3),grid_size(3))+extra;

[space.X,space.Y,space.Z] = meshgrid(space.Xline, space.Yline, space.Zline);

space.radialN = sqrt(space.X.^2 + space.Y.^2 + space.Z.^2); %[m^3]

space.volume = prod(2.*world_range./grid_size); % [m^3]


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
                                
        [Akoun(nP).HxAkoun(nX,nY,nZ), Akoun(nP).HyAkoun(nX,nY,nZ), Akoun(nP).HzAkoun(nX,nY,nZ)] = Jannsen(space(nP).Rx(nX,nY,nZ),space(nP).Ry(nX,nY,nZ),space(nP).Rz(nX,nY,nZ),cell_size);
        
        [MagH(nP).HxMagH(nX,nY,nZ), MagH(nP).HyMagH(nX,nY,nZ), MagH(nP).HzMagH(nX,nY,nZ)] = Magstat(space(nP).Rx(nX,nY,nZ),space(nP).Ry(nX,nY,nZ),space(nP).Rz(nX,nY,nZ),cell_size);
        
            end 
        end
    end
    
    [Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun] = multiply(Msat*mu0/4/pi,Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun);
    Akoun(nP).modBAkoun = sqrt(Akoun(nP).HxAkoun.^2 + Akoun(nP).HyAkoun.^2 + Akoun(nP).HzAkoun.^2);
    
    [MagH(nP).HxMagH, MagH(nP).HyMagH, MagH(nP).HzMagH] = multiply (mu0*Msat, MagH(nP).HxMagH, MagH(nP).HyMagH, MagH(nP).HzMagH);
    MagH(nP).modBMagH = sqrt(MagH(nP).HxMagH.^2 + MagH(nP).HyMagH.^2 + MagH(nP).HzMagH.^2);
    
end 

% Just to debug the system
debug.odd = Akoun(1).HzAkoun + Akoun(2).HzAkoun + Akoun(3).HzAkoun + Akoun(4).HzAkoun;
debug.even = Akoun(5).HzAkoun + Akoun(6).HzAkoun + Akoun(7).HzAkoun + Akoun(8).HzAkoun;
debug.btot = debug.even + debug.odd;

toc


figure(1)
subplot(1,2,1)
slice(space(1).X,space(1).Y,space(1).Z, debug.odd, 0,0,0)
colorbar 
caxis([-0.0001,0.0001])
title 'Positive charges'
subplot(1,2,2)
slice(space(1).X,space(1).Y,space(1).Z, debug.even, 0,0,0)
colorbar 
caxis([-0.0001,0.0001])
title 'Negative charges'
