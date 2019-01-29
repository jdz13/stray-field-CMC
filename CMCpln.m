function [space,Akoun] = CMCpln(grid_size, world_range, unimagV, Msat)

%The old version of this code did not take into account the fact that this needed to be for a dipole. This takes this into account. 

% SI units unless stated. 

% JDZ 04/09/2018.

tic
% Start looking at the spatial parameters. 

% Think about how this might want to be done, could be changed for a GUI.

%%
mu0 = 4* pi * 10^-7; % [H/m] SI.

cell_size = 2.*world_range./grid_size;

space.Xline = linspace(-world_range(1)+cell_size(1)/2,world_range(1)-cell_size(1)/2,grid_size(1));
space.Yline = linspace(-world_range(2)+cell_size(2)/2,world_range(2)-cell_size(2)/2,grid_size(2));
space.Zline = linspace(-world_range(3)+cell_size(3)/2,world_range(3)-cell_size(3)/2,grid_size(3));

[space.X,space.Y,space.Z] = meshgrid( space.Yline, space.Xline, space.Zline);

space.radialN = sqrt(space.X.^2 + space.Y.^2 + space.Z.^2); %[m^3]

space.volume = prod(2.*world_range./grid_size); % [m^3]

% Now making the magnetic vector.

Mag = zeros(grid_size(1),grid_size(2), grid_size(3));

% Add in the dirac delta function somewhere

Mag(grid_size(1)/2:grid_size(1)/2+1,grid_size(2)/2:grid_size(2)/2+1,grid_size(3)/2)=1;
Mag(grid_size(1)/2:grid_size(1)/2+1,grid_size(2)/2:grid_size(2)/2+1,grid_size(3)/2+1)=1;

scaling = Msat*space.volume; % To be factored into the Dirac delta functions 

Mag = Mag.*scaling; % provides the moment from each individual volume.

% find each magnetic particle, and denote it's indicies. 
[space(1).find,space(2).find,space(3).find]=ind2sub(size(Mag), find(Mag));


% lines for preallocation of variables, should speed up the code. 
xl = zeros(grid_size);
Akoun(2).HxAkoun = xl;
for xx = 1:length(find(Mag))
Akoun(xx).HxAkoun = xl; Akoun(xx).HyAkoun = xl; Akoun(xx).HzAkoun = xl; Akoun(xx).modBAkoun = xl; 
space(xx).Rx = xl; space(xx).Ry = xl; space(xx).Rz = xl; space(xx).modR = xl;
end

Akoun(1).totB = 0;

%%
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
                
        % rvec = [space(nP).Rx(nX,nY,nZ) space(nP).Ry(nX,nY,nZ) space(nP).Rz(nX,nY,nZ)];
         
        [Akoun(nP).HxAkoun(nX,nY,nZ), Akoun(nP).HyAkoun(nX,nY,nZ), Akoun(nP).HzAkoun(nX,nY,nZ)] = Jannsen(space(nP).Rx(nX,nY,nZ),space(nP).Ry(nX,nY,nZ),space(nP).Rz(nX,nY,nZ),cell_size);
        
            end 
        end
    end
    
    [Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun] = multiply(Msat*mu0/4/pi,Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun);
    Akoun(nP).modBAkoun = sqrt(Akoun(nP).HxAkoun.^2 + Akoun(nP).HyAkoun.^2 + Akoun(nP).HzAkoun.^2);
       
    Akoun(1).totB = Akoun(1).totB + Akoun(nP).modBAkoun;
    
    
end 




toc
end
