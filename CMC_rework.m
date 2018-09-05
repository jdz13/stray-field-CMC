%% The old version of this code did not take into account the fact that this needed to be for a dipole. This takes this into account. 

%clc 
clear

tic
% Start looking at the spatial parameters. 

% Think about how this might want to be done, could be changed for a GUI.

grid_size = [4,4,4];

world_range = [10^-3, 10^-3, 10^-3];

extra = 10^-10;

space.Xline = linspace(-world_range(1),world_range(1),grid_size(1))+extra;
space.Yline = linspace(-world_range(2),world_range(2),grid_size(2))+extra;
space.Zline = linspace(-world_range(3),world_range(3),grid_size(3))+extra;

[space.X,space.Y,space.Z] = meshgrid(space.Xline, space.Yline, space.Zline);

space.radialN = sqrt(space.X.^2 + space.Y.^2 + space.Z.^2);

space.volume = prod(world_range./grid_size); % [m^3]


% Now making the magnetic vector.

Mag = zeros(2*grid_size(1),2*grid_size(2), 2*grid_size(3));

% Add in the dirac delta function somewhere

Mag(grid_size(1):grid_size(1)+1,grid_size(2):grid_size(2)+1,grid_size(3))=1;
Mag(grid_size(1):grid_size(1)+1,grid_size(2):grid_size(2)+1,grid_size(3)+1)=-1;

Msat = 10^6; % [Am^2]

scaling = Msat*space.volume; % To be factored into the Dirac delta functions 

Mag = Mag.*scaling;

[space(1).find,space(2).find,space(3).find]=ind2sub(size(Mag), find(Mag));

toc