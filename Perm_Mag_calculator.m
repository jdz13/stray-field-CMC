clear

tic

Msat = 1.27e6; % M_sat of the particle - Neodymium has Msat 1.27e-6 [A/m];

mu0 = 4* pi * 10^-7; % [H/m] SI.

grid_size = [50,50,50];

world_range =  [89.5*10^-2, 89.5*10^-2, 89.5*10^-2];
cell_size = 2.*world_range./grid_size;

px = linspace(-world_range(1)+cell_size(1)/2,world_range(1)-cell_size(1)/2,grid_size(1));
py = linspace(-world_range(2)+cell_size(2)/2,world_range(2)-cell_size(2)/2,grid_size(2));
pZ = linspace(-world_range(3)+cell_size(3)/2,world_range(3)-cell_size(3)/2,grid_size(3)); % If only one plane is wanted, fill this here [m]

cl = 1.79e-2; % If cuboidal, this is the dimension [m]

mag_size = [cl,cl,cl];

n = 50; % How many planes to evaluate over

% Initialise the variables that need to be
nP = 0; mP = 0; max = zeros(1,n);

for pz = pZ % Loop over a set of planes in z [m]
    nP = nP+1; % Counting variable for indexing
for nX = 1:length(px)
    for nY = 1:length(py)
        for nZ = 1:length(pz)
            % Running the function for each variable
            [Akoun(nP).HxAkoun(nX,nY,nZ), Akoun(nP).HyAkoun(nX,nY,nZ), ...
                Akoun(nP).HzAkoun(nX,nY,nZ)] = Jannsen(px(nX),py(nY),...
                pz(nZ),mag_size);
            
        end 
    end 
end
if rem(nP,100) == 0
nP % checker for longer sets, outputs every hundred so know where it's at. 
end
[Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun] = multiply(...
    Msat*mu0/4/pi,Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun);

Akoun(nP).modBAkoun = sqrt(Akoun(nP).HxAkoun.^2 + Akoun(nP).HyAkoun.^2 + Akoun(nP).HzAkoun.^2);
pzplot(nP) = pz; % Multiplying in the constants.
max(nP) = MiddleVariableLine(Akoun(nP).modBAkoun); % Finds the middle value for each plane
end 
toc

%Plot this up. Show the field dependence at a distance. 
figure(1)
clf
semilogy(pzplot,max)
title (['Field - distance dependance for a (' num2str(cl*100) 'cm)^3' ...
   ' cuboidal Nd magnet'])
xlabel 'Distance from the defined centre of the magnet'
ylabel 'Field value [T]'
xL = xlim;
line(xL, [0 0],'k');  % zero on the y-axis