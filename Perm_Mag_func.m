function [Akoun, pzplot, max] = Perm_Mag_func(Msat,XYmax,Nxy,cl,...
    Nzplanes,Zmin,Zmax)

tic

mu0 = 4*pi*1e-7; % [H/m]

px = linspace(-XYmax,XYmax,Nxy); % Distance to evaluate field over in x [m]
py = linspace(-XYmax,XYmax,Nxy); % Distance to evaluate field over in y [m] 

mag_size = [cl,cl,cl];

% Initialise the variables that need to be
nP = 0; max = zeros(1,Nzplanes);

for pz = linspace(Zmin,Zmax,Nzplanes) % Loop over a set of planes in z [m]
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
pzplot(nP) = pz; % Multiplying in the constants.
max(nP) = (Akoun(nP).HzAkoun(round(Nxy/2),round(Nxy/2))); % Finds the max value for each plane
end 
toc

%Plot this up. Show the field dependence at a distance. 
figure(10)
plot(pzplot,max)
title (['Field - distance dependance for a (' num2str(cl*100) 'cm)^3' ...
   ' cuboidal Nd magnet'])
xlabel 'Distance from the defined centre of the magnet'
ylabel 'Field value [T]'
line(xlim, [0 0],'color','k');  % zero on the y-axis
hold on
end
