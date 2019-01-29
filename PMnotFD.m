function [Akoun, pzplot,max] = PMnotFD(mag_size, px,py,pZ, Msat, displacement)

tic

px = px+displacement(1); py = py+displacement(2); pZ = pZ+displacement(3);

mu0 = 4* pi * 10^-7; % [H/m] SI.

% Initialise the variables that need to be
nP = 0; max = zeros(1,length(pZ));
pzplot = zeros(1,length(pZ));

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
%disp(nP) % checker for longer sets, outputs every hundred so know where it's at. 
end
[Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun] = multiply(...
    Msat*mu0/4/pi,Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun);

Akoun(nP).modBAkoun = sqrt(Akoun(nP).HxAkoun.^2 + Akoun(nP).HyAkoun.^2 + Akoun(nP).HzAkoun.^2);
pzplot(nP) = pz; % Multiplying in the constants.
max(nP) = MiddleVariableLine(Akoun(nP).modBAkoun); % Finds the middle value for each plane
end 
toc

end
