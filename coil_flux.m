clear

tic

Msat = 1e6;

mu0 = 4*pi*1e-7;

px = linspace(-5e-3,5e-3,100);
py = linspace(-5e-3,5e-3,100);
pz = 1e-6;

area = ( px(2)-px(1))*( py(2)-py(1));

mag_size = [2e-5,2e-5,2e-8];

mask = ones(length(px), length(py));
rmask = 5e-3;

nP = 0; 


for pz = linspace(1e-6,1e-4,101)
    nP = nP+1;
for nX = 1:length(px)
    for nY = 1:length(py)
        for nZ = 1:length(pz)
            
            [Akoun(nP).HxAkoun(nX,nY,nZ), Akoun(nP).HyAkoun(nX,nY,nZ), ...
                Akoun(nP).HzAkoun(nX,nY,nZ)] = Jannsen(px(nX),py(nY),...
                pz(nZ),mag_size);
            
            if (px(nX)^2)+(py(nY)^2) >= rmask^2
                mask(nX,nY) = 0;
            end
            
        end 
    end 
end
if rem(nP,100) == 0
nP
end
[Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun] = multiply(...
    Msat*mu0/4/pi,Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun);
Akoun(nP).coilB = mask.*Akoun(nP).HzAkoun*area;
Akoun(nP).flux = sum(sum(Akoun(nP).coilB));
flux(nP) = Akoun(nP).flux;
pzplot(nP) = pz;
end 


toc

figure(1)
clf
semilogy(pzplot,flux)
title 'Flux captured by a coil with changing coil height d_z'
xlabel 'd_z [m]', ylabel 'Flux captured (\phi) [wb]'

% Akoun(nP).modBAkoun = sqrt(Akoun(nP).HxAkoun.^2 + Akoun(nP).HyAkoun.^2...
% + Akoun(nP).HzAkoun.^2);
