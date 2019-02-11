clear

tic

Msat = 1e6;

mu0 = 4*pi*1e-7;

px = linspace(-5e-3,5e-3,500);
py = linspace(-5e-3,5e-3,500);
pz = 1e-6;

area = ( px(2)-px(1))*( py(2)-py(1));

mag_size = [2e-5,2e-5,2e-8];

mask = ones(length(px), length(py));
rmask = 5e-3;

nP = 0; 
mP = 0;

for pz = linspace(0,1e-4,100)
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
disp(nP)
end
[Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun] = multiply(...
    Msat*mu0/4/pi,Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun);
Akoun(nP).coilB = mask.*Akoun(nP).HzAkoun*area;
Akoun(nP).flux = sum(sum(Akoun(nP).coilB));
Akoun(nP).max = (Akoun(nP).coilB(50,50));
max(nP) = Akoun(nP).max;
flux(nP) = Akoun(nP).flux;
pzplot(nP) = pz;
end 


for rmask = linspace(5e-4,5e-3,21)
    mP = mP+1;
    Mask_fac(mP).mask = zeros(length(px), length(py));
    for nX = 1:length(px)
        for nY = 1:length(py)
            for nZ = 1:length(pz)     
                if (px(nX)^2)+(py(nY)^2) <= rmask^2
                    Mask_fac(mP).mask(nX,nY) = 1;
                end
            end
        end
    end
    rmaskplot(mP) = rmask;
end 

parasw = zeros(length(Akoun), length(Mask_fac));

for xP = 1:length(Akoun)
    for yP = 1:length(Mask_fac)
        Akoun(nP).coilB = mask.*Akoun(nP).HzAkoun*area;
        Akoun(nP).flux = sum(sum(Akoun(nP).coilB));
        flux(nP) = Akoun(nP).flux;

        parasw(xP,yP)= sum(sum(Mask_fac(yP).mask.*Akoun(xP).HzAkoun*area));
        
    end 
end



toc

figure(1)
clf
plot(pzplot,flux)
title 'Flux captured by a coil with changing coil height d_z'
xlabel 'd_z [m]', ylabel 'Flux captured (\phi) [wb]'


figure(2)
clf
for mP = 1:length(rmaskplot)
imagesc(Mask_fac(mP).mask)
pause(0.3)
end

figure(3)
imagesc(rmaskplot, pzplot, parasw)
xlabel 'coil radius [m]', ylabel 'distance from particle d_z [m]'
title 'flux captured by a coil (\phi) [wb]'
colorbar

figure(4)
n = 1;
for t = [1,6,13,21]
subplot(2,2,n)
imagesc(px,py,Mask_fac(t).mask)
n = n+1;
title (['mask for radii r_c_o_i_l =  ',num2str(rmaskplot(t)),'m'])
end

% Akoun(nP).modBAkoun = sqrt(Akoun(nP).HxAkoun.^2 + Akoun(nP).HyAkoun.^2...
% + Akoun(nP).HzAkoun.^2);
