
tic

Msat = 1e6;

mu0 = 4*pi*1e-7;

px = linspace(-5e-3,5e-3,1000);
py = linspace(-5e-3,5e-3,1000);
pz = 1e-6;

[Px,Py,Pz] = meshgrid(px,py,pz);

mag_size = [2e-5,2e-5,2e-8];

mask = zeros(size(Px));
rmask = 5e-3;

nP = 1;

for nX = 1:length(px)
    for nY = 1:length(py)
        for nZ = 1:length(pz)
            
            [Akoun(nP).HxAkoun(nX,nY,nZ), Akoun(nP).HyAkoun(nX,nY,nZ), Akoun(nP).HzAkoun(nX,nY,nZ)] = Jannsen(px(nX),py(nY),pz(nZ),mag_size);
            
            if (px(nX)^2)+(py(nY)^2) <= rmask^2
                mask(nX,nY) = 1;
            end
            
        end 
    end 
end

[Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun] = multiply(Msat*mu0/4/pi,Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun);

Akoun(nP).modBAkoun = sqrt(Akoun(nP).HxAkoun.^2 + Akoun(nP).HyAkoun.^2 + Akoun(nP).HzAkoun.^2);

area = ( px(2)-px(1))*( py(2)-py(1));
coilB = mask.*Akoun(1).HzAkoun*area;

toc