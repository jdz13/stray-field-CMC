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
cell_size = world_range./grid_size;

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

field(1).totX = zeros(size(Mag));   field(1).totY = zeros(size(Mag));
field(1).totZ = zeros(size(Mag));   field(1).totB = zeros(size(Mag));

ffeq(1).totXffeq = zeros(size(Mag));   ffeq(1).totYffeq = zeros(size(Mag));
ffeq(1).totZffeq = zeros(size(Mag));   ffeq(1).totBffeq = zeros(size(Mag));

crl(1).totXcrl = zeros(size(Mag));   crl(1).totYcrl = zeros(size(Mag));
crl(1).totZcrl = zeros(size(Mag));   crl(1).totBcrl = zeros(size(Mag));   

for nP = 1:length(find(Mag))
    space(nP).par = [space(1).find(nP) space(2).find(nP) space(3).find(nP)];
    space(nP).mom = Mag(space(nP).par(1),space(nP).par(2),space(nP).par(3));
    space(nP).momV = space(nP).mom.*unimagV;
    space(nP).Rx = space(1).X + space(1).X(space(nP).par(1),space(nP).par(2),space(nP).par(3));
    space(nP).Ry = space(1).Y + space(1).Y(space(nP).par(1),space(nP).par(2),space(nP).par(3));
    space(nP).Rz = space(1).Z + space(1).Z(space(nP).par(1),space(nP).par(2),space(nP).par(3));
    space(nP).modR = sqrt(space(nP).Rx.^2 + space(nP).Ry.^2 + space(nP).Rz.^2);
    crl(nP).totXcrl = zeros(size(Mag));
    
    % finding the field
    
    for nX = 1:length(space(1).Xline)
        for nY = 1:length(space(1).Yline)
            for nZ = 1:length(space(1).Zline)
                
        rvec = [space(nP).Rx(nX,nY,nZ) space(nP).Ry(nX,nY,nZ) space(nP).Rz(nX,nY,nZ)];
        B = mu0/4/pi*((3.*dot(space(nP).momV,rvec).*rvec/(space(nP).modR(nX,nY,nZ).^5))-(space(nP).momV/space(nP).modR(nX,nY,nZ).^3));
        [field(nP).Bx(nX,nY,nZ), field(nP).By(nX,nY,nZ), field(nP).Bz(nX,nY,nZ)] = field_comps(B(1), B(2), B(3));  
        
        Bffeq = mu0/4/pi*space(nP).mom.*rvec./(space(nP).modR(nX,nY,nZ)^3);
        [ffeq(nP).Bxffeq(nX,nY,nZ), ffeq(nP).Byffeq(nX,nY,nZ), ffeq(nP).Bzffeq(nX,nY,nZ)] = field_comps(Bffeq(1), Bffeq(2), Bffeq(3)); 
        
        A = mu0/4/pi.*cross(space(nP).momV, rvec)./(space(nP).modR(nX,nY,nZ)^3);
        [crl(nP).Ax(nX,nY,nZ), crl(nP).Ay(nX,nY,nZ), crl(nP).Az(nX,nY,nZ)] = field_comps(A(1), A(2), A(3));
        
        [Akoun(nP).HxAkoun(nX,nY,nZ), Akoun(nP).HyAkoun(nX,nY,nZ), Akoun(nP).HzAkoun(nX,nY,nZ)] = Jannsen(space(nP).Rx(nX,nY,nZ),space(nP).Ry(nX,nY,nZ),space(nP).Rz(nX,nY,nZ),cell_size);
        
        [MagH(nP).HxMagH(nX,nY,nZ), MagH(nP).HyMagH(nX,nY,nZ), MagH(nP).HzMagH(nX,nY,nZ)] = Jannsen(space(nP).Rx(nX,nY,nZ),space(nP).Ry(nX,nY,nZ),space(nP).Rz(nX,nY,nZ),cell_size);
        
            end 
        end
    end
    
    field(nP).modB = sqrt(field(nP).Bx.^2 + field(nP).By.^2 + field(nP).Bz.^2);
    field(1).totX = field(1).totX + field(nP).Bx;
    field(1).totY = field(1).totY + field(nP).By;
    field(1).totZ = field(1).totZ + field(nP).Bz;
    field(1).totB = field(1).totB + field(nP).modB;
    
    ffeq(nP).modBffeq = sqrt(ffeq(nP).Bxffeq.^2 + ffeq(nP).Byffeq.^2 + ffeq(nP).Bzffeq.^2);
    ffeq(1).totXffeq = ffeq(1).totXffeq + ffeq(nP).Bxffeq;
    ffeq(1).totYffeq = ffeq(1).totYffeq + ffeq(nP).Byffeq;
    ffeq(1).totZffeq = ffeq(1).totZffeq + ffeq(nP).Bzffeq;
    ffeq(1).totBffeq = ffeq(1).totBffeq + ffeq(nP).modBffeq;
    
    [crl(nP).Ax(nX,nY,nZ), crl(nP).Ay(nX,nY,nZ), crl(nP).Az(nX,nY,nZ)] = field_comps(A(1), A(2), A(3));
    [crl(nP).curlx,crl(nP).curly,crl(nP).curlz,crl(nP).cav] = curl(space(1).X,space(1).Y,space(1).Z,crl(nP).Ax, crl(nP).Ay, crl(nP).Az);
    crl(1).totXcrl = crl(1).totXcrl + crl(nP).curlx;
    crl(1).totYcrl = crl(1).totYcrl + crl(nP).curly;
    crl(1).totZcrl = crl(1).totZcrl + crl(nP).curlz;
    crl(nP).modBcrl = sqrt(crl(nP).curlx.^2 + crl(nP).curly.^2 + crl(nP).curlz.^2); 
    crl(1).totBcrl = crl(1).totBcrl + crl(nP).modBcrl;
    
    
    [Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun] = multiply(Msat*mu0/4/pi,Akoun(nP).HxAkoun, Akoun(nP).HyAkoun, Akoun(nP).HzAkoun);
    Akoun(nP).modBAkoun = sqrt(Akoun(nP).HxAkoun.^2 + Akoun(nP).HyAkoun.^2 + Akoun(nP).HzAkoun.^2);
    
    [MagH(nP).HxMagH, MagH(nP).HyMagH, MagH(nP).HzMagH] = multiply (mu0*Msat/4/pi, MagH(nP).HxMagH, MagH(nP).HyMagH, MagH(nP).HzMagH);
    MagH(nP).modBMagH = sqrt(MagH(nP).HxMagH.^2 + MagH(nP).HyMagH.^2 + MagH(nP).HzMagH.^2);
    
end 

% Just to debug the system
debug.odd = field(1).Bz + field(2).Bz + field(3).Bz + field(4).Bz;
debug.even = field(5).Bz + field(6).Bz + field(7).Bz + field(8).Bz;
debug.btot = debug.even + debug.odd;

debug.oddffeq = ffeq(1).Bzffeq + ffeq(2).Bzffeq + ffeq(3).Bzffeq + ffeq(4).Bzffeq;
debug.evenffeq = ffeq(5).Bzffeq + ffeq(6).Bzffeq + ffeq(7).Bzffeq + ffeq(8).Bzffeq;
debug.bffeqtot = debug.evenffeq + debug.oddffeq;


%%
filter = [10^-4, 10^-4, 10^-4];
[field(1).totXfil, field(1).totYfil, field(1).totZfil] = quiver_filter(filter, field(1).totX, field(1).totX, field(1).totX);

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

figure(2)
slice(space(1).X,space(1).Y,space(1).Z, field(1).totZ, 0,0,0)
caxis([-0.0001,0.0001])
colorbar
%hold on 
%quiver3(space(1).X,space(1).Y,space(1).Z,field(1).totXfil,field(1).totYfil,field(1).totZfil);


toc