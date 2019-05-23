function [norm_vol_comp, thetad, diphis] = rot_V2(inputPZ,swfield,theta,spread,PM_cl)

mu0 = 4* pi * 10^-7; % [H/m] SI.

MsatPM = 1.27e6; %1.27e6; % M_sat of the particle - Neodymium has Msat 1.27e-6 [A/m];
Msat_part = 1e6; % [A/m]

PMdisplacement = [0,0,0]; % Only 0,0,0 in this model. Becomes difficult in rotated space. 
sample.displacement = [0,0,0]; 

sample.plane_N = [200,200];
sample.plane_lengths =  [5e-3,5e-3]; % [m] plus/minus extents to 
sample.cell_size = 2.*sample.plane_lengths./sample.plane_N;

sample.px = linspace(-sample.plane_lengths(1)+sample.cell_size(1)/2,sample.plane_lengths(1)-sample.cell_size(1)/2,sample.plane_N(1))+PMdisplacement(1);
sample.py = linspace(-sample.plane_lengths(2)+sample.cell_size(2)/2,sample.plane_lengths(2)-sample.cell_size(2)/2,sample.plane_N(2))+PMdisplacement(2);
sample.pZ = inputPZ+PMdisplacement(3); % If only one plane is wanted, fill this here [m]

[sample.X,sample.Y,sample.Z] = meshgrid(sample.px,sample.py,sample.pZ);

PM_mag_size = [PM_cl,PM_cl,PM_cl]; % [m]

[PMfield(1), pzplot(1), Bmax] = PMnotFD(PM_mag_size, sample.px,sample.py,sample.pZ, MsatPM, PMdisplacement);

sample.r_spread = spread;
particle_loc = plane_mask(sample.px,sample.py,sample.r_spread); control = sum(sum(particle_loc));

sample_mag_size = [2e-5, 2e-5, 2e-8];

coil.plane_N = [200,200];
coil.plane_lengths =  [5e-2,5e-2]; % [m] plus/minus extents to 
coil.cell_size = 2.*sample.plane_lengths./sample.plane_N;

coil.px = linspace(-coil.plane_lengths(1)+coil.cell_size(1)/2,coil.plane_lengths(1)-coil.cell_size(1)/2,coil.plane_N(1));
coil.py = linspace(-coil.plane_lengths(2)+coil.cell_size(2)/2,coil.plane_lengths(2)-coil.cell_size(2)/2,coil.plane_N(2));
coil.pZ = 1e-4; % Only one plane is wanted, fill this here [m]

Rcoil = 5e-3;

[coilfield(1), coilpz(1), coilBmax] = PMnotFD(sample_mag_size, coil.px,coil.py,coil.pZ, Msat_part, sample.displacement);

coil_mask = plane_mask(coil.px,coil.py,Rcoil); area = (coil.px(2)-coil.px(1))*(coil.py(2)-coil.py(1));

flux_cap = sum(sum(coil_mask.*coilfield.HzAkoun*area));

vol_compa = zeros(length(theta),length(swfield));


Npart = 1000; % Number of particles, scalar factor

np = 1;
%%
for pull= theta

xprime = sample.X.*cos(pull) - sample.Z.*sin(pull);
yprime = sample.Y;
zprime = sample.X.*sin(pull) + sample.Z.*cos(pull);

[Akoun] = Janssen_with_meshgrid(xprime,yprime,zprime,PM_mag_size);

% rotfield.fbx = Akoun.HxAkoun*cos(-pull) - Akoun.HzAkoun.*sin(-pull);
% rotfield.fby = Akoun.HyAkoun;
rotfield.fbz = Akoun.HxAkoun.*sin(-pull) + Akoun.HzAkoun.*cos(-pull);
% rotfield.fbtot = sqrt(rotfield.fbx.^2 + rotfield.fby.^2 + rotfield.fbz.^2);

%-------------------------------------------------------------------------
% can put in different particles here. 

for zzz = 1:length(swfield)
rotfield.bzmask_sf1(:,:,zzz) = (rotfield.fbz >= swfield(zzz)) - (rotfield.fbz <= -swfield(zzz));
rotfield.corrmask_sf1(:,:,zzz) = rotfield.bzmask_sf1(:,:,zzz) .* particle_loc;
vol_compa(np,zzz) = sum(sum(rotfield.corrmask_sf1(:,:,zzz)));
end

np = np + 1;
end

norm_vol_comp = vol_compa./control;
thetad = rad2deg(theta);
diphis = diff(vol_compa);


end