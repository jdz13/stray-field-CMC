clear
tic

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
sample.pZ = 3e-2+PMdisplacement(3); % If only one plane is wanted, fill this here [m]

[sample.X,sample.Y,sample.Z] = meshgrid(sample.px,sample.py,sample.pZ);

PM_cl = 1.79e-2;  % If cuboidal, this is the dimension [m]
PM_mag_size = [PM_cl,PM_cl,PM_cl]; % [m]

[PMfield(1), pzplot(1), Bmax] = PMnotFD(PM_mag_size, sample.px,sample.py,sample.pZ, MsatPM, PMdisplacement);

sample.r_spread = 1e-3;
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
      

theta_end = 4*pi; theta_n =  361;
theta = linspace(0,theta_end, theta_n); 
vol_compa = zeros(1,length(theta));


swfield = [0.4,0.3,0.2,0.1]; % [T]


Npart = 1000; % Number of particles, scalar factor

np = 1;
%%
for pull= theta

xprime = sample.X.*cos(pull) - sample.Z.*sin(pull);
yprime = sample.Y;
zprime = sample.X.*sin(pull) + sample.Z.*cos(pull);

[Akoun(np)] = Janssen_with_meshgrid(xprime,yprime,zprime,PM_mag_size);

rotfield(np).fbx = Akoun(np).HxAkoun*cos(-pull) - Akoun(np).HzAkoun.*sin(-pull);
rotfield(np).fby = Akoun(np).HyAkoun;
rotfield(np).fbz = Akoun(np).HxAkoun.*sin(-pull) + Akoun(np).HzAkoun.*cos(-pull);
rotfield(np).fbtot = sqrt(rotfield(np).fbx.^2 + rotfield(np).fby.^2 + rotfield(np).fbz.^2);

%-------------------------------------------------------------------------
% can put in different particles here. 

for zzz = 1:length(swfield)
rotfield(np).bzmask_sf1(:,:,zzz) = (rotfield(np).fbz >= swfield(zzz)) - (rotfield(np).fbz <= -swfield(zzz));
rotfield(np).corrmask_sf1(:,:,zzz) = rotfield(np).bzmask_sf1(:,:,zzz) .* particle_loc;
vol_compa(np, zzz) = sum(sum(rotfield(np).corrmask_sf1(:,:,zzz)));
end

np = np + 1;
end

norm_vol_comp = vol_compa./control;
thetad = rad2deg(theta);
diphis = diff(vol_compa);


figure(75); clf; hold on; figure(76); clf; hold on; 
for cx = 1:size(swfield,2)
figure(75)
plot(thetad, norm_vol_comp(:,cx))
figure(76)
plot(thetad(2:length(thetad)), diphis(:,cx))
end
figure(75); clear cx
legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
xlabel 'Angle (degrees)'; ylabel 'Normalised number of particles switched'
title 'Normalised number of particles switched for different angles' 
figure(76);
legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
xlabel 'Angle (degrees)'; ylabel 'Differentiated N particles switched'
title 'Differentiated number of switched for different angles'



% npartflux = flux_cap*Npart;
% angular_flux = norm_vol_comp.*npartflux;
% dphi = diff(angular_flux);
% rpm = 1200; 
% rotspeed = 2*pi * rpm/60; 
% dt = (theta(2)-theta(1))/rotspeed;
% EMF = -dphi./dt;

% t= 0:dt:((length(dphi))*dt);
% figure(1); plot(t(2:length(t)),EMF); xlabel 'Time [s]'; ylabel '\epsilon EMF (V)';
toc

