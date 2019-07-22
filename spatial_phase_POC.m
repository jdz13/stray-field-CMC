clear 
tic

PZ = 0.05;
swfield = [2005,1703,1385, 1099, 798, 579, 427, 261, 150, 40]./1e4; %[linspace(0.4, 0.1,7),linspace(0.08,0.02,4)];
theta_end = 6*pi; theta_n =  1441;
theta = linspace(0,theta_end, theta_n); 
spread = 1e-3; % sampspac =linspace(2.5e-4,2e-3,8);
PM_cl = 1e-2;

inputPZ = PZ + PM_cl/2;


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
tic

for pull= theta

xprime = sample.X.*cos(pull) - sample.Z.*sin(pull);
yprime = sample.Y;
zprime = sample.X.*sin(pull) + sample.Z.*cos(pull);

[Akoun(np)] = Janssen_with_meshgrid(xprime,yprime,zprime,PM_mag_size);

% rotfield.fbx = Akoun.HxAkoun*cos(-pull) - Akoun.HzAkoun.*sin(-pull);
% rotfield.fby = Akoun.HyAkoun;
rotfield(np).fbz = Akoun(np).HxAkoun.*sin(-pull) + Akoun(np).HzAkoun.*cos(-pull);
% rotfield.fbtot = sqrt(rotfield.fbx.^2 + rotfield.fby.^2 + rotfield.fbz.^2);

%-------------------------------------------------------------------------
% can put in different particles here. 

for zzz = 1:length(swfield)
rotfield(np).bzmask_sf1(:,:,zzz) = (rotfield(np).fbz >= swfield(zzz)) - (rotfield(np).fbz <= -swfield(zzz));
rotfield(np).corrmask_sf1(:,:,zzz) = rotfield(np).bzmask_sf1(:,:,zzz) .* particle_loc;
vol_compa(np,zzz) = sum(sum(rotfield(np).corrmask_sf1(:,:,zzz)));
end

np = np + 1;
end

norm_vol_comp = vol_compa./control;
thetad = rad2deg(theta);
diphis = diff(vol_compa);

toc

%%

%points = [67,113; 157,163; 13,97; 29,129];

% points = [10,113; 30,113; 50,113; 70,113; 90,113; 110,113; 130,113; 150,113; 170,113; 190,113];
 points = [113,10; 113,30; 113,50; 113,70; 113,90; 113,110; 113,130; 113,150; 113,170; 113,190];


Bztheta = zeros(size(points,1), length(theta));

 for kk = 1:length(theta)
for jj = 1:size(points,1)
    
Bztheta(jj,kk) = rotfield(kk).fbz(points(jj,1),points(jj,2));

end
end


for kk = 1:size(points,1)
figure(20) ; hold on
plot(thetad, Bztheta(kk,:)*1000,'--')
xlabel 'Angle (degrees)'
ylabel 'Field value at that point (mT)'
end

%legend ('[67,113]','[157,163]','[13,97]','[29,129]', 'location','Southeast')
%legend ('[10,113]','[30,113]', '[50,113]', '[70,113]', '[90,113]', '[110,113]', '[130,113]', '[150,113]', '[170,113]', '[190,113]')
legend ('[113,10]','[113,30]', '[113,50]', '[113,70]', '[113,90]', '[113,110]', '[113,130]', '[113,150]', '[113,170]', '[113,190]')
title 'Field dependence on angle for different locations on the sample plane'


%%

figure(100); clf;

for kk = 1:size(points,1)

    ftans = fft(Bztheta(kk,:));
    
subplot(2,2,1); hold on
plot(abs(ftans))
xlabel 'arbitrary'
ylabel 'Amplitude (not normalised)'
title 'Amplitude of FFT results'

subplot(2,2,3); hold on
plot(abs(ftans))
xlabel 'arbitrary'
ylabel 'Amplitude (not normalised)'
title 'Zoomed section of above'

subplot(2,2,2); hold on
plot(rad2deg(angle(ftans)))
xlabel 'arbitrary'
ylabel 'phase (not normalised)'
title 'Phase of FFT results'

subplot(2,2,4); hold on
plot(rad2deg(angle(ftans)))
xlabel 'arbitrary'
ylabel 'phase (not normalised)'
title 'Zoomed section of above'
end
