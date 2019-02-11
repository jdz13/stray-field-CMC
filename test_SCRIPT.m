clear

Msat = 1e6; %1.27e6; % M_sat of the particle - Neodymium has Msat 1.27e-6 [A/m];

mu0 = 4* pi * 10^-7; % [H/m] SI.

displacement = [0,0,0]; %[-2e-5,-2e-5,-2e-5];

grid_size = [50,50,500];

world_range =  [5e-2,5e-2,1e-1]; % [1e-3,1e-3,1e-3];
cell_size = 2.*world_range./grid_size;

px = linspace(-world_range(1)+cell_size(1)/2,world_range(1)-cell_size(1)/2,grid_size(1));
py = linspace(-world_range(2)+cell_size(2)/2,world_range(2)-cell_size(2)/2,grid_size(2));
pZ = linspace(-world_range(3)+cell_size(3)/2,world_range(3)-cell_size(3)/2,grid_size(3)); % If only one plane is wanted, fill this here [m]

cl = 1.79e-2; % 4e-5; % If cuboidal, this is the dimension [m]



pZ = linspace(0,world_range(3), grid_size(3));

loe = 1;

Magcl = linspace(4e-3,2e-2,10);

for cl = Magcl

mag_size = [cl,cl,cl];
[Akoun, pzplot, field(loe).Bmax] = PMnotFD(mag_size, px,py,pZ, Msat, displacement);

figure(12); plot(pzplot,field(loe).Bmax); hold on

figure(13); semilogy(pzplot,field(loe).Bmax); hold on

loe = loe+1;
end
figure(12); xlabel 'Z distance (m)'; ylabel 'Max field (T)';
title 'Field - distance relationship for a cuboidal magnet of different sizes'
figure(13); xlabel 'Z distance (m)'; ylabel 'Max field (T)';
title 'Field - distance relationship for a cuboidal magnet of different sizes'

for count = 1:length(field)
halfMax = (min(field(count).Bmax) + max(field(count).Bmax)) / 2;
% Find where the data first drops below half the max.
index1 = find(field(count).Bmax >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(field(count).Bmax >= halfMax, 1, 'last');
fwhm = index2-index1 + 1; % FWHM in indexes.
% OR, if you have an x vector
fwhmx(count) = pzplot(index2) - pzplot(index1);
end
sanity = fwhmx - (Magcl./2);

figure(123)
plot(Magcl, sanity)
ylabel 'Decay length (m)'; xlabel 'Magnet cuboidal length(m)'

%%
%[Bx, By, Bz,range] = Mumax_data_tool();
%BtotMumax = sqrt(Bx.^2 + By.^2  + Bz.^2);
%e = 50;
%line1 = ['Field from a magnet of size (' num2str(mag_size*1e6) ')\mum'];
%line2 = [' displaced (' num2str(displacement*1e6) ')\mum with Msat ' num2str(Msat) 'A/m'];
%line3 = ['Percentage difference between PMnotFD and Mumax at a height (' num2str(pzplot(e)*1000) 'mm)'];
%limits = ([-world_range(1),world_range(1),-world_range(2),world_range(2)]);

%Test_struc.TESTVI = ((Akoun(e).modBAkoun)-BtotMumax(:,:,e))./BtotMumax(:,:,e);
%Test_struc.Testproof1pc = abs(Test_struc.TESTVI)>= 0.01;
%Test_struc.s1pc = nonzeros(Test_struc.Testproof1pc);
%figure(21); subplot(1,2,2); imagesc(px,py, 100.*Test_struc.TESTVI)
%title ({line3,[]}); axis(limits);xlabel 'x space (m)'; ylabel 'yspace (m)'; colorbar
%subplot(1,2,1); imagesc(px,py, Akoun(e).modBAkoun)
%title ({line1,line2,[]}); axis (limits);xlabel 'x space (m)'; ylabel 'yspace(m)'; colorbar