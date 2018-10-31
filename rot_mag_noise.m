% plots defining the noise from a misaligned permanent magnet, and the feild from a magnet of finite volume.

% second needs some thought as this is wrong. On how to define the magnet.  

clear

mu0 = 4*pi*10^-7; %[H/m]
Ms = 10^6; % [A/m]

ptdif = (1/1.01);
freq = linspace(1,200,200);
crad = linspace(0.5,5,200).*10^-3;
carea = pi.* (crad.^2);
[F,A] = meshgrid(freq,carea);
delta = 1 - ((ptdif)^3);
emf1 = 0.2.*F.*A.*delta;
imagesc(freq,crad,emf1)
title 'EMF_n_o_i_s_e (\epsilon [V]) for B_m_a_x=0.2T and \deltal=1%'
xlabel 'Frequency (Hz)', ylabel 'Coil radius [m]' 
colorbar


radm = linspace(2.5*10^-4,25*10^-3,100); %[m]
vol = 4/3*pi .*radm.^3; % [m^3]
mu = Ms .*vol; % [Am^2]
R = linspace(10^-3,100*10^-3,100); %[m]
[rmag, rdis] = meshgrid (mu,R);
B = mu0.*rmag./2./pi./(rdis.^3); %[T]
figure(2)
imagesc(radm,R,B)
xlabel 'magnet radius [m]',ylabel 'Distance (magnet to coil) [m]'
title ' B-distance relationship for different magnet radius'' - circle'
colorbar
caxis([0,0.2])