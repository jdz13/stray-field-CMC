tic

pltr1 = zeros(1,500);
pltr2 = pltr1;
pltr3 = pltr1;
pltr4 = pltr1;
pltr5 = pltr1;




swinit = 0.25; % What's the max channel value?
KRV = 5;

pltr1(1) = swinit;
sstart = 0.2500;

pulse = 1;

pm = 4;
% Find where that sits in space (Pz)
pzcut =  find(MxB(pm,:) <= sstart, 1, 'first')-1;
e = 1e-16; % tolerance - numerical rounding 



while abs(pltr1(pulse+1) - pltr1(pulse)) >= 50e-4
        
[tmp2,RK] = gamble_search(sstart, swinit,KRV,theta,Mdl_dtl, PZ, Bobj, particle_loc,control,MxB);

pltr1(pulse+1) = tmp2(1);
pltr2(pulse+1) = tmp2(2);
pltr3(pulse) = RK;
pltr4(pulse) = sstart;
pltr5(pulse) = swinit;


swinit = tmp2(1);
sstart = (pltr1(pulse+1) + pltr1(pulse))/2;

disp(['Count = ', num2str(pulse), ', range = ', num2str(tmp2), ', RK = ', num2str(RK)])

pulse = pulse+1;

end  

toc
