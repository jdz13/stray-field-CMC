clear

grid_size = [50,50,50];

world_range = [1e-3,1e-3,1e-3];

unimagV = [ 0 0 1 ];

Msat = 1e6; % [Am^2]

displacement = [0,0,0];

[space, Akoun] = CMC(grid_size, world_range, unimagV, Msat,displacement);

figure(2)
slice(space(1).X,space(1).Y,space(1).Z, Akoun(1).modBAkoun, 0,0,0)
caxis([-0.0001,0.0001])
colorbar




[Bx, By, Bz,range] = Mumax_data_tool();
BtotMumax = sqrt(Bx.^2 + By.^2  + Bz.^2);
%%
plane = 50;
Test_struc.TESTVI = ((Akoun(plane).modBAkoun)-BtotMumax(:,:,plane))./BtotMumax(:,:,plane);
Test_struc.Testproof1pc = abs(Test_struc.TESTVI)>= 0.01;
Test_struc.s1pc = nonzeros(Test_struc.Testproof1pc);

imagesc(px,py,Test_struc.TESTVI*100)
title (['Percentage errot comparison of PM code with Mumax for a plane height of ', num2str(pzplot(plane)*1e3),'mm'])
colorbar
