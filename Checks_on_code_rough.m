 n% spy is also a great tool for 2D ploits, not so good for 3D

BtotMumax = sqrt(Bx.^2 + By.^2  + Bz.^2);



Test_struc.TESTVI = ((Akoun(1).modBAkoun)-BtotMumax)./BtotMumax;% ((Akoun(1).totB)-BtotMumax)./BtotMumax;
Test_struc.Testproof = Test_struc.TESTVI<= 0.001;

Test_struc.Testproof1pc = abs(Test_struc.TESTVI)>= 0.01;
Test_struc.s1pc = nonzeros(Test_struc.Testproof1pc);

Test_struc.Testproofp2pc = abs(Test_struc.TESTVI)>= 0.002;
Test_struc.sp2pc = nonzeros(Test_struc.Testproofp2pc);

Test_struc.Testproofp3pc = abs(Test_struc.TESTVI)>= 0.003;
Test_struc.sp3pc = nonzeros(Test_struc.Testproofp3pc);

Test_struc.Testproofp4pc = abs(Test_struc.TESTVI)>= 0.004;
Test_struc.sp4pc = nonzeros(Test_struc.Testproofp4pc);

Test_struc.Testproofp5pc = abs(Test_struc.TESTVI)>= 0.005;
Test_struc.sp5pc = nonzeros(Test_struc.Testproofp5pc);


ind = find(Test_struc.Testproof1pc);
[fx, fy, fz] = ind2sub(size(Test_struc.Testproof1pc), ind);

for trail = 1:length(fx)
values(trail) = Test_struc.TESTVI(fx(trail),fy(trail),fz(trail));
loc (trail) = {[num2str(fx(trail)),',',num2str(fy(trail)),',',num2str(fz(trail))]};
end
clear trail 

Test_struc.finders.fx = fx; Test_struc.finders.fy = fy; Test_struc.finders.fz = fz; 
Test_struc.finders.ind = ind; Test_struc.graph.values = values; Test_struc.graph.loc = loc;

clear ind fx fy fz values loc

Test_struc.graph.cats = categorical(Test_struc.graph.loc);
%%
figure(12)
subplot(1,3,1)
bar(Test_struc.graph.cats,Test_struc.graph.values*100); xlabel 'Location';  ylabel 'Percentage error'
title 'Percentage error locations over 1%' ;
subplot(1,3,[2,3])
slice(space(1).X,space(1).Y,space(1).Z, Test_struc.TESTVI*100, 0,0,[])
title 'Percentage error over all space'; colorbar; 
caxis([-1,1]);
