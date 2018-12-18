% spy is also a great tool for 2D ploits, not so good for 3D


TEST66 = ((Akoun(partno).modBAkoun)-BtotMumax)./BtotMumax;
Testproof = TEST66<= 0.001;

Testproof1pc = TEST66>= 0.01;
s1pc = nonzeros(Testproof1pc);

Testproofp2pc = TEST66>= 0.002;
sp2pc = nonzeros(Testproofp2pc);
Testproofp3pc = TEST66>= 0.003;
sp3pc = nonzeros(Testproofp3pc);

Testproofp4pc = TEST66>= 0.004;
sp4pc = nonzeros(Testproofp4pc);
Testproofp5pc = TEST66>= 0.005;
sp5pc = nonzeros(Testproofp5pc);


for trail = 1:length(findx)
values(trail) = TEST66(findx(trail),findy(trail),findz(trail));
loc (trail) = {[num2str(findx(trail)),',',num2str(findy(trail)),',',num2str(findz(trail))]};
end

cats = categorical(loc);

figure(12)
bar(cats,values*100)
xlabel 'Location'
ylabel 'Percentage error'
title 'Percentage error locations over 1%'

