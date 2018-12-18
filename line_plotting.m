
midAkoun = MiddleVariableLine(Akoun(partno).modBAkoun);
midMagH = MiddleVariableLine(MagH(partno).modBMagH);
midB = MiddleVariableLine(field(1).modB);
midCrl = MiddleVariableLine(crl(1).modBcrl);

figure(17)
plot(space(1).Xline,midAkoun, space(1).Xline, midMagH,'ro',space(1).Xline, midB, 'bx', space(1).Xline, midCrl,'g--')
legend ('Akoun','MagH','Dipole eq', 'Curl')
ylabel 'Field B (T)', xlabel 'Distance from the particle (m)', title 'A comparison of codes'

figure(18)
semilogy(space(1).Xline,midAkoun, space(1).Xline, midMagH,'ro',space(1).Xline, midB, 'bx', space(1).Xline, midCrl,'g--')
legend ('Akoun','MagH','Dipole eq', 'Curl')
ylabel 'Field B (T)', xlabel 'Distance from the particle (m)', title 'A comparison of codes'

lineMumax = linspace(-0.001,0.001,50);
BtotMumax = sqrt(Bxmumax.^2+Bymumax.^2+ Bzmumax.^2);
midMumax = MiddleVariableLine(BtotMumax);


figure(20)
subplot(1,2,1)
semilogy(space(1).Xline,midAkoun, space(1).Xline, midMagH,'ro',space(1).Xline, midB, 'bx', space(1).Xline, midCrl,'g--',lineMumax,midMumax,'mo-')
legend ('Akoun','MagH','Dipole eq', 'Curl','Mumax')
ylabel 'Field B (T)', xlabel 'Distance from the particle (m)', title 'A comparison of codes'

midMumaxold = midMumax;
midMumax = midMumax;
 
%lower = linspace(1,length(midB),length(midB));

test11 = (midMumax-midAkoun)./midMumax; test22 = (midMumax-midMagH)./midMumax; test33 = (midMumax-midCrl)./midMumax; test44 = (midMumax-midB)./midMumax; test55 = (midMumax-midMumax)./midMumax;

subplot(1,2,2)
plot(space(1).Xline, test11,'bo',space(1).Xline,test22,'rx',space(1).Xline,test44,'r--', space(1).Xline, test55,'k--')
legend('Akoun', 'MagH', 'Dipole eq', 'Mumax')
title 'A comparison of codes, percentage difference'
xlabel 'Distance from the particle (m)'
ylabel '(Method - Mumax)/Mumax'

%figure(21)
%plot(space(1).Xline, test11,'bo',space(1).Xline,test22,'rx',space(1).Xline,test33,'k',space(1).Xline,test44,'r--', space(1).Xline, test55,'k--')
%legend('Akoun', 'MagH', 'Crl', 'Dipole eq', 'Mumax')

TEST66 = ((Akoun(partno).modBAkoun)-BtotMumax)./BtotMumax;
