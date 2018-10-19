%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the MagH Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partno = 1;

figure(2)
clf
subplot(2,2,1)
slice(space(1).X,space(1).Y,space(1).Z,MagH(partno).HxMagH , 0.0002,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(space(1).X,space(1).Y,space(1).Z,MagH(partno).HyMagH , 0,0.0002,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Y component'

subplot(2,2,3)
slice(space(1).X,space(1).Y,space(1).Z,MagH(partno).HzMagH , 0,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Z component'


subplot(2,2,4)
slice(space(1).X,space(1).Y,space(1).Z,MagH(partno).modBMagH , 0,0,[])
caxis([0,0.0001])
polarmap
colorbar
title 'tot B MagH'

colormap (parula)

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the Janssen Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partno = 1;

checker = Akoun(1).modBAkoun./field(1).modB;
checker2 = ffeq(1).totBffeq./field(1).modB;
checker3 = MagH(1).modBMagH./field(1).modB;

figure(3)
clf
subplot(2,2,1)
slice(space(1).X,space(1).Y,space(1).Z,Akoun(partno).HxAkoun , 0.0002,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(space(1).X,space(1).Y,space(1).Z,Akoun(partno).HyAkoun , 0,0.0002,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Y component'

subplot(2,2,3)
slice(space(1).X,space(1).Y,space(1).Z,Akoun(partno).HzAkoun , 0,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Z component'


subplot(2,2,4)
slice(space(1).X,space(1).Y,space(1).Z,Akoun(partno).modBAkoun , 0,0,[])
caxis([0,0.0001])
polarmap
colorbar
title 'tot B Janssen'

colormap (parula)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the far field equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
clf
subplot(2,2,1)
slice(space(1).X,space(1).Y,space(1).Z,ffeq(partno).totXffeq , 0.0002,0,[])
caxis([-0.00000001,0.00000001])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(space(1).X,space(1).Y,space(1).Z,ffeq(partno).totYffeq , 0,0.0002,[])
caxis([-0.00000001,0.00000001])
polarmap
colorbar
title 'Y component'

subplot(2,2,3)
slice(space(1).X,space(1).Y,space(1).Z,ffeq(partno).totZffeq , 0,0,[])
caxis([-0.00000001,0.00000001])
polarmap
colorbar
title 'Z component'

subplot(2,2,4)
slice(space(1).X,space(1).Y,space(1).Z,ffeq(1).totBffeq , 0,0,[])
caxis([0,0.0000001])
polarmap
colorbar
title 'tot B Ffeq'
colormap (parula)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the dipole equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
clf
subplot(2,2,1)
slice(space(1).X,space(1).Y,space(1).Z,field(partno).Bx , 0.0002,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(space(1).X,space(1).Y,space(1).Z,field(partno).By , 0,0.0002,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Y component'

subplot(2,2,3)
slice(space(1).X,space(1).Y,space(1).Z,field(partno).Bz , 0,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Z component'

subplot(2,2,4)
slice(space(1).X,space(1).Y,space(1).Z,field(1).modB , 0,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Tot B Dipole eq'

colormap(parula)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the vector field equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
clf
subplot(2,2,1)
slice(space(1).X,space(1).Y,space(1).Z,crl(partno).curlx , 0.0002,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(space(1).X,space(1).Y,space(1).Z,crl(partno).curly, 0,0.0002,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Y component'


subplot(2,2,3)
slice(space(1).X,space(1).Y,space(1).Z,crl(partno).curlz , 0,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Z component'

subplot(2,2,4)
slice(space(1).X,space(1).Y,space(1).Z,crl(1).modBcrl , 0,0,[])
caxis([0,0.0001])
polarmap
colorbar
title 'Tot B Curl'

colormap(parula)

