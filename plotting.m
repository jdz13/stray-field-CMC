%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the Janssen Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partno = 1;

figure(3)
clf
subplot(2,2,1)
slice(space(1).X,space(1).Y,space(1).Z,mu0.*space(partno).HxAkoun , 0.0002,0,[])
caxis([-0.00000001,0.00000001])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(space(1).X,space(1).Y,space(1).Z,mu0.*space(partno).HyAkoun , 0,0.0002,[])
caxis([-0.00000001,0.00000001])
polarmap
colorbar
title 'Y component'

subplot(2,2,3)
slice(space(1).X,space(1).Y,space(1).Z,mu0.*space(partno).HzAkoun , 0,0,[])
caxis([-0.00000001,0.00000001])
polarmap
colorbar
title 'Z component'


subplot(2,2,4)
slice(space(1).X,space(1).Y,space(1).Z,sqrt(space(partno).HxAkoun.^2+space(partno).HyAkoun.^2+space(partno).HzAkoun.^2)  , 0,0,[])
caxis([0,0.001])
polarmap
colorbar
title 'tot B Janssen'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the far field equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
clf
subplot(2,2,1)
slice(space(1).X,space(1).Y,space(1).Z,space(partno).totXffeq , 0.0002,0,[])
caxis([-0.00000001,0.00000001])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(space(1).X,space(1).Y,space(1).Z,space(partno).totYffeq , 0,0.0002,[])
caxis([-0.00000001,0.00000001])
polarmap
colorbar
title 'Y component'

subplot(2,2,3)
slice(space(1).X,space(1).Y,space(1).Z,space(partno).totZffeq , 0,0,[])
caxis([-0.00000001,0.00000001])
polarmap
colorbar
title 'Z component'

subplot(2,2,4)
slice(space(1).X,space(1).Y,space(1).Z,space(1).totBffeq , 0,0,[])
caxis([0,0.0000001])
polarmap
colorbar
title 'tot B Ffeq'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the dipole equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
clf
subplot(2,2,1)
slice(space(1).X,space(1).Y,space(1).Z,space(partno).Bx , 0.0005,0,[])
caxis([-0.00001,0.00001])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(space(1).X,space(1).Y,space(1).Z,space(partno).By , 0,0.0005,[])
caxis([-0.00001,0.00001])
polarmap
colorbar
title 'Y component'

subplot(2,2,3)
slice(space(1).X,space(1).Y,space(1).Z,space(partno).Bz , 0,0,[])
caxis([-0.00001,0.00001])
polarmap
colorbar
title 'Z component'

subplot(2,2,4)
slice(space(1).X,space(1).Y,space(1).Z,space(1).modB , 0,0,[])
caxis([-0.0001,0.0001])
polarmap
colorbar
title 'Tot B Dipole eq'

%colormap(parula)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the vector field equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
clf
subplot(2,2,1)
slice(space(1).X,space(1).Y,space(1).Z,space(partno).curlx , 0.0002,0,[])
caxis([-0.00001,0.00001])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(space(1).X,space(1).Y,space(1).Z,space(partno).curly, 0,0.0002,[])
caxis([-0.00001,0.00001])
polarmap
colorbar
title 'Y component'


subplot(2,2,3)
slice(space(1).X,space(1).Y,space(1).Z,space(partno).curlz , 0,0,[])
caxis([-0.00001,0.00001])
polarmap
colorbar
title 'Z component'

subplot(2,2,4)
slice(space(1).X,space(1).Y,space(1).Z,space(1).modBcrl , 0,0,[])
caxis([0,0.0001])
polarmap
colorbar
title 'Tot B Curl'

%colormap(parula)

%subplot(2,2,3)
%hold on
%quiver3(space(1).X,space(1).Y,space(1).Z,space(1).curlx,space(1).curly,space(1).curlz)
%%

plane = size(space(nP).curly,2)/2+6;

dat = space(nP).modB (:,:,plane);

quivdatX = space(nP).Bx(:,:,plane);

reszemat = size(quivdatX);
reszemat(reszemat == 1) = [];
quivdatX = reshape (quivdatX, reszemat);


quivdatZ = space(nP).Bz(:,:,plane);

reszemat = size(quivdatZ);
reszemat(reszemat == 1) = [];
quivdatZ = reshape (quivdatZ, reszemat);

quivdatY = space(nP).By(:,:,plane);

reszemat = size(quivdatY);
reszemat(reszemat == 1) = [];
quivdatY = reshape (quivdatY, reszemat);

quivx = space(1).X(:,:,plane);
quivz = space(1).Z(:,:,plane);
quivy = space(1).Y(:,:,plane);

figure(7)

clf
colorbar
imagesc(space(1).Xline, space(1).Yline,dat)
hold on 
quiver(quivx,quivy,quivdatX,quivdatZ)
title ('2D plane of 3D output')

%%

plane = size(space(nP).curly,2)/2;

dat = space(1).totBcrl(:,:,plane);
reszemat = size(dat);
reszemat(reszemat == 1) = [];
dat = reshape (dat, reszemat);

quivdatX = space(nP).curlx(:,:,plane);

reszemat = size(quivdatX);
reszemat(reszemat == 1) = [];
quivdatX = reshape (quivdatX, reszemat);


quivdatZ = space(nP).curlz(:,:,plane);

reszemat = size(quivdatZ);
reszemat(reszemat == 1) = [];
quivdatZ = reshape (quivdatZ, reszemat);

quivdatY = space(nP).curly(:,:,plane);

reszemat = size(quivdatY);
reszemat(reszemat == 1) = [];
quivdatY = reshape (quivdatY, reszemat);

quivx = space(1).X(:,:,plane);
quivz = space(1).Z(:,:,plane);
quivy = space(1).Y(:,:,plane);

figure(7)

clf
colorbar
imagesc(space(1).Xline, space(1).Yline,dat)
hold on 
quiver(quivx,quivy,quivdatX,quivdatZ)
title ('2D plane of 3D output')