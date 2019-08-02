

option1 = 0;
option2 = 0;
option3 = 1;

secondary_option = 2;


% first thing - do you need to import the data?
% 1 = yes
% 0 = no
if option1 == 1

    % first select X data, then Y, then Z;
    [Bobj.BXobject,Bobj.BXx,Bobj.BXy,Bobj.BXz] = Mumax_data_extractor();
    [Bobj.BYobject,Bobj.BYx,Bobj.BYy,Bobj.BYz] = Mumax_data_extractor();
    [Bobj.BZobject,Bobj.BZx,Bobj.BZy,Bobj.BZz] = Mumax_data_extractor();

elseif option1 == 0
        
else
        disp 'Please input either 1(yes) or 0(no) for option1'
end



% second thing - having a look at the original data 
% 1 for replot
% 0 to turn off

if option2 == 1
% lets have an initial look at the mumax results - need all three
szr = size(Bobj.BXx); rndhlfszr = round(szr./2);
fno = 81; n = 10;
figure(fno)
subplot(3,3,1); imagesc(squeeze(Bobj.BXx(:,:,n))); title 'X - x comp'
subplot(3,3,2); imagesc(squeeze(Bobj.BXy(:,:,n))); title 'X - y comp'
subplot(3,3,3); imagesc(squeeze(Bobj.BXz(:,:,n))); title 'X - z comp'

subplot(3,3,4); imagesc(squeeze(Bobj.BYx(:,:,n))); title 'Y - x comp'
subplot(3,3,5); imagesc(squeeze(Bobj.BYy(:,:,n))); title 'Y- y comp'
subplot(3,3,6); imagesc(squeeze(Bobj.BYz(:,:,n))); title 'Y - z comp'

subplot(3,3,7); imagesc(squeeze(Bobj.BZx(:,:,n))); title 'Z - x comp'
subplot(3,3,8); imagesc(squeeze(Bobj.BZy(:,:,n))); title 'Z - y comp'
subplot(3,3,9); imagesc(squeeze(Bobj.BZz(:,:,n))); title 'Z - z comp'

elseif option2 == 0
        
else
        disp 'Please input either 1(yes) or 0(no) for option1'
end


%% add in a theta component 

%fno = fno+1;
n = 20; % veriable to choose what plane to look at in the plots
%for theta = linspace(0,2*pi,37)

theta = 67*pi/180;

Bxnew = Bobj.BXx.*sin(theta) + Bobj.BZx.*cos(theta);
Bynew = Bobj.BXy.*sin(theta) + Bobj.BZy.*cos(theta);
Bznew = Bobj.BXz.*sin(theta) + Bobj.BZz.*cos(theta);


figure(fno); clf;
subplot(1,3,1); imagesc(squeeze(Bxnew(:,n,:))'); title 'New X component'; colorbar
subplot(1,3,2); imagesc(squeeze(Bynew(:,n,:))'); title 'New Y component'; colorbar
subplot(1,3,3); imagesc(squeeze(Bznew(:,n,:))'); title 'New Z component'; colorbar
pause (0.5)

%end 

%%

if option3 == 1

    if secondary_option == 1
 [B67.B67object,B67.Bx,B67.By,B67.Bz] = Mumax_data_extractor();
 [B23.B23object,B23.Bx,B23.By,B23.Bz] = Mumax_data_extractor();

    elseif secondary_option == 2
    
        checkBz67 = (Bznew-B67.Bz)./(Bznew);


        Test_struc.TESTVI = checkBz67;
        
        Test_struc.thresh1 = 0.01;
        Test_struc.Testproof1pc = abs(Test_struc.TESTVI)>= Test_struc.thresh1;
        Test_struc.s1pc = nonzeros(Test_struc.Testproof1pc);
        
        Test_struc.thresh2 = 1e-6;
        additional_test = Test_struc.Testproof1pc .* (abs(Bznew)> Test_struc.thresh2);
         
        ind = find(additional_test);
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

        figure(12)
        bar(Test_struc.graph.cats,Test_struc.graph.values*100); xlabel 'Location';  ylabel 'Percentage error'
        title ([ 'Percentage error locations over ',num2str(Test_struc.thresh1*100),'% with ', num2str(Test_struc.thresh2*1e6), '\muT threshold']) ; 
    elseif secondary_option == 0
        
    else 
        disp 'Please input either 1(yes) or 0(no) for secondary_option'

    end 
elseif option3 == 0
        
    else 
        disp 'Please input either 1(yes) or 0(no) for option3'
       
end 

%%

% Start looking at the probe plane. This will then allow us to do the
% logical matrix manipulation.

Mdl_dtl.CellNo = size(Bobj.BXx);
Mdl_dtl.gridsize = Bobj.BXobject.GridSize';
Mdl_dtl.wrldSz = Mdl_dtl.gridsize.* Mdl_dtl.CellNo;

Mdl_dtl.extents = Mdl_dtl.wrldSz./2;

Mdl_dtl.purelinex = linspace(-Mdl_dtl.extents(1)+Mdl_dtl.gridsize(1)/2,Mdl_dtl.extents(1)-Mdl_dtl.gridsize(1)/2,Mdl_dtl.CellNo(1));
Mdl_dtl.pureliney = linspace(-Mdl_dtl.extents(2)+Mdl_dtl.gridsize(2)/2,Mdl_dtl.extents(2)-Mdl_dtl.gridsize(2)/2,Mdl_dtl.CellNo(2));
Mdl_dtl.purelinez = linspace(-Mdl_dtl.extents(3)+Mdl_dtl.gridsize(3)/2,Mdl_dtl.extents(3)-Mdl_dtl.gridsize(3)/2,Mdl_dtl.CellNo(3));
[Mdl_dtl.Xp,Mdl_dtl.Yp,Mdl_dtl.Zp] = meshgrid(Mdl_dtl.purelinex,Mdl_dtl.pureliney,Mdl_dtl.purelinez);

Mdl_dtl.ID = 6e-3;
Mdl_dtl.OD = 2e-2;

Mdl_dtl.displace = [0,0,0.5*(-Mdl_dtl.OD + Mdl_dtl.wrldSz(3))];

Mdl_dtl.cntrmagLinex = Mdl_dtl.purelinex + Mdl_dtl.displace(1);
Mdl_dtl.cntrmagLiney = Mdl_dtl.pureliney + Mdl_dtl.displace(2);
Mdl_dtl.cntrmagLinez = Mdl_dtl.purelinez + Mdl_dtl.displace(3);
[Mdl_dtl.Xcm,Mdl_dtl.Ycm,Mdl_dtl.Zcm] = meshgrid(Mdl_dtl.cntrmagLinex,Mdl_dtl.cntrmagLiney,Mdl_dtl.cntrmagLinez);

Mdl_dtl.topmagLinez = Mdl_dtl.cntrmagLinez - Mdl_dtl.OD/2 + Mdl_dtl.gridsize(3)/2;

maxline = zeros(1,size(Bznew,3));
for Zin = 1:size(Bznew,3)
    tempvar = squeeze(Bznew(:,:,Zin));
    maxline(Zin) = max(abs(tempvar(:)));
end
fno = fno + 1; figure(fno);
plot((Mdl_dtl.cntrmagLinez), maxline*1e4)
xlabel 'Z distance from magnet centre'; ylabel 'Field maximum (Oe)';
title 'Maximum field profile for a Magnet'

%%

% Start defining and manipulating the probe field 

% Begin with the definition

Pzset = 50e-3;

e = 1e-16;
idx = find(abs(Mdl_dtl.topmagLinez - Pzset)<=e);

lkatpln = Bznew(:,:,idx);

% Now look at applying the logical matricies. 

sample.r_spread = 1e-3;
particle_loc = plane_mask(Mdl_dtl.cntrmagLinex,Mdl_dtl.cntrmagLiney,sample.r_spread); control = sum(sum(particle_loc));

swfield = 0.001;

vol_compa = zeros(length(theta),length(swfield));

np = 1;

for pull = theta
    
    
for zzz = 1:length(swfield)
rotfield.bzmask_sf1(:,:,zzz) = (lkatpln >= swfield(zzz)) - (lkatpln <= -swfield(zzz));
rotfield.corrmask_sf1(:,:,zzz) = rotfield.bzmask_sf1(:,:,zzz) .* particle_loc;
vol_compa(np,zzz) = sum(sum(rotfield.corrmask_sf1(:,:,zzz)));
end

np = np + 1;
end

norm_vol_comp = vol_compa./control;
thetad = rad2deg(theta);
diphis = diff(vol_compa);

