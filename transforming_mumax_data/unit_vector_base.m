
option1 = 0;

% first thing - do you need to import the data?
% 1 = yes
% 0 = no
if option1 == 1
<<<<<<< HEAD
    
    clear
    
=======

>>>>>>> f8b824162267bde2f5d97a0a169a1ec31eab8856
    % first select X data, then Y, then Z;
    [Bobj.BXobject,Bobj.BXx,Bobj.BXy,Bobj.BXz] = Mumax_data_extractor();
    [Bobj.BYobject,Bobj.BYx,Bobj.BYy,Bobj.BYz] = Mumax_data_extractor();
    [Bobj.BZobject,Bobj.BZx,Bobj.BZy,Bobj.BZz] = Mumax_data_extractor();

elseif option1 == 0
        
else
        disp 'Please input either 1(yes) or 0(no) for option1'
end

%%

% Bringing in the world definitions 

Mdl_dtl.CellNo = size(Bobj.BXx);
Mdl_dtl.gridsize = Bobj.BXobject.GridSize';
Mdl_dtl.wrldSz = Mdl_dtl.gridsize.* Mdl_dtl.CellNo;

Mdl_dtl.extents = Mdl_dtl.wrldSz./2;

Mdl_dtl.purelinex = linspace(-Mdl_dtl.extents(1)+Mdl_dtl.gridsize(1)/2,Mdl_dtl.extents(1)-Mdl_dtl.gridsize(1)/2,Mdl_dtl.CellNo(1));
Mdl_dtl.pureliney = linspace(-Mdl_dtl.extents(2)+Mdl_dtl.gridsize(2)/2,Mdl_dtl.extents(2)-Mdl_dtl.gridsize(2)/2,Mdl_dtl.CellNo(2));
Mdl_dtl.purelinez = linspace(-Mdl_dtl.extents(3)+Mdl_dtl.gridsize(3)/2,Mdl_dtl.extents(3)-Mdl_dtl.gridsize(3)/2,Mdl_dtl.CellNo(3));


Mdl_dtl.ID = 6e-3;
Mdl_dtl.OD = 2e-2;

Mdl_dtl.displace = [0,0,0.5*(-Mdl_dtl.OD + Mdl_dtl.wrldSz(3))];

Mdl_dtl.cntrmagLinex = Mdl_dtl.purelinex + Mdl_dtl.displace(1);
Mdl_dtl.cntrmagLiney = Mdl_dtl.pureliney + Mdl_dtl.displace(2);
Mdl_dtl.cntrmagLinez = Mdl_dtl.purelinez + Mdl_dtl.displace(3);

Mdl_dtl.topmagLinez = Mdl_dtl.cntrmagLinez - Mdl_dtl.OD/2 + Mdl_dtl.gridsize(3)/2;

%%
tic

% Start defining and manipulating the probe field 

% Begin with the definition

<<<<<<< HEAD

% Now look at applying the logical matricies. 

pm_cl = 2e-2;
sampspac = 1e-3;
particle_loc = plane_mask(Mdl_dtl.cntrmagLinex,Mdl_dtl.cntrmagLiney,sampspac); 
control = sum(sum(particle_loc));

swfield = [2005,1703,1385, 1099, 798,700, 579, 427, 261, 150, 40]./1e4; 
theta = linspace(0,2*pi,361);
PZ = linspace(1e-3,100e-3,100);

thetad = rad2deg(theta);


% pre allocate the structures 
 for ms = size(pm_cl,2):-1:1
    for sc = size(sampspac,2):-1:1
         for k = size(PZ,2):-1:1    % Backwards!
            norm_vol_comp(k,sc,ms).data = zeros(length(theta), size(swfield,2));
            diphis(k,sc,ms).data = zeros(length(theta), size(swfield,2));
         end
    end
 end 
clear k ms sc

for pm = 1:length(pm_cl)
for sp = 1:length(sampspac)
for nm = 1:length(PZ)
    
    e = 1e-16;
    idx = find(abs(Mdl_dtl.topmagLinez - PZ(nm))<=e);
    
for pull = 1:length(theta)
   %Bxnew = Bobj.BXx.*sin(theta) + Bobj.BZx.*cos(theta);
   %Bynew = Bobj.BXy.*sin(theta) + Bobj.BZy.*cos(theta);
   Bznew = Bobj.BXz.*sin(theta(pull)) + Bobj.BZz.*cos(theta(pull)); 
=======
Pzset = 20e-3;

e = 1e-16;
idx = find(abs(Mdl_dtl.topmagLinez - Pzset)<=e);

% Now look at applying the logical matricies. 

sample.r_spread = 1e-3;
particle_loc = plane_mask(Mdl_dtl.cntrmagLinex,Mdl_dtl.cntrmagLiney,sample.r_spread); control = sum(sum(particle_loc));

swfield = 0.0547; 
theta = linspace(0,2*pi,361);

vol_compa = zeros(length(theta),length(swfield));

np = 1;



for pull = theta
   %Bxnew = Bobj.BXx.*sin(theta) + Bobj.BZx.*cos(theta);
   %Bynew = Bobj.BXy.*sin(theta) + Bobj.BZy.*cos(theta);
   Bznew = Bobj.BXz.*sin(pull) + Bobj.BZz.*cos(pull); 
>>>>>>> f8b824162267bde2f5d97a0a169a1ec31eab8856
    
   lkatpln = Bznew(:,:,idx);
        for zzz = 1:length(swfield)
            rotfield.bzmask_sf1(:,:,zzz) = (lkatpln >= swfield(zzz)) - (lkatpln <= -swfield(zzz));
            rotfield.corrmask_sf1(:,:,zzz) = rotfield.bzmask_sf1(:,:,zzz) .* particle_loc;
<<<<<<< HEAD
            vol_compa = sum(sum(rotfield.corrmask_sf1(:,:,zzz)));
            norm_vol_comp(nm,sp,pm).data(pull,zzz) = vol_compa./control;
        end

end
diphis(nm,sp,pm).data = diff(norm_vol_comp(nm).data);
end
end 
end

toc

%%

option2 = 1;

chosen_plane = 56;

if option2 == 1
figure(3)
plot(thetad(2:length(thetad)),squeeze(diphis(chosen_plane).data()))
xlabel 'Angle (degrees)'; ylabel 'Differentiated signal (\deltaN)/(\delta\theta)'
end
=======
            vol_compa(np,zzz) = sum(sum(rotfield.corrmask_sf1(:,:,zzz)));
        end

    np = np + 1;
end

norm_vol_comp = vol_compa./control;
thetad = rad2deg(theta);
diphis = diff(vol_compa);
 
figure(2)
plot(thetad(2:length(thetad)),diphis)
xlabel 'Angle (degrees)'; ylabel 'Differentiated signal (\deltaN)/(\delta\theta)'

toc
>>>>>>> f8b824162267bde2f5d97a0a169a1ec31eab8856
