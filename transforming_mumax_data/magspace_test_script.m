
option1 = 0;

% first thing - do you need to import the data?
% 1 = yes
% 0 = no 

       
       
if option1 == 1
       tic
       
       clear
    
       how_many = 4;
       mag_IDs = [6,6,6,6]*1e-3;
       mag_ODs = [10,20,30,40]*1e-3;% this should be the amount of different magnet sizes that you want to include in the code 
        
       if  length(mag_IDs)~=length(mag_ODs) || length(mag_IDs)~= how_many
          disp ('Inner and outer magnet diameter variables inconsistent')
       end
       
       for count = 1:how_many
    % first select X data, then Y, then Z;
    disp (['X',num2str(count)])
    [Bobj(count).BXobject,Bobj(count).BXx,Bobj(count).BXy,Bobj(count).BXz] = Mumax_data_extractor();
    disp (['Y',num2str(count)])
    [Bobj(count).BYobject,Bobj(count).BYx,Bobj(count).BYy,Bobj(count).BYz] = Mumax_data_extractor();
    disp (['Z',num2str(count)])
    [Bobj(count).BZobject,Bobj(count).BZx,Bobj(count).BZy,Bobj(count).BZz] = Mumax_data_extractor();
       end 
      
%% ------------------------------------------------------------------------
       
elseif option1 == 0

% Bringing in the world definitions 

Mdl_dtl.CellNo = size(Bobj(1).BXx); % Finding out how many cells there are in the world in each dimension
Mdl_dtl.gridsize = Bobj(1).BXobject.GridSize'; % Finding out how big each cell is in each dimension;
Mdl_dtl.wrldSz = Mdl_dtl.gridsize.* Mdl_dtl.CellNo; % Total worldsize
% This should be consistent for all data, so the first structure variable
% is used (there should always be a first one!). This also holds for all of
% the purelines below as well as the extents.  

Mdl_dtl.extents = Mdl_dtl.wrldSz./2;

Mdl_dtl.purelinex = linspace(-Mdl_dtl.extents(1)+Mdl_dtl.gridsize(1)/2,Mdl_dtl.extents(1)-Mdl_dtl.gridsize(1)/2,Mdl_dtl.CellNo(1));
Mdl_dtl.pureliney = linspace(-Mdl_dtl.extents(2)+Mdl_dtl.gridsize(2)/2,Mdl_dtl.extents(2)-Mdl_dtl.gridsize(2)/2,Mdl_dtl.CellNo(2));
Mdl_dtl.purelinez = linspace(-Mdl_dtl.extents(3)+Mdl_dtl.gridsize(3)/2,Mdl_dtl.extents(3)-Mdl_dtl.gridsize(3)/2,Mdl_dtl.CellNo(3));

%This is where it may need changing though. Will have this update from
%above - it's easiest to make any user input all in one place. 
Mdl_dtl.ID = mag_IDs;
Mdl_dtl.OD = mag_ODs; % This variable definitely changes.

for count = 1:how_many
    
Mdl_dtl(count).displace = [0,0,0.5*(-Mdl_dtl(1).OD(count) + Mdl_dtl(1).wrldSz(3))];

Mdl_dtl(count).cntrmagLinex = Mdl_dtl(1).purelinex + Mdl_dtl(count).displace(1);
Mdl_dtl(count).cntrmagLiney = Mdl_dtl(1).pureliney + Mdl_dtl(count).displace(2);
Mdl_dtl(count).cntrmagLinez = Mdl_dtl(1).purelinez + Mdl_dtl(count).displace(3);

Mdl_dtl(count).topmagLinez = Mdl_dtl(count).cntrmagLinez - Mdl_dtl(1).OD(count)/2 + Mdl_dtl(1).gridsize(3)/2;
end 

%%
tic

% Start defining and manipulating the probe field 

% Begin with the definition


% Now look at applying the logical matricies. 
% Supply variables telling it what went in. This is where the real data
% get's created. 

pm_cl = mag_ODs;
sampspac = 1e-3;
particle_loc = plane_mask(Mdl_dtl(1).cntrmagLinex,Mdl_dtl(1).cntrmagLiney,sampspac);  % I think this is okay - X and Y should be the same as all should be using the same zero datum.
control = sum(sum(particle_loc));

swfield = linspace(2500,250,10)./1e4;
%swfield = [2005,1703,1385, 1099, 798,700, 579, 427, 261, 150, 40]./1e4;
theta = linspace(0,2*pi,361);
PZ = linspace(1e-3,60e-3,237);

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
    idx = find(abs(Mdl_dtl(pm).topmagLinez - PZ(nm))<=e);
    
for pull = 1:length(theta)
   %Bxnew = Bobj.BXx.*sin(theta) + Bobj.BZx.*cos(theta);
   %Bynew = Bobj.BXy.*sin(theta) + Bobj.BZy.*cos(theta);
   Bznew = Bobj(pm).BXz.*sin(theta(pull)) + Bobj(pm).BZz.*cos(theta(pull)); 
    
   lkatpln = Bznew(:,:,idx);
        for zzz = 1:length(swfield)
            rotfield.bzmask_sf1(:,:,zzz) = (lkatpln >= swfield(zzz)) - (lkatpln <= -swfield(zzz));
            rotfield.corrmask_sf1(:,:,zzz) = rotfield.bzmask_sf1(:,:,zzz) .* particle_loc;
            vol_compa = sum(sum(rotfield.corrmask_sf1(:,:,zzz)));
            norm_vol_comp(nm,sp,pm).data(pull,zzz) = vol_compa./control;
        end

end
diphis(nm,sp,pm).data = diff(norm_vol_comp(nm,sp,pm).data);
end
end 
end

else
        disp 'Please input either 1(yes) or 0(no) for option1'
end

toc