function [particle_loc, control, Bobj, Mdl_dtl] = get_data_in (how_many, mag_IDs, mag_ODs) 
    
    sampspac = 1e-3;

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


    particle_loc = plane_mask(Mdl_dtl(1).cntrmagLinex,Mdl_dtl(1).cntrmagLiney,sampspac);  % I think this is okay - X and Y should be the same as all should be using the same zero datum.
    control = sum(sum(particle_loc));

end

