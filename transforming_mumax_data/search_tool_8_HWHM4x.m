function [varst, SWres,Bset,FWHMres,ind1res,ind2res,MxB] = search_tool_8_HWHM4x(KRV,PM,RES,pm_cl,Mdl_dtl,Bobj,theta,particle_loc,control,con)

tic

%% ------------------------------------------------------------------------

clear SWres SH0 swinit SWnext count count2 NVC DNVC FWHMX MLOC SWnextpos

% We need to have defined the maxima's for each part
MxB = zeros(length(pm_cl),length(Mdl_dtl(1).topmagLinez));

for pm = 1:length(pm_cl)
    for ww = 1:length(Mdl_dtl(pm).topmagLinez)
        
        if Mdl_dtl(pm).topmagLinez(ww) > 0
            lkpln = Bobj(pm).BZz(:,:,ww);
            MxB(pm,ww) = max(lkpln,[],'all');    
        else 
            
        end 
    end
end

ntestmax = 100;

SWres = zeros(size(KRV,2),ntestmax,length(PM),length(RES)); 
FWHMres = SWres; ind1res = SWres; ind2res = SWres; Bset = SWres;

for gg = 1:length(RES)
    SWres(:,1,:,gg) = RES(gg);
end

for rescount = 1:length(RES)

    res = RES(rescount);

    for pmcount = 1:length(PM)

        pm = PM(pmcount);

        for  count2 = 1:size(KRV,2)

            % Define the initial conditions
            SH0 = 1.05*res;
            swinit = res; % What's the max channel value?
            count = 1;
            tmps = [0,1];

            while abs(tmps(1) - tmps(2)) > 1e-4 && SH0 > MxB(pm,length(MxB(pm,:))) && tmps(2) ~= 0

                % Find where that sits in space (Pz)
               pzcut =  find(MxB(pm,:) >= SH0, 1, 'last');
               
                clear variable % fill this with whatever needs refreshing through each run. 

                % Initialise the variables you're going to want 
                NVC = zeros(2,length(theta)); % Normalised volume comparison
                FWHMX = [0,0]; MLOC = [0,0]; % Full width half max and Maxima location

                %% ------------------------------------------------------------------------
                % Look for information about this initial condition


                    for pull = 1:length(theta) 
                       % Use rotation matricies to find the Z component 
                       Bznew = Bobj(pm).BXz.*sin(theta(pull)) + Bobj(pm).BZz.*cos(theta(pull)); 
                       % Look at the plane where the sample sits 
                       lkatpln = Bznew(:,:,pzcut);
                       %Find out how much of this areas is above or below the threshold
                       BZM = (lkatpln >= swinit) - (lkatpln <= -swinit);
                       % Correlate with where the particles actually are in the world
                       CM = BZM .* particle_loc;
                       % Find a qualitative number for how much is 'on'
                       vc = sum(sum(CM));
                       % Compare this to how many are in the sample space 
                       NVC(1,pull) = vc./control;

                    end

                [FWHMX(1),MLOC(1),indout] = FWHMNVC(NVC(1,:),theta,con);
                %RKout = zeros(1,100); 

                SWnextpos = MLOC(1)+(KRV(count2)*FWHMX(1));  

                SWnext = swinit*cos(SWnextpos);

                SWres(count2,count+1,pmcount,rescount) = SWnext;
                FWHMres(count2,count+1,pmcount,rescount) = FWHMX(1); % if +1 not needed then can use FWHMres(:,1,:,:) = [];
                ind1res(count2,count+1,pmcount,rescount)= indout(1);
                ind2res(count2,count+1,pmcount,rescount) = indout(2);
                Bset(count2,count,pmcount,rescount) = MxB(pm,pzcut);
                % manipulate the results to run the next leg

                tmps(1) = swinit; 

                    if SWnextpos >= pi/2
                         tmps(2) = 0;
                    else
                        tmps(2) = SWnext;
                    end

                SH0 = (tmps(1)+tmps(2))/2;
                swinit = tmps(2); 

                count = count +1;

                disp (['count = ', num2str(count),', range = ', num2str(tmps)])

            end
        end
    end
end 

varst.KRV = KRV;
varst.PM = PM;
varst.RES = RES;
varst.CON = con;
varst.theta = theta;
varst.MxB = MxB;
varst.timer = toc;

end 