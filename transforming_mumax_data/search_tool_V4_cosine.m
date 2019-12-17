 tic

%% ------------------------------------------------------------------------

clear SWres SH0 swinit SWnext count count2 NVC DNVC FWHMX MLOC SWnextpos

% We need to have defined the maxima's for each part
MxB = zeros(length(pm_cl),size(PZ,2));

for pm = 1:length(pm_cl)
    for ww = 1:size(PZ,2)
    
        e = 1e-16;
        idx = find(abs(Mdl_dtl(pm).topmagLinez - PZ(ww))<=e);
        lkpln = Bobj(pm).BZz(:,:,idx);
        MxB(pm,ww) = max(max(lkpln));
        
        
    end
end

%KRV = [5,4,3,2.5,2];
%PM = [1,2,3,4];
%RES = [0.15,0.2,0.25,0.3,0.35,0.4];%,0.45,0.5,0.55,0.6];
%f = [6,7];
ntestmax = 100;

ppp=2;

SWres = zeros(size(KRV,2),ntestmax,length(PM),length(RES)); 
FWHMres = SWres; Bset = SWres;

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
            SWnext = 0; % initial condition
            count = 1;
            tmps = [0,1];

            while abs(tmps(1) - tmps(2)) > 1e-3 && SH0 > MxB(pm,length(MxB(pm,:))) && tmps(2) ~= 0

                % Find where that sits in space (Pz)
                pzcut =  find(MxB(pm,:) <= SH0, 1, 'first')-1;
                e = 1e-14; % tolerance - numerical rounding 
                idx = find(abs(Mdl_dtl(pm).topmagLinez - PZ(pzcut))<=e);

                clear variable % fill this with whatever needs refreshing through each run. 

                % Initialise the variables you're going to want 
                NVC = zeros(2,length(theta)); % Normalised volume comparison
                DNVC = zeros(2,length(theta)-1); % Differentiated NVC 
                FWHMX = [0,0]; MLOC = [0,0]; % Full width half max and Maxima location

                %% ------------------------------------------------------------------------
                % Look for information about this initial condition


                    for pull = 1:length(theta) 
                       % Use rotation matricies to find the Z component 
                       Bznew = Bobj(pm).BXz.*sin(theta(pull)) + Bobj(pm).BZz.*cos(theta(pull)); 
                       % Look at the plane where the sample sits 
                       lkatpln = Bznew(:,:,idx);
                       %Find out how much of this areas is above or below the threshold
                       BZM = (lkatpln >= swinit) - (lkatpln <= -swinit);
                       % Correlate with where the particles actually are in the world
                       CM = BZM .* particle_loc;
                       % Find a qualitative number for how much is 'on'
                       vc = sum(sum(CM));
                       % Compare this to how many are in the sample space 
                       NVC(1,pull) = vc./control;
                    end

                % Differentiate the value to see the impulse (what we'd measure)
                DNVC(1,:) = diff (NVC(1,:));

                [FWHMX(1),MLOC(1)] = FWHM_algo(DNVC(1,:),theta);
                %RKout = zeros(1,100); 

                SWnextpos = MLOC(1)+(KRV(count2)*FWHMX(1));  

                SWnext = swinit*cos(SWnextpos);

                SWres(count2,count+1,pmcount,rescount) = SWnext;
                FWHMres(count2,count+1,pmcount,rescount) = FWHMX(1); % if +1 not needed then can use FWHMres(:,1,:,:) = [];

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

% SWres(:,nnz(SWres(size(SWres,1),:))+1:size(SWres,2)) = [];
% SWres(SWres==0)=nan; 

trialler = SWres;
 trialler(:,nnz(trialler(size(trialler,1),:))+1:size(trialler,2)) = [];
 trialler(trialler==0)=nan; 
%%
%figure(5); plot(SWres','x--')
%xlabel 'Channel number (n)'; ylabel 'Channel strength (T)'
%legendCell = cellstr(num2str(KRV', 'KRV =%-g ')); legend(legendCell)

toc
%figure(50); plot(trialler','x--')
%xlabel 'Channel number (n)'; ylabel 'Channel strength (T)'
%legendCell = cellstr(num2str(KRV', 'KRV =%-g ')); legend(legendCell)

plotter14 = zeros(length(KRV), length(PM),length(RES));

for i = 1:length(KRV)
    for j = 1:length(PM)
        for k = 1:length(RES)
               plotter14(i,j,k) = nnz(SWres(i,:,j, k));
        end 
    end 
end

plotyn = 0;

if plotyn == 0

elseif plotyn == 1 

    figure(f(1)); clf;
    for jj = 1:length(RES)
        subplot(2,2,jj); imagesc(PM, KRV,plotter14(:,:,jj));
        xlabel 'PM size [cm]'; ylabel 'KRV'; title (['Start field = ', num2str(RES(jj)),'T'])
        colorbar
    end
    %%
    figure(f(2)); clf;
    for jj = 1:length(RES)
        subplot(2,5,jj); imagesc(PM, KRV,plotter14(:,:,jj));
        xlabel 'PM size [cm]'; ylabel 'KRV'; title (['Start field = ', num2str(RES(jj)),'T'])
        caxis([min(plotter14,[],'all'),max(plotter14,[],'all')]); colorbar
        ax = gca; trial = linspace(ax.YLim(1),ax.YLim(2),length(KRV)+1);
        trial = trial - (trial(2)-trial(1))/2;    trial(1) = [];
        yticks(trial);    yticklabels((fliplr(KRV)));
    end

else 
    disp 'Put in a proper value for if you want a plot (plotyn)'

end 
