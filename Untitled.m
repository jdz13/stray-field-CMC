                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    %% multi use analysis code for testing_rot_multi_distance.m

% option 1 -> Number of channels above threshold vs magnet size plotter
% option 2 -> Get a single scan - see which channels turn on when
% option 3 -> Proper analysis - with the search through the data - needs
% the analysis.keyratio etc - option 3 on testing_rot_multi_distance.m
% option 4 -> Proper analysis again - but this time for 'read first' data
% option 5 -> 
% option 6 ->
%






%%

<<<<<<< HEAD
mainoption = 6;
secondaryoption = 3;
=======
mainoption = 5;
secondaryoption = 5;
>>>>>>> f8b824162267bde2f5d97a0a169a1ec31eab8856
thirdoption = 0; 

%-------------------------------------------------------------------------
if mainoption == 1
    
    figno = 20; h = 16; % need to know which height to take this from.
    thrshld = linspace(0.1,0.5,5);
    figure(figno); clf; hold on; 
    xlabel 'Magnet size(m)' ; ylabel 'No of channels above threshold' ;
    
    for m = 1:length(thrshld)
        B = analysis.keyratio <= thrshld(m);
        adder = sum(B,3);
        plot(pm_cl, reshape(adder(:,h,:,:),[1,length(pm_cl)]))
       
    end
    
    legendCell = cellstr(num2str(thrshld', 'Threshold =%-g')); 
    legend(legendCell, 'Location', 'Northwest')
    title(['Sample to magnet surface distance = ',num2str(PZset(h)*100),'cm'])
    clear figno h m 
    
%-------------------------------------------------------------------------       
elseif mainoption == 2
    
    figno = 21; 
    nm = 15; % z height location in array;
    sp = 1; % sample space in array 
    pm = 4; % magnet size in array
    figure(figno); clf; subplot(2,1,1);hold on; subplot(2,1,2);hold on;
    for cx = 1:size(swfield,2)
    subplot(2,1,1); plot(thetad(nm,sp,pm).data, norm_vol_comp(nm,sp,pm).data(:,cx));
    subplot(2,1,2); plot(thetad(nm,sp,pm).data(2:length(thetad(nm,sp,pm).data)), diphis(nm,sp,pm).data(:,cx));
    end
    legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
    xlabel 'Angle (degrees)'; ylabel 'Differentiated N particles switched'
    title (['Differentiated number of switched for different angles at ' ,num2str(PZset(nm)*1000) ,'mm'])
    subplot(2,1,1)
    legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
    xlabel 'Angle (degrees)'; ylabel 'Normalised number of particles switched'
    title (['Normalised number of particles switched for different angles at ', num2str(PZset(nm)*1000), 'mm with sample space ', num2str(sampspac(sp)*1000), 'mm and magnet size ', num2str(pm_cl(pm)*100),'cm'])
    
    clear cx nm figno sp pm
    
    
%-------------------------------------------------------------------------    
elseif mainoption == 3
                         
    thrshld = linspace(0.1,0.5,5);
    sp = 1;
    pzindex = zeros(length(sp),length(swfield),length(pm_cl));
    
    for sf = 1:length(swfield)
        for pm = 1:length(pm_cl)
           ko = find(analysis.maxloc(sp,:,sf,pm) > 0,1,'last');
           if isempty(ko) == 1
               pzindex(sp,sf,pm) = 0;
           else 
               pzindex(sp,sf,pm) = ko;
           end
        end        
    end 
    
    clear sf pm ko
    
    % Need to choose which switching field to use, almost like a threshold.
    sfc = 4;
    plotter = zeros(1,length(pm_cl));
    figno = 22;
    
    
    for ko = 1:length(swfield)
    sfc = ko;  
    
    figure(figno+ko); clf; hold on;
    xlabel 'Magnet size(m)' ; ylabel 'No of channels above threshold' ;
    
    for m = 1:length(thrshld)
            for pm = 1:length(pm_cl)
                if pzindex(sp, sfc,pm) == 0
                    plotter(pm) = 0;
                    
                else
                    testline = reshape(analysis.keyratio(sp, pzindex(sp, sfc,pm), :, pm), [1,length(swfield)]); 
                    B = testline <= thrshld(m);
                    adder = sum(B);
                    plotter(pm) = adder;
                                                         
                end
                clear testline B adder
            end        
             plot(pm_cl,plotter)
             clear plotter
    end
    
    clear m pm
    
    legendCell = cellstr(num2str(thrshld', 'Threshold =%-g')); 
    legend(legendCell, 'Location', 'Northwest')
    title(['Using first instance of fields above ',num2str(swfield(sfc)),'T'])
    
    secondaryoption = 0;
    
    if secondaryoption == 1
    
        stravlin = ['C:\Users\Jake\Documents\MATLAB\Git\stray-field-CMC-cloned\Output_plots\Mag_size_param_sweep\3rd run threshold plots\max field = ' num2str(swfield(sfc)*1000) '_T.bmp'];
        fig = figure(figno+ko);
        set(gcf, 'Position', get(0, 'Screensize')); % Makes the figure full screen
        saveas(fig,stravlin);
        close
    
    elseif secondaryoption == 0
        
    else
        disp 'Invalid secondaryoption, please choose another'
    end
    
    
    end 
    
    clear option2 ko sp sfc figno
    
%-------------------------------------------------------------------------    
elseif mainoption == 4
   
     
    sp = 1;
    pzindex = zeros(length(sp),length(swfield),length(pm_cl));
    
    % secondaryoption = 1;
    % secondary option 1 - display graphs seperately for each sf
    % secondary option 2 - show all graphs overlaid
    % secondary option 3 - 2D graph - all data.
    
    for sf = 1:length(swfield)
        for pm = 1:length(pm_cl)
           ko = find(analysis.maxloc(sp,:,sf,pm) > 0,1,'last');
           if isempty(ko) == 1
               pzindex(sp,sf,pm) = 0;
           else 
               pzindex(sp,sf,pm) = ko;
           end
        end        
    end 
    
    clear sf pm ko
    
    % Need to choose which switching field to use, almost like a threshold.
    plotter = zeros(length(swfield),length(pm_cl));
    figno = 43;
    figure(figno); clf;
    
    for ko = 1:length(swfield)
    sfc = ko;    
    
     for pm = 1:length(pm_cl)
                if pzindex(sp, sfc,pm) == 0
                    plotter(sfc,pm) = 0;
                    
                else
                    
                    plotter(sfc,pm) = analysis.fwhmxxd(sp,pzindex(sp, sfc,pm),sfc,pm);
                                                         
                end
                
     end       
            clear pm
            
    end 

            
        if secondaryoption == 1 
            for sfc = 1:length(swfield)
                figure(figno+sfc-1); clf; hold on;
                xlabel 'Magnet size(m)' ; ylabel 'No of channels seperated below a threshold' ;
                plot(pm_cl,plotter(sfc,:))
                title(['Using first instance of fields above ',num2str(swfield(sfc)),'T'])
            end 
            clear sfc
       
        elseif secondaryoption == 2
             figure(figno); clf; hold on;
             for sfc = 1:length(swfield)
                  plot(pm_cl,plotter(sfc,:))
             end
             xlabel 'Magnet size(m)' ; ylabel 'Key ratio (FWHM/near-neighbour-distance)' ;
             title(['Using first instance of fields above ',num2str(swfield(sfc)),'T'])     
             legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
             clear sfc legendCell
                         
        elseif secondaryoption == 3
            figure(figno); 
            imagesc(pm_cl*100,swfield,flipud(plotter))
            xlabel 'Magnet size (cm)'
            ylabel 'Switching field (T)'
            title '2D map of ''read first'' configuration data'
            colorbar
            ax = gca; trial = linspace(ax.YLim(1),ax.YLim(2),length(swfield)+1);
            trial = trial - (trial(2)-trial(1))/2;    trial(1) = [];
            yticks(trial);    yticklabels((swfield));
              
        elseif secondaryoption == 4
            thrs =0.2; 
            B = sum((plotter <= thrs).*( ne(plotter,0)), 1); lims = [0,max(B)+1];
            figure(figno); plot(pm_cl*100, B); xlabel 'Magnet size (cm)';  
            ylabel(['Number of channels <', num2str(thrs)]); ylim(lims);
            title (['Number of channels under a threshold of ', num2str(thrs), ' for a read first configuration'])
      
        elseif secondaryoption == 5
            figure(figno); hold on; ylabel 'Number of channels <threshold'; xlabel 'Magnet size (cm)';
            title ('Number of channels over a threshold for a read first configuration')
            thrs = [2,3,4,5];
            B = zeros(size(pm_cl,2), size(thrs,2));
            for pm = 1:size(pm_cl,2)
                for count = 1:size(thrs,2) 
                    B(pm,count) = sum(plotter(:,pm) >= thrs(count)); 
                end 
            end 
            plot(pm_cl*100, B);   
                   % lims = [0,max(B)+1];        ylim(lims);
            legendCell = cellstr(num2str(thrs', 'Threshold =%-g')); legend(legendCell,'Location','Northwest')
            
        else 
            disp 'invalid secondaryoption - please choose another'
        end 
       
    % thirdoption = 0;
    % thirdoption 1 - save files (needs filename chosen)
    % thirdoption 0 - dont.
    
    if thirdoption == 1
    
        stravlin = ['C:\Users\Jake\Documents\MATLAB\Git\stray-field-CMC-cloned\Output_plots\Mag_size_param_sweep\3rd run threshold plots\max field = ' num2str(swfield(sfc)*1000) '_T.bmp'];
        fig = figure(figno+ko);
        set(gcf, 'Position', get(0, 'Screensize')); % Makes the figure full screen
        saveas(fig,stravlin);
        close
    
    elseif thirdoption == 0
        
    else
        disp 'Invalid thirdoption, please choose another'
    end
    
    clear option2 ko sp sfc figno secondaryoption
    
%-------------------------------------------------------------------------    
elseif mainoption == 5
    
    thrshld = linspace(0.1,0.5,5);
    sp = 1;
    pzindex = zeros(length(sp),length(swfield),length(pm_cl));
    
    plotter = zeros(length(pm_cl),length(swfield),length(thrshld));
    figno = 44;
    %figure(figno); clf;
    
    for sf = 1:length(swfield)
        for pm = 1:length(pm_cl)
           ko = find(analysis.maxloc(sp,:,sf,pm) > 0,1,'last');
           if isempty(ko) == 1
               pzindex(sp,sf,pm) = 0;
           else 
               pzindex(sp,sf,pm) = ko;
           end
           
           for m = 1:length(thrshld)
            
                if pzindex(sp, sf,pm) == 0
                    plotter(pm,sf,m) = 0;
                    
                else
                    testline = reshape(analysis.keyratio(sp, pzindex(sp, sf,pm), :, pm), [1,length(swfield)]); 
                    B = testline <= thrshld(m);
                    adder = sum(B);
                    plotter(pm,sf,m) = adder;                              
                
                clear testline B adder
                end        
             
            end
        end        
    end 
    
    clear sf pm ko m 
    
    
    % secondaryoption = 1;
    % secondary option 1 - display graphs seperately for each magsize (all
    % thresholds)
    % secondary option 2 - show all magsize graphs overlaid for one
    % threshold
    % secondary option 3 - 2D graph - all data. (one threshold
    
    
    if secondaryoption == 0
        
    elseif secondaryoption == 1
       
        for mgs = 1:length(pm_cl)
        
        figure(figno); clf; hold on;
        xlabel 'Field (T)'; ylabel 'Number of channels under threshold'
        
        for thr = 1:length(thrshld)
            plot(swfield, reshape(plotter(mgs, :,thr),[1,length(swfield)]))
            title (['number of channels below threshold for magsize ',num2str(pm_cl(mgs))])
        end
        figno = figno+1;
        end
        clear thr mgs
        
    elseif secondaryoption == 2 
        
        thrv = 2;
        figure(figno); clf; hold on;
        xlabel 'Field (T)'; ylabel (['Number of channels under threshold ', num2str(thrshld(thrv))]);
        
        for mgs = 1:length(pm_cl)
                       
            plot(swfield, reshape(plotter(mgs, :,thrv),[1,length(swfield)]))
            title 'number of channels below threshold for differing magnet size '
        
        end
        clear thr mgs
         legendCell = cellstr(num2str((pm_cl*100)', 'Mag size = %-g (cm)')); legend(legendCell)
        
    elseif secondaryoption == 3
        thrv = 2;
        figure(figno); clf;
        imagesc(swfield,pm_cl, reshape(plotter(:,:,thrv),[length(swfield),length(pm_cl)]))
        xlabel 'Field (T)'; ylabel 'magnet size (m)';
        title (['no of channels below threshold of ',num2str(thrshld(thrv))])
        clear thrv
      
        
    elseif secondaryoption == 4 
        
        thrv = 2; fieldpos = zeros(1, length(pm_cl));
        
        for mgs = 1:length(pm_cl)
                       
            line = reshape(plotter(mgs, :,thrv),[1,length(swfield)]);
            tet = max(line); ff = find(line == tet);
            fieldpos(mgs) = swfield(ff(1));
        clear line tet ff
        end
        
        figure(figno); clf; plot(pm_cl*100, fieldpos);
        xlabel 'Magnet size (cm)'; ylabel 'Field at max number of channels';
        title (['Field at max channel number under threshold ', num2str(thrshld(thrv)), ' for differing magnet size'])
        
        
    elseif secondaryoption == 5
        
        fieldpos = zeros(1, length(pm_cl)); figure(figno); clf; hold on;
        xlabel 'Magnet size (cm)'; ylabel 'Field at max number of channels';
        title 'Field at max channel number under different thresholds for differing magnet size'
        
        for thrv = 1:length(thrshld)
        for mgs = 1:length(pm_cl)
                       
            line = reshape(plotter(mgs, :,thrv),[1,length(swfield)]);
            tet = max(line); ff = find(line == tet);
            fieldpos(mgs) = swfield(ff(1));
        clear line tet ff
        end
        plot(pm_cl*100, fieldpos);
        end
        legendCell = cellstr(num2str(thrshld', 'Threshold =%-g')); legend(legendCell,'Location','Northwest')
        
        clear mgs thrv
         
    else     
        disp 'Invalid secondaryoption, please try another'
    end 
%-------------------------------------------------------------------------
elseif mainoption == 6
       tic
pm = 4;
swinit  = 0.25; % [T]
KRV = 5;

pzcut =  find(MxB(pm,:) <= swinit, 1, 'first')-1;

steps = [250,25,5,1];
clear swtestv

 e = 1e-16;
    idx = find(abs(Mdl_dtl(pm).topmagLinez - PZ(pzcut))<=e);
    
    NVC = zeros(2,length(theta));
    DNVC = zeros(2,length(theta)-1);
    FWHMX = [0,0]; MLOC = [0,0];
    
for pull = 1:length(theta)
   %Bxnew = Bobj.BXx.*sin(theta) + Bobj.BZx.*cos(theta);
   %Bynew = Bobj.BXy.*sin(theta) + Bobj.BZy.*cos(theta);
   Bznew = Bobj(pm).BXz.*sin(theta(pull)) + Bobj(pm).BZz.*cos(theta(pull)); 
    
   lkatpln = Bznew(:,:,idx);
   BZM = (lkatpln >= swinit) - (lkatpln <= -swinit);
   CM = BZM .* particle_loc;
   vc = sum(sum(CM));
   NVC(1,pull) = vc./control;
      
end

DNVC(1,:) = diff (NVC(1,:));
    
    ppp = 1;
    temp.testline = abs(DNVC(ppp,:));
    
    temp.maxfield = max(temp.testline);
    temp.TESTMAT = temp.testline == temp.maxfield;
    [K,L] = find(temp.TESTMAT);
    
    temp.indthet1 = [0,0];
    temp.indthet1(1) = find(temp.testline >= 0.01*temp.maxfield ,1,'first'); 
    % Takes 1% of the max values so that only the switch is looked at in
    % the fwhm indexing - don't worry about the value here.
    temp.indthet1(2) = temp.indthet1(1)+find(temp.testline(temp.indthet1(1):length(temp.testline)) <= 0.01*temp.maxfield ,1,'first');
        
             temp.tlfwhm = temp.testline(temp.indthet1(1):temp.indthet1(2));
             if any(temp.tlfwhm) == 0 
                FWHMX(ppp) = 0;
             else
             temp.fractionx = temp.testline(L(1)) * (1/(2)); % currently a half,  1/sqrt(2) if on log plot
             temp.index1 = find(temp.tlfwhm >= temp.fractionx, 1, 'first');
             % Find where the data last rises above half the max.
             temp.index2 = find(temp.tlfwhm >= temp.fractionx, 1, 'last');
             temp.fwhm = temp.index2-temp.index1 + 1;% FWHM in indexes.
             % OR, if you have an x vector
             FWHMX =  theta(temp.fwhm) + theta(2) - theta(1);
             end
        
    
        MLOC(ppp) = theta(L(1));
%-------------------------------------------------------------------------
   ppp = 2;
   swtestv = zeros(length(steps)+1, ((swinit/steps(1))*1e4)+1);
   swtest(1).v = [swinit,steps(1)./1e4];
   ww = 2;
for ii = 2:length(steps)+1
        if ii == 2
            writer = linspace((swtest(ii-1).v(ww-1)*1e4),swtest(ii-1).v(ww)*1e4,((swtest(ii-1).v(ww-1)-swtest(ii-1).v(ww))*1e4/steps(ii-1))+1)./1e4;
        else
            writer = linspace((swtest(ii-1).v(ww-1)*1e4),swtest(ii-1).v(ww)*1e4,((swtest(ii-1).v(ww-1)-swtest(ii-1).v(ww))*1e4/steps(ii-1))+2)./1e4;        
        end 
            
        swtest(ii).v = writer; 
            
for ww = 2:length(swtest(ii).v)
for pull = 1:length(theta)
   %Bxnew = Bobj.BXx.*sin(theta) + Bobj.BZx.*cos(theta);
   %Bynew = Bobj.BXy.*sin(theta) + Bobj.BZy.*cos(theta);
   Bznew = Bobj(pm).BXz.*sin(theta(pull)) + Bobj(pm).BZz.*cos(theta(pull)); 
    
   lkatpln = Bznew(:,:,idx);
   BZM = (lkatpln >= swtest(ii-1).v(ww)) - (lkatpln <= -swtest(ii-1).v(ww));
   CM = BZM .* particle_loc;
   vc = sum(sum(CM));
   NVC(2,pull) = vc./control;
      
end

    DNVC(2,:) = diff (NVC(2,:));
    
    temp.testline = abs(DNVC(ppp,:));
    
    temp.maxfield = max(temp.testline);
    temp.TESTMAT = temp.testline == temp.maxfield;
    [K,L] = find(temp.TESTMAT);
    MLOC(ppp) = theta(L(1));
 
    NN = abs(MLOC(2) - MLOC(1));
    RK = NN/FWHMX(1);
    
    if RK >= KRV
        break
    end 

end

end

disp (['Values for a key ratio over 5 are between ', num2str(swtest(ii).v(ww)),' and ', num2str(swtest(ii).v(ww-1)),' [T]'])

    toc
%-------------------------------------------------------------------------
else 
    disp 'Invalid mainoption - please choose another'   

end