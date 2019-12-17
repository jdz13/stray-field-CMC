%% Both a running and analysis script for rotating magnetic driving field. 

% Load in the 'multichannel_multisize_signal_results_05_03_19.mat' for
% a base data. Else run your own option 1 scan. 

% Option 1 allows you to run the simulation for whatever parameters you choose.
% ** runs off the rot_V2.m function 
% Option 2 plots the channel scan for each height PZ. 
% Option 3 peforms the FWHM and distance to next channel analysis.
% Option 4 outputs the FWHM vs PZ plots for each different sample space.  
% **requires option 3 to be run first
% Option 5 displays FWHM for each channel - imagesc plots. PZ/swfield axes
% **requires option 3 to be run first
% **can be adapted for ratio/max-loc/separation etc. Any analysis variable.
% Option 6 - not really sure what I was attempting here.... or if it works
% Option 7 produces a graph showing the field ranges possible for a given
% set of magnet sizes over a particular distance. Needs these inputs. 
% Option 8 allows you to plot all of the imagesc plots of ratio_key on one
% set of subplots 
% ** requires option 3 to be ran first
% ** needs inputs - check the notation at the beginning of the code

% JDZ March 2019 (jdz25@cam.ac.uk)

tic

%-------------------------------------------------------------------------
option = 3;
%-------------------------------------------------------------------------

if option == 1 

    clear

    howlong(1) = datetime;
    count = 0;
    PZset = [linspace(1e-3,8.5e-3,6),linspace(1e-2,50e-2,50)];
    
PZ = PZset; %0.042; %linspace(5e-2,3e-2,21);n
swfield = [2005,1703,1385, 1099, 798, 579, 427, 261, 150, 40]./1e4; %[linspace(0.4, 0.1,7),linspace(0.08,0.02,4)];
theta_end = pi; theta_n =  361;
theta = linspace(0,theta_end, theta_n); 
sampspac = 1e-3; % sampspac =linspace(2.5e-4,2e-3,8);
pm_cl = linspace(1e-2,4e-2,7); %1.79e-2;  % If cuboidal, this is the dimension [m]


for sc = 1:length(sampspac)
for k = size(PZ,2):-1:1    % Backwards!
  norm_vol_comp(k,sc).data = zeros(theta_n, size(swfield,2));
  thetad(k,sc).data = zeros(1,theta_n);
  diphis(k,sc).data = zeros(theta_n-1, size(swfield,2));
end
end

for pm = 1:length(pm_cl)
    PZ = PZset + (pm_cl(pm)/2);
for sp = 1:length(sampspac)
for nm = 1:size(PZ,2)
    [norm_vol_comp(nm,sp,pm).data, thetad(nm,sp,pm).data, diphis(nm,sp,pm).data] = rot_V2(PZ(nm),swfield,theta,sampspac(sp),pm_cl(pm));
    
    count = count+1;
    if rem(count, 20) == 0
        disp(count)
    end
end
end 
end
toc
clear nm sp pm sc k
    howlong(2) = datetime;
d = date; curfol = pwd; str = [curfol,'\magsize_results_',d];
save(str);

%% ------------------------------------------------------------------------
tic

elseif option == 2
    
    % option2 selection criteria.
    % 1 = PM distance scans
    % 2 = magsize scans
    
    option2 = 2;
    
        if option2 == 1
            
figno = 75;
for nm = 1:size(PZ,2)
    figure(figno); clf; subplot(2,1,1);hold on; subplot(2,1,2);hold on;
    for cx = 1:size(swfield,2)
    subplot(2,1,1); plot(thetad(nm).data, norm_vol_comp(nm).data(:,cx));
    subplot(2,1,2); plot(thetad(nm).data(2:length(thetad(nm).data)), diphis(nm).data(:,cx));
    end
    legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
    xlabel 'Angle (degrees)'; ylabel 'Differentiated N particles switched'
    title (['Differentiated number of switched for different angles at ' ,num2str(PZ(nm)*1000) ,'mm'])
    subplot(2,1,1)
    legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
    xlabel 'Angle (degrees)'; ylabel 'Normalised number of particles switched'
    title (['Normalised number of particles switched for different angles at ', num2str(PZ(nm)*1000), 'mm'])

   figno = figno+1;
end
        elseif option2 == 2
            figno = 75;
            nm = 30; % choose the z height 
            sp = 1; % choose the sample size 
            
            for pm = 1:length(pm_cl)
                figure(figno); clf; subplot(2,1,1);hold on; subplot(2,1,2);hold on;
                for sf = 1:size(swfield,2)
                    subplot(2,1,1); plot(thetad(nm,sp,pm).data, norm_vol_comp(nm,sp,pm).data(:,sf));
                    subplot(2,1,2); plot(thetad(nm,sp,pm).data(2:length(thetad(nm,sp,pm).data)), diphis(nm,sp,pm).data(:,sf));
                end
                legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
                xlabel 'Angle (degrees)'; ylabel 'Differentiated N particles switched'
                title (['Differentiated number of switched for different angles for a magnet of ' ,num2str(pm_cl(pm)*1000) ,'mm'])
                subplot(2,1,1)
                legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
                xlabel 'Angle (degrees)'; ylabel 'Normalised number of particles switched'
                title (['Normalised number of particles switched for different angles for a magnet of  ', num2str(pm_cl(pm)*1000), 'mm'])
                
                figno = figno+1;
            end 
            
        end 
clear nm cx mgs sf pm option2


toc

%% ------------------------------------------------------------------------

elseif option == 3
    
    tic
    
    analysis.fwhmxx = zeros(length(sampspac),size(PZ,2),length(swfield),length(pm_cl));
    analysis.maxloc = zeros(length(sampspac),size(PZ,2),length(swfield),length(pm_cl));
    
    if length(pm_cl) == 1
        analysis.fwhmxx(:,:,:,2) = 0;
        analysis.maxloc(:,:,:,2) = 0;       
    end 
    
for sp = 1:length(sampspac)
for nm = 1:size(PZ,2)
for sf = 1:length(swfield)
for pm = 1:length(pm_cl)
    
     temp.testline = abs(diphis(nm,sp,pm).data(:,sf));
    
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
                analysis.fwhmxx(sp,nm,sf,pm) = 0;
             else
             temp.fractionx = diphis(nm,sp,pm).data(K(1),sf) * (1/(2)); % currently a half,  1/sqrt(2) if on log plot
             temp.index1 = find(temp.tlfwhm >= temp.fractionx, 1, 'first');
             % Find where the data last rises above half the max.
             temp.index2 = find(temp.tlfwhm >= temp.fractionx, 1, 'last');
             temp.fwhm = temp.index2-temp.index1 + 1;% FWHM in indexes.
             % OR, if you have an x vector
             analysis.fwhmxx(sp,nm,sf,pm) =  theta(temp.fwhm) + theta(2) - theta(1);
             end
        
    
    analysis.maxloc(sp,nm,sf,pm) = theta(K(1));
    
end
end 
end 
end

clear sp nm sf K L pm
analysis.fwhmxx(analysis.fwhmxx==90) = 0; %removes all the large switches to 90deg fwhm to 0.
analysis.fwhmxxd = rad2deg(analysis.fwhmxx);

toc

analysis.dmax = diff(analysis.maxloc,1,3); analysis.dmax(analysis.dmax > pi/2) = 0 ;
analysis.nneig = nearneigh(analysis.maxloc,3); analysis.nneig(analysis.nneig > pi/2) = 0 ;

analysis.keyratio = analysis.fwhmxx./analysis.nneig; analysis.keyratio(analysis.keyratio == inf) = NaN;

analysis.ratiokey = analysis.nneig./analysis.fwhmxx;
%% ------------------------------------------------------------------------
elseif option == 4
    mgs = 1; % choose the magnet size variable to lock
    
    for jkk = 1:length(sampspac)
    figure(jkk+5); clf; hold on
    for hk = 1:length(swfield)
    plot(PZ, analysis.fwhmxxd(jkk,:,hk,mgs))
    end
    legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
    title (['FWHM plot vs PZ for sample radius = ', num2str(sampspac(jkk)*1000), 'mm']); axis([0.03,0.05,0,10]); % hline(0,'k--','0 deg'); 
    ylabel 'FWHM (degrees)'; xlabel 'Sample - driver magnet distance (m)'
    end
 
clear jkk hk mgs

%% ------------------------------------------------------------------------
elseif option == 5
    
    sel = 1; % the sample space that you would like to view
    mgs = 4;
    
    % Need to choose both axes
    % Option 1 = Sample space
    % Option 2 = Pz - magnet-sample height
    % Option 3 = Switching field 
    % Option 4 = PM cubic length
    
    temp.axtitles = {'Sample Space (m)','PZ (m)', 'switching field (T)', 'Magsize (cubic length) (m)'};
    temp.lengths = [length(sampspac),length(PZ),length(swfield),length(pm_cl)];
    
    Xop = 3;       Yop = 4;   
    
    figure(16); imagesc(swfield,PZ, (squeeze((analysis.ratiokey(sel,:,:,mgs))))); 
    ax = gca; trial = linspace(ax.XLim(1),ax.XLim(2),length(swfield)+1);
    trial = trial - (trial(2)-trial(1))/2;    trial(1) = [];
    xticks(trial);    xticklabels(fliplr(swfield));
    xlabel 'swfield (T)';ylabel 'Pz (m)';
    title (['Ratio of FWHM to Near Neighbour plot for each channel with sample radius = ', num2str(sampspac(sel)*1000), 'mm'])
    colorbar

clear sel mgs trial
   
%% ------------------------------------------------------------------------
elseif option == 6
    
    % Need to choose both axes
    % Option 1 = Sample space
    % Option 2 = Pz - magnet-sample height
    % Option 3 = Switching field 
    % Option 4 = PM cubic length
    
    temp.axtitles = {'Sample Space (m)','PZ (m)', 'switching field (T)', 'Magsize (cubic length) (m)'};
    temp.lengths = [length(sampspac),length(PZ),length(swfield),length(pm_cl)];
    
    Xop = 3;       Yop = 4;   
    for fig = 6:10
    figure(fig); 
    imagesc(swfield, (pm_cl), (rot90(reshape(rad2deg(analysis.nneig(1,fig-5,:,:)),[temp.lengths(Xop),temp.lengths(Yop)]))));
    title (['Ratio of fwhm/nearest-neighbour-distance at sample to magnet surface distance = ', num2str(PZset(fig-5))])
    xlabel (temp.axtitles(Xop)); ylabel(temp.axtitles(Yop));
    ax = gca; trial = linspace(ax.XLim(1),ax.XLim(2),length(swfield)+1);
    trial = trial - (trial(2)-trial(1))/2;    trial(1) = [];
    xticks(trial);    xticklabels(fliplr(swfield));
    
    ax = gca; trial = linspace(ax.YLim(1),ax.YLim(2),length(pm_cl)+1);
    trial = trial - (trial(2)-trial(1))/2;    trial(1) = [];
    yticks(trial);    yticklabels(fliplr(pm_cl));
    colorbar
    
    end
    clear Xop Yop
 
%% ------------------------------------------------------------------------
elseif option == 7
    
    graphonly = 1;
    
    if graphonly == 0
    
    op7.magsizes = linspace(1e-2,4e-2,7); % should I just use pm_cl? to keep parameters
    op7.mxds = 1e-1; % define the max magnet surface to sample distance - lowest field
    op7.mnds = 5e-3; % same but minimum - highest field 
    op7.PZuset = [op7.mxds, op7.mnds];
    
    
    op7.plane_N = [200,200];
    op7.plane_lengths =  [5e-3,5e-3]; % [m] plus/minus extents to 
    op7.cell_size = 2.*op7.plane_lengths./op7.plane_N;

    op7.px = linspace(-op7.plane_lengths(1)+op7.cell_size(1)/2,op7.plane_lengths(1)-op7.cell_size(1)/2,op7.plane_N(1));
    op7.py = linspace(-op7.plane_lengths(2)+op7.cell_size(2)/2,op7.plane_lengths(2)-op7.cell_size(2)/2,op7.plane_N(2));    
    
    op7.displa = [0,0,0]; 
    op7.Msatpm = 1.27e6;
    l = 1;
    mxmdplt = zeros(1,length(op7.magsizes));
    mnmdplt = mxmdplt;
    
    mxmxplt = mxmdplt;
    mnmxplt = mxmdplt;
    
    mxmnplt = mxmdplt;
    mnmnplt = mxmdplt;
       
    
        for kkl = 1:length(op7.magsizes)
            mgzs = [op7.magsizes(kkl),op7.magsizes(kkl),op7.magsizes(kkl)];
            
            op7.PZuser = op7.PZuset + (op7.magsizes(kkl)/2);
                                   
            [op7o(l).akoun, pzplot, op7o(l).max] = PMnotFD(mgzs, op7.px, op7.py, op7.PZuset,op7.Msatpm, op7.displa);
                    
            mxmdplt(l) = op7o(l).max(1); mnmdplt(l) = op7o(l).max(2);
            
            mxmxplt(l) = max(max(op7o(l).akoun(1).HzAkoun)); mnmxplt(l) = max(max(op7o(l).akoun(2).HzAkoun));
            
            mxmnplt(l) = mean(mean(op7o(l).akoun(1).HzAkoun)); mnmnplt(l) = mean(mean(op7o(l).akoun(2).HzAkoun));
                        
            l = l+1;
                        
        end 
    
        elseif graphonly == 1
        secop = 4;
        figno = 67;
        
        %secop1 - plot maxfield
        %secop2 - plot midfield
        %secop3 - plot avfield 
        %secop4 - plot all max options to compare
        
        if secop == 1
            figure(figno); clf; hold on;
            plot(op7.magsizes.*100, mxmxplt, op7.magsizes.*100, mnmxplt)
            xlabel 'Magnet cuboidal length (cm)'
            ylabel 'Field (T)'
            legendCell = cellstr(num2str((op7.PZuset.*100)', 'PZ =%-g (cm)')); legend(legendCell,'Location','Northwest')
            title 'Maximum and minimum fields for different magnet sizes - max method'
             
        elseif secop == 2
            figure(figno); clf; hold on;
            plot(op7.magsizes.*100, mxmdplt, op7.magsizes.*100, mnmdplt)
            xlabel 'Magnet cuboidal length (cm)'
            ylabel 'Field (T)'
            legendCell = cellstr(num2str((op7.PZuset.*100)', 'PZ =%-g (cm)')); legend(legendCell,'Location','Northwest')
            title 'Maximum and minimum fields for different magnet sizes - middle method'
            
        elseif secop == 3
            figure(figno); clf; hold on;
            plot(op7.magsizes.*100, mxmnplt, op7.magsizes.*100, mnmnplt)
            xlabel 'Magnet cuboidal length (cm)'
            ylabel 'Field (T)'
            legendCell = cellstr(num2str((op7.PZuset.*100)', 'PZ =%-g (cm)')); legend(legendCell,'Location','Northwest')
            title 'Maximum and minimum fields for different magnet sizes - mean method'
            
        elseif secop == 4
            figure(figno); clf; hold on;
            plot(op7.magsizes.*100, mnmxplt,'ro', op7.magsizes.*100, mnmdplt,'bx', op7.magsizes.*100, mnmnplt,'g--')
            xlabel 'Magnet cuboidal length (cm)'
            ylabel 'Field (T)'
            title (['Maximum fields for different magnet sizes - all methods at ', num2str(op7.mnds*1000),'mm'])
            legend ('Max','Middle','Mean', 'Location', 'Northwest')
            
        else 
            disp 'Invalid secop - choose another'
        end 
    end 
   
%% ------------------------------------------------------------------------
elseif option == 8
    
    % Things that need to be input
    
    figno1 = 16; % what figure number do you want this to be?
    sel = 1; % choose which sample space to look over (index)
    hm = [2,2]; % How many - what size subplot space do you want? 
    if hm(1)*hm(2) < length(pm_cl)
        disp 'Check you''re number of subplots - it doesn''t appear to be enough!'
    else 
    
    figure(figno1);clf
     
 for mgs = 1:length(pm_cl)
    
     subplot(hm(1),hm(2),mgs);
    
    imagesc(swfield,PZ, (squeeze((analysis.ratiokey(sel,:,:,mgs))))); 
    ax = gca; trial = linspace(ax.XLim(1),ax.XLim(2),length(swfield)+1);
    trial = trial - (trial(2)-trial(1))/2;    trial(1) = [];
    xticks(trial);    xticklabels(fliplr(swfield));
    xlabel 'swfield (T)';ylabel 'Pz (m)';
    title (['Ratio plot for magnet OD = ', num2str(pm_cl(mgs)*1000), 'mm'])
    colorbar; %caxis([2,5])
 end 
    end
    
%% ------------------------------------------------------------------------
else 
    disp 'Invalid option - please choose another'
    
end 
   