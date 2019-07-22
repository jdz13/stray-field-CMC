
% Option 1 - Takes the data that was produced by the vsm_data_plotting_tool
% and plots all raw data on one graph. It then slope corrects the whole
% graph, using user input.

% Option 2 - Now has a user input data mechanism. First UI select the 
% SAMPLE data folder. Second select NO SAMPLE data folder. This will again 
% plot all of the difference data (the real sample data) and again correct
% for the slope using user input.


% option2 toggles the slope correction on/off. 1 = on, 0 = off.

% JDZ 20/06/2019

option = 2;
option2 = 1;

if option == 1
    

    clear datastr
    figno = 112;
    
    figure(figno); hold on 
    xlabel 'Applied field (Oe)'
    ylabel 'Raw moment (X) (memu)'
    title 'Raw Data plots'
    
    KP= 1;

        for tr = 1:size(zz,1)-2
    
            if isempty(zz(tr).data) == 1 
                continue
                
            else
                vhd = '.VHD';
                newSt = erase(zz(tr).name,vhd);
                newStr = strrep(newSt, '_', ' ');
    

                if size(zz(tr).data,2) == 13 
                    sec = 12;
                elseif size(zz(tr).data,2) == 11 
                    sec = 11;
                end
                
                plot(zz(tr).data(:,8), zz(tr).data(:,10))
                datastr(KP,1:size(zz(tr).data(:,10),1)) = zz(tr).data(:,10)';
                fieldstr(KP,1:size(zz(tr).data(:,8),1)) = zz(tr).data(:,8)';    
                KP = KP+1;

            end
        end 
    clear KP 

    if option2 == 0
        
    
    elseif option2 == 1
        hline = gline();
        w(1) = waitforbuttonpress; w(2) = waitforbuttonpress;
        coefficients = polyfit([hline.XData(2), hline.XData(1)], [hline.YData(2), hline.YData(1)], 1);
        gradient = coefficients (1);
        intercept = coefficients (2);
        checkergrad = (hline.YData(2)- hline.YData(1))/(hline.XData(2)- hline.XData(1));
        datacstr = (datastr) - (gradient.*fieldstr) ;
        
        figure(figno+1); hold on;
        xlabel 'Applied field (Oe)'
        ylabel 'Raw moment (X) (memu)'
        title 'Slope corrected raw data plots'
        
            for lp = 1:size(datacstr,1)
    
                plot((fieldstr(lp,1:length(nonzeros(datacstr(lp,:))))), nonzeros(datacstr(lp,:)))
    
            end

    else
        disp 'Invalid option 2, please input either 1 (Yes) or 0 (No).'
        
    end 
%%
elseif option == 2
    Sample_data = fol_data_ext_function();
    No_Sample = fol_data_ext_function();
    
    %%
    figno = 122;
    
    A = 6:2:34;
    B = [6:4:30,38:2:52]; %6:2:34; %
    figure(figno); clf; hold on 
    xlabel 'Applied field (Oe)'
    ylabel 'Raw moment (X) (memu)'
    
    for kpp = 1:length(A)
        
        subt_data(kpp,1:size(Sample_data(B(kpp)).data(:,10),1)) = (Sample_data(B(kpp)).data(:,10) - No_Sample(A(kpp)).data(:,10))';
        
        plot((Sample_data(B(kpp)).data(:,8))', subt_data(kpp,1:size(Sample_data(B(kpp)).data(:,10),1)))
        
        fldst(kpp,1:size(Sample_data(B(kpp)).data(:,10),1)) = (Sample_data(B(kpp)).data(:,8));
    end
    
    if option2 == 0
    
    elseif option2 == 1
        
        hline = gline();
        w(1) = waitforbuttonpress; w(2) = waitforbuttonpress;
        coefficients = polyfit([hline.XData(2), hline.XData(1)], [hline.YData(2), hline.YData(1)], 1);
        gradient = coefficients (1);
        intercept = coefficients (2);
        checkergrad = (hline.YData(2)- hline.YData(1))/(hline.XData(2)- hline.XData(1));
        
        subt_data = (subt_data) - (gradient.*fldst) ;
        
        figure(figno+1); hold on;
        xlabel 'Applied field (Oe)'
        ylabel 'Raw moment (X) (memu)'
        title 'Slope corrected (sample - no sample) consistency plots'
        
        for lp = 1:size(A,2)
            
            plot((fldst(lp,1:length(nonzeros(subt_data(lp,:))))), nonzeros(subt_data(lp,:)))
            
        end
        
    else
        disp 'Invalid option 2, please input either 1 (Yes) or 0 (No).'
    end
    
    
else
    disp 'Invalid option, please input another'
    
end 
