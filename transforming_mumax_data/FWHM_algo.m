function [FWHMX,MLOC] = FWHM_algo(input,theta)
%FWHM Summary of this function goes here
%   Detailed explanation goes here
%%
temp.testline = abs(input);
    
    temp.maxfield = max(temp.testline);
    temp.TESTMAT = temp.testline == temp.maxfield;
    [K,L] = find(temp.TESTMAT);
    
    [temp.pks, temp.pkloc] = findpeaks(temp.testline);
    m = length(temp.pkloc);
    
    temp.indthet1 = [0,0];
    
    if input(1) ~= 0 
        
         % find the first zero before the halfway start (from algorithm)
       
        nn1 = ((m+1)/2) - ((m-1)/4);
        nn11 = nn1;
        
        while input(nn11) ~= 0
            nn11 = nn11-1;
        end 
        temp.indthet1(1) = nn11;
        
        % and find last after middle peak (simpler)
        nn2 = ((m+1)/2) + ((m-1)/4);
        temp.indthet1(2) = temp.pkloc(nn2) + find(temp.testline(nn2+1:length(temp.testline)) == 1 ,1,'first');
        
                 
    elseif input((length(input)/4)+1) ~= 0
        
        % find the first incidence of non zero
        temp.indthet1(1) = find(temp.testline >= 0.01*temp.maxfield ,1,'first'); 
        
        % and find first 0 after first full peak (simpler)
        temp.indthet1(2) = temp.pkloc(m/2) + find(temp.testline(temp.pkloc(m/2)+1:length(temp.testline)) == 0 ,1,'first') - 1;
    
    else
        
        temp.indthet1(1) = find(temp.testline >= 0.01*temp.maxfield ,1,'first'); 
        % Takes 1% of the max values so that only the switch is looked at in
        % the fwhm indexing - don't worry about the value here.

        %------------------------------------------------------ 
        %New bit! 
        n = temp.pkloc(m/4);
        temp.ttt = temp.testline == 0;
        temp.indthet1(2) = n + find(temp.ttt(n+1:length(temp.ttt)) == 1 ,1,'first');
        
        %temp.indthet1(2) = temp.indthet1(1)+find(temp.testline(temp.indthet1(1):length(temp.testline)) <= 0.01*temp.maxfield ,1,'first');       
        
    end
    %------------------------------------------------------    
         temp.tlfwhm = temp.testline(temp.indthet1(1):temp.indthet1(2));
         if any(temp.tlfwhm) == 0 
            FWHMX = 0;
         else
             temp.fractionx = temp.testline(L(1)) * (1/(2)); % currently a half,  1/sqrt(2) if on log plot
             temp.index1 = find(temp.tlfwhm >= temp.fractionx, 1, 'first');
             % Find where the data last rises above half the max.
             temp.index2 = find(temp.tlfwhm >= temp.fractionx, 1, 'last');
             temp.fwhm = temp.index2-temp.index1 + 1;% FWHM in indexes.
             % OR, if you have an x vector
             FWHMX =  theta(temp.fwhm) + theta(2) - theta(1);
         end
        
    
        MLOC = theta(L(1));
end

