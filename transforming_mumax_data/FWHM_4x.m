function [FWHMX,MLOC] = FWHM_4x(input,theta)
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
        
        % find last after final peak (simpler)
        
        temp.indthet1(2) = temp.pkloc(m) + find(temp.testline(m+1:length(temp.testline)) == 0 ,1,'first');
        
        temp.fwhm = 2 * (temp.indthet1(2));% FWHM in indexes.
        % OR, if you have an x vector
        FWHMX =  theta(temp.fwhm) + theta(2) - theta(1);
        
                 
    elseif input(length(input)) ~= 0
        
        % find the first incidence of non zero
        temp.indthet1(1) = find(temp.testline >= 0.01*temp.maxfield ,1,'first'); 
        
        temp.fwhm = 2 * (length(input) - temp.indthet1(1));% FWHM in indexes.
        % OR, if you have an x vector
        FWHMX =  theta(temp.fwhm) + theta(2) - theta(1);
    
    else
        % find the first incidence of non zero
        temp.indthet1(1) = find(temp.testline >= 0.01*temp.maxfield ,1,'first'); 
        % find last after final peak (simpler)
        temp.indthet1(2) = temp.pkloc(m) + find(temp.testline(m+1:length(temp.testline)) == 0 ,1,'first');
    
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
        
    end
        MLOC = theta(L(1));
end

