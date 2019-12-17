function [FWHMX,MLOC] = FWHM(input,theta)
%FWHM Summary of this function goes here
%   Detailed explanation goes here
%%
temp.testline = abs(input);
    
    temp.maxfield = max(temp.testline);
    temp.TESTMAT = temp.testline == temp.maxfield;
    [K,L] = find(temp.TESTMAT);
    
    temp.indthet1 = [0,0];
    temp.indthet1(1) = find(temp.testline >= 0.01*temp.maxfield ,1,'first'); 
    % Takes 1% of the max values so that only the switch is looked at in
    % the fwhm indexing - don't worry about the value here.
    
    %------------------------------------------------------ 
    %New bit! 
    n = temp.indthet1(1);
    while any(temp.testline(n+1:n+4)) ~= 0
        n = n+1;
    end 
    temp.indthet1(2) = n;
    
    %temp.indthet1(2) = temp.indthet1(1)+find(temp.testline(temp.indthet1(1):length(temp.testline)) <= 0.01*temp.maxfield ,1,'first');
    
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

