   tic
pm = 4;
swinit  = 0.25; % [T]
KRV = 5;

pzcut =  find(MxB(pm,:) <= swinit, 1, 'first')-1;

steps = [250,25,5,1];
many = [10,11,6,6];
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
        
            writer = linspace((swtest(ii-1).v(ww-1)*1e4),swtest(ii-1).v(ww)*1e4,many(ii-1))./1e4;        

        swtest(ii).v = writer; 
            
for ww = 2:length(swtest(ii).v)
    
    %%
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
    %%
    if RK >= KRV
        break
    end 

end

end

disp (['Values for a key ratio over 5 are between ', num2str(swtest(ii).v(ww)),' and ', num2str(swtest(ii).v(ww-1)),' [T]'])
toc
