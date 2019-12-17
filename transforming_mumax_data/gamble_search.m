function [swchg, RK] = gamble_search(sstart, swinit,KRV,theta,Mdl_dtl, PZ, Bobj, particle_loc,control,MxB)
%GAMBLE_SEARCH Summary of this function goes here
%   Detailed explanation goes here


 tic

%% ------------------------------------------------------------------------
% Define the initial condition
pm = 4;
% Find where that sits in space (Pz)
pzcut =  find(MxB(pm,:) <= sstart, 1, 'first')-1;
e = 1e-16; % tolerance - numerical rounding   

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

[FWHMX(1),MLOC(1)] = FWHM(DNVC(1,:),theta);


%% ------------------------------------------------------------------------
% Now start to find next peak 
    
   ppp = 2;
   swchg(1) = sstart;
   swchg(2) = sstart/2;
   count = 1;
    
while abs(swchg(1) - swchg(2)) >= 1e-4
    
    if count == 1
        
    else
        tmps = swchg;
       
    
        if RK >= KRV 
            swchg(2) = (tmps(1)+tmps(2))/2;
        else 
            swchg(1) = tmps(2) - ((tmps(1)-tmps(2))/2);
        end 
    end 
    
for pull = 1:length(theta)
   %Bxnew = Bobj.BXx.*sin(theta) + Bobj.BZx.*cos(theta);
   %Bynew = Bobj.BXy.*sin(theta) + Bobj.BZy.*cos(theta);
   Bznew = Bobj(pm).BXz.*sin(theta(pull)) + Bobj(pm).BZz.*cos(theta(pull)); 
    
   lkatpln = Bznew(:,:,idx);
   BZM = (lkatpln >= swchg(2)) - (lkatpln <= -swchg(2));
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
  
       
    disp (['count = ', num2str(count),', range = ', num2str(swchg), ', RK = ', num2str(RK)])
    count = count +1;
end

toc


end

