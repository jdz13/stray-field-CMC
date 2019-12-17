 tic

%% ------------------------------------------------------------------------

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

% Define the initial conditions
swinit = 0.2500; % What's the max channel value?
KRV = 5;
pm = 4;


% Find where that sits in space (Pz)
pzcut =  find(MxB(pm,:) <= swinit, 1, 'first')-1;
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

[FWHMX(1),MLOC(1)] = FWHM(DNVC(1,:),theta);
RKout = zeros(1,100); 

%% ------------------------------------------------------------------------
% Now start to find next peak 

   swchg(1) = swinit;
   swchg(2) = swinit/2;
   count = 1;
   ppp = 2;
     
while abs(swchg(1) - swchg(2)) >= 1e-4
    
    if count == 1
        
    else
        tmps = swchg;
       
    
        if RK - KRV >= -e
            swchg(2) = (tmps(1)+tmps(2))/2;
        else 
            swchg(2) = tmps(2) - ((tmps(1)-tmps(2))/2);
            swchg(1) = tmps(2);
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

    RKout(count) = RK; 
    count = count +1;
end



toc
