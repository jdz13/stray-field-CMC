

% what to test? 
% 1 = KRV
% 2 = pm
% 3 = start field 

PMno = 3;
SFno = 4;
KRVno = 2;

plotter12 = zeros(1,100);

testop = 3; 

if testop == 1

    for i = 1:length(KRV)
        plotter12(i) = nnz(SWres(i,:,PMno, SFno));
    end
    
elseif testop == 2

    for i = 1:length(PM)
        plotter12(i) = nnz(SWres(KRVno,:,i, SFno));
    end
    
        
elseif testop == 3
    
    for i = 1:length(RES)
        plotter12(i) = nnz(SWres(KRVno,:,PMno, i));
    end

end 

plotter12(plotter12==0)=nan; 



%%


plotter19 = zeros(length(KRV), length(PM),length(RES));

    for i = 1:length(KRV)
        for j = 1:length(PM)
            for k = 1:length(RES)
                   plotter19(i,j,k) = nnz(SaveVar1_SWres(i,:,j, k));
            end 
        end 
    end

    
%%
    
figure(57); clf; subplot(2,2,1)
plot(VSM_t1,VSM_M1*1e6); title 'zero field - with probe'
ylabel 'Moment (\muemu)' ; xlabel 'time (s)'
subplot(2,2,2)
plot(VSM_t1,VSM_T1)
ylabel 'Temperature (C)' ; xlabel 'time (s)'
subplot(2,2,3)
plot(VSM_t2,VSM_M2*1e6); title '1.5T field - with probe'
ylabel 'Moment (\muemu)' ; xlabel 'time (s)'
subplot(2,2,4)
plot(VSM_t2,VSM_T2)
ylabel 'Temperature (C)' ; xlabel 'time (s)'


figure(58); clf; subplot(2,2,1)
plot(VSM_t3,VSM_M3*1e6); title 'zero field - without probe'
ylabel 'Moment (\muemu)' ; xlabel 'time (s)'
subplot(2,2,2)
plot(VSM_t3,VSM_T3)
ylabel 'Temperature (C)' ; xlabel 'time (s)'
subplot(2,2,3)
plot(VSM_t4,VSM_M4*1e6); title '1.5T field - without probe'
ylabel 'Moment (\muemu)' ; xlabel 'time (s)'
subplot(2,2,4)
plot(VSM_t4,VSM_T4)
ylabel 'Temperature (C)' ; xlabel 'time (s)'


figure(59); clf; subplot(2,2,1)
plot(VSM_t5,VSM_M5*1e6); title 'zero field - no vibration'
ylabel 'Moment (\muemu)' ; xlabel 'time (s)'
subplot(2,2,2)
plot(VSM_t5,VSM_T5)
ylabel 'Temperature (C)' ; xlabel 'time (s)'
subplot(2,2,3)
plot(VSM_t6,VSM_M6*1e6); title '1.5T field - no vibration'
ylabel 'Moment (\muemu)' ; xlabel 'time (s)'
subplot(2,2,4)
plot(VSM_t6,VSM_T6)
ylabel 'Temperature (C)' ; xlabel 'time (s)'


figure(59); clf; subplot(2,2,1)
plot(VSM_t6./3600,VSM_M6*1e6); title 'zero field - no vibration'
ylabel 'Moment [\muemu]' ; xlabel 'time [Hrs]'
subplot(2,2,2)
plot(VSM_t6./3600,VSM_T6)
ylabel 'Temperature [degC]' ; xlabel 'time [Hrs]'
subplot(2,2,[3,4])
plot(VSM_T6,VSM_M6*1e6); 
ylabel 'Moment [\muemu]' ; xlabel 'Temperature [degC]'

%%



    res = RES(4);

    pm = PM(1);
    
for  count2 = 1:size(KRV,2)
    
% Define the initial conditions
SH0 = 0.2500;
swinit = res; % What's the max channel value?
SWnext = 0; % initial condition
count = 1;

while abs(swinit - swnext) > 1e-4 && SH0 > MxB(pm,length(MxB(pm,:)))

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

[FWHMX(1),MLOC(1)] = FWHM(DNVC(1,:),theta);
%RKout = zeros(1,100); 

SWnextpos = MLOC(1)+(KRV(count2)*FWHMX(1));  

SWnext = swinit*cos(SWnextpos);

SWres(count2,count+1,pmcount,rescount) = SWnext;
FWHMres(count2,count+1,pmcount,rescount) = FWHMX(1); % if +1 not needed then can use FWHMres(:,1,:,:) = [];

% manipulate the results to run the next leg

tmps(1) = swinit; tmps(2) = SWnext;

SH0 = (tmps(1)+tmps(2))/2;
swinit = tmps(2); 

count = count +1;

disp (['count = ', num2str(count),', range = ', num2str(tmps)])

end
end

%%


figure(10)
subplot(2,1,1)
plot(thetad,NVC(1,:))
xlabel 'Angle (deg)'
ylabel 'Normalised volume component'
%xlim([0,20])
%ylim([-1.1,1.1])

subplot(2,1,2)
plot(thetad(2:length(thetad)),DNVC(1,:))
xlabel 'Angle (deg)'
ylabel 'Differentiated signal'
%xlim([0,20])
%ylim([-0.17,0.17])