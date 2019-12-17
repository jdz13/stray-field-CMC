timer12 = [0,0,0,0,0];
plotyn = 0;

%-------------------------------------------------------------------------
% set up the 4x model parameters and see if it runs nicely. /how fast too 
tic

theta = linspace(0,pi/2,91);

f = [11,12];

KRV = [5,4,3,2.5,2];
PM = [1,2,3,4];
RES = [0.15,0.2,0.25,0.3,0.35,0.4];%,0.45,0.5,0.55,0.6]

%Save outputs
[varst, SaveVar7_SWres,SaveVar7_Bset,SaveVar7_FWHMres] = search_tool_V6_4x(KRV,PM,RES,f,pm_cl, PZ,Mdl_dtl,Bobj,theta,particle_loc,control,plotyn);

timer12(1) = toc;
%% -------------------------------------------------------------------------
tic
KRV = [5,4,3,2.5,2];
PM = [2,3,4];
RES = [0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6];

f = f+2;

[varst, SaveVar8_SWres,SaveVar8_Bset,SaveVar8_FWHMres] = search_tool_V6_4x(KRV,PM,RES,f,pm_cl, PZ,Mdl_dtl,Bobj,theta,particle_loc,control,plotyn);
timer12(2) = toc;

%%

tic

theta = linspace(0,pi/2,91);

KRV = [5,4,3,2.5,2];
PM = [1,2,3,4];
RES = [0.15,0.2,0.25,0.3,0.35,0.4];
con = 0.9;

%Save outputs
[varst, SaveVar9_SWres,SaveVar9_Bset,SaveVar9_FWHMres,SaveVar9_ind1res,SaveVar9_ind2res] = search_tool_V7_HWHM4x(KRV,PM,RES,pm_cl, PZ,Mdl_dtl,Bobj,theta,particle_loc,control,con);
timer12(3) = toc;
%%
tic
con = 0.7;
%Save outputs
[varst, SaveVar10_SWres,SaveVar10_Bset,SaveVar10_FWHMres,SaveVar10_ind1res,SaveVar10_ind2res] = search_tool_V7_HWHM4x(KRV,PM,RES,pm_cl, PZ,Mdl_dtl,Bobj,theta,particle_loc,control,con);
timer12(4) = toc;




%%

% RUN WITH 5.1 DATA (FULL WORLD).

%clear
%run(get_data_in)
tic

theta = linspace(0,pi/2,91);

KRV = [5,4,3,2.5,2];
PM = [2,3,4];
RES = [0.15,0.2,0.25,0.3,0.35,0.4];
con = 0.7;

%Save outputs
[SaveVar11b.varst, SaveVar11b.SWres,SaveVar11b.Bset,SaveVar11b.FWHMres,SaveVar11b.ind1res,SaveVar11b.ind2res] = search_tool_V7_HWHM4x(KRV,PM,RES,pm_cl, PZ,Mdl_dtl,Bobj,theta,particle_loc,control,con);
SaveVar11b.timer = toc; SaveVar11b.comments = "same as SV11 but this time with 0.7 condition - reloaded 5.1 data";

tic
[SaveVar12.varst, SaveVar12.SWres,SaveVar12.Bset,SaveVar12.FWHMres] = search_tool_V6_4x(KRV,PM,RES,pm_cl, PZ,Mdl_dtl,Bobj,theta,particle_loc,control);
SaveVar12.timer = toc; SaveVar12.comments = "Same as SV11 though this time using FWHM, not HWHM";

%%

% ADD THE VARIABLE TO THE SAVED ONES BEFORE RUNNING THIS. 

% RUN WITH 5.2 DATA (SHORTENED WORLD).
%clear
%run(get_data_in)
tic

theta = linspace(0,pi/2,91);

KRV = [5,4,3,2.5,2];
PM = [2,3,4];
RES = [0.15,0.2,0.25,0.3,0.35,0.4];
con = 0.9;
plotyn = 0;

%Save outputs
[SaveVar13.varst, SaveVar13.SWres,SaveVar13.Bset,SaveVar13.FWHMres,SaveVar13.ind1res,SaveVar13.ind2res] = search_tool_V7_HWHM4x(KRV,PM,RES,pm_cl, PZ,Mdl_dtl,Bobj,theta,particle_loc,control,con);
SaveVar13.timer = toc;
%%
tic
[SaveVar14.varst, SaveVar14.SWres,SaveVar14.Bset,SaveVar14.FWHMres] = search_tool_V6_4x(KRV,PM,RES,f,pm_cl, PZ,Mdl_dtl,Bobj,theta,particle_loc,control,plotyn);
SaveVar14.timer = toc;

%%

% RUN WITH 5.1 DATA (FULL WORLD).

%clear
%run(get_data_in)
tic

theta = linspace(0,pi/2,91);

KRV = [5,4,3,2.5,2];
PM = [2,3,4];
RES = [0.15,0.2,0.25,0.3,0.35,0.4];
con = 0.7;

%Save outputs
[SaveVar15.varst, SaveVar15.SWres,SaveVar15.Bset,SaveVar15.FWHMres,SaveVar15.ind1res,SaveVar15.ind2res] = search_tool_8_HWHM4x(KRV,PM,RES,pm_cl,Mdl_dtl,Bobj,theta,particle_loc,control,con);
SaveVar15.timer = toc; SaveVar15.comments = "same as SV11 but this time with search tool V* - changed PZ/MxB creation - reloaded 5.1 data";
