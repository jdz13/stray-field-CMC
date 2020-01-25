% RUN WITH 5.1 DATA (FULL WORLD).
clear

how_many = 4;
mag_IDs = [6,6,6,6]*1e-3;
mag_ODs = [10,20,30,40]*1e-3;% this should be the amount of different magnet sizes that you want to include in the code 
[particle_loc, control, Bobj, Mdl_dtl] = get_data_in (how_many, mag_IDs, mag_ODs);
pm_cl = mag_ODs;

tic

theta = linspace(0,pi/2,91);

KRV = [5,4,3,2.5,2];
PM = [1,2,3,4];
RES = [0.15,0.2,0.25,0.3,0.35,0.4];
con = 0.7;
plotyn = 0;

%Save outputs
[SaveVar25.varst, SaveVar25.SWres,SaveVar25.Bset,SaveVar25.FWHMres,SaveVar25.ind1res,SaveVar25.ind2res, SaveVar25.MxB] = search_tool_8_HWHM4x(KRV,PM,RES,pm_cl,Mdl_dtl,Bobj,theta,particle_loc,control,con);
SaveVar25.timer = toc; SaveVar25.comments = "Same as SV21 using all for PM size - using 5.1 data";

con = 0.9;
%Save outputs
[SaveVar26.varst, SaveVar26.SWres,SaveVar26.Bset,SaveVar26.FWHMres,SaveVar26.ind1res,SaveVar26.ind2res,SaveVar26.MxB] = search_tool_8_HWHM4x(KRV,PM,RES,pm_cl,Mdl_dtl,Bobj,theta,particle_loc,control,con);
SaveVar26.timer = toc;  SaveVar26.comments = "Same as SV22 using all PM size - using 5.1 data";