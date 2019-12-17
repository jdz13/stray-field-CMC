function [] = VSM_stab_plotting(t,sig,Tc,Ta,H,figno)

    time = t - t(1);
    timeh = time/3600;
    HkOe = H/1000;
    
    sd.sig = std(sig); sd.Tc = std(Tc); sd.Ta = std(Ta); sd.H = std(H);
    pd.sig = (max(sig) - min(sig))/max(sig)*100; pd.Tc = (max(Tc) - min(Tc))/max(Tc)*100;
    pd.Ta = (max(Ta) - min(Ta))/max(Ta) *100; pd.H = (max(H) - min(H))/max(H)*100;
    
    figure (figno); 
    subplot(4,1,1)
    plot(timeh, sig); xlabel 'time [hrs]'; ylabel 'LIA X signal [V]'
    legend(compose("stdev  = " + num2str(sd.sig) + " \n% change = " + num2str(pd.sig) + "%"))
    subplot(4,1,2)
    plot(timeh, Tc); xlabel 'time [hrs]'; ylabel 'Control temperature [degC]'
    legend(compose("stdev  = " + num2str(sd.Tc) + " \n% change = " + num2str(pd.Tc) + "%"))
    subplot(4,1,3)
    plot(timeh, Ta); xlabel 'time [hrs]'; ylabel 'Ambient temperature [degC]'
    legend(compose("stdev  = " + num2str(sd.Ta) + " \n% change = " + num2str(pd.Ta) + "%"))
    subplot(4,1,4)
    plot(timeh, HkOe); xlabel 'time [hrs]'; ylabel 'Applied field [kOe]'
    legend(compose("stdev  = " + num2str(sd.H) + " \n% change = " + num2str(pd.H) + "%"))

end