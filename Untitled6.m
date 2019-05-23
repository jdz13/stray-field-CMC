clear
OOP =  fol_data_ext_function();
IP =  fol_data_ext_function();

ipc = [24,8,12,20,16]; oopc = [4,20,8,16,12];
titvar={'No sample','Clean sample holder','PTFE','Tape I VSMok','Tape II'};

for jol = 1:length(ipc)
    
    figno = jol+6;    figure(figno);    subplot(2,2,1); hold on
    plot(IP(ipc(jol)).data(:,5),IP(ipc(jol)).data(:,10), ...
        OOP(oopc(jol)).data(:,5),OOP(oopc(jol)).data(:,10))
    xlabel 'Applied Field (Oe)'; ylabel 'Raw Signal (memu)'
    title (sprintf('Raw signal for both IP and OOP settings - %s',...
        string(titvar(jol))))
    legend ('IP settings','OOP settings'); subplot(2,2,2); hold on
    plot(IP(ipc(jol)).data(:,5),IP(ipc(jol)).data(:,12),...
        OOP(oopc(jol)).data(:,5),OOP(oopc(jol)).data(:,12))
    xlabel 'Applied Field (Oe)'; ylabel 'Moment (emu)'a
    title 'Processed signal for both IP and OOP settings'
    legend ('IP settings','OOP settings'); subplot(2,2,3); hold on
    vsmplot(IP(ipc(jol)).data(:,5),IP(ipc(jol)).data(:,10)); 
    vsmplot( OOP(oopc(jol)).data(:,5),OOP(oopc(jol)).data(:,10))
    xlabel 'Applied Field (Oe)'; ylabel 'Raw Signal (memu)'
    title 'Corrected raw signal for both IP and OOP settings '
    legend ('IP settings','OOP settings'); subplot(2,2,4); hold on
    vsmplot(IP(ipc(jol)).data(:,5),IP(ipc(jol)).data(:,12)); 
    vsmplot(OOP(oopc(jol)).data(:,5),OOP(oopc(jol)).data(:,12))
    xlabel 'Applied Field (Oe)'; ylabel 'Moment (emu)'
    title 'Corrected processed signal for both IP and OOP settings '
    legend ('IP settings','OOP settings'); subplot(2,2,2); hold on
end


for sl = 1:length(ipc)
stravlin = join(['C:\Users\Jake\Documents\MATLAB\Git\stray-field-CMC-cloned\Output_plots\VSM_runs_april_2019\' string(titvar(sl)) '.bmp']);
fig = figure(sl+6);
set(gcf, 'Position', get(0, 'Screensize')); % Makes the figure full screen
saveas(fig,stravlin);
end