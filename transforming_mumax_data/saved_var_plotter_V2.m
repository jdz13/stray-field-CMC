function [] = saved_var_plotter_V2 (ist,figno)

plotter = zeros(length(ist.varst.KRV), length(ist.varst.PM), length(ist.varst.RES));
str = inputname(1);

for i = 1:length(ist.varst.KRV)
    for j = 1:length(ist.varst.PM)
        for k = 1:length(ist.varst.RES)
               plotter(i,j,k) = nnz(ist.SWres(i,:,j, k));
        end 
    end 
end

nox = 2;
if rem(length(ist.varst.RES),2) ==1
    noy = length(ist.varst.RES)+1/2;
else 
    noy = length(ist.varst.RES)/2;
end

figure(figno); clf;
    for jj = 1:length(ist.varst.RES)
        subplot(nox,noy,jj); imagesc(ist.varst.PM, ist.varst.KRV,plotter(:,:,jj));
        xlabel 'PM size [cm]'; ylabel 'KRV'; title (['Start field = ', num2str(ist.varst.RES(jj)),'T'])
        caxis([min(plotter,[],'all'),max(plotter,[],'all')]); colorbar
        ax = gca; trial = linspace(ax.YLim(1),ax.YLim(2),length(ist.varst.KRV)+1);
        trial = trial - (trial(2)-trial(1))/2;    trial(1) = [];
        yticks(trial);    yticklabels((fliplr(ist.varst.KRV)));
        
        if jj == ceil(noy/2)
            title (compose("Data for " + str + "\nKRV  = [" + num2str(ist.varst.KRV) + "]\nPM sizes = [" + num2str(ist.varst.PM) + "], Con = " + num2str(ist.varst.CON) + "\n \nStart field = " + num2str(ist.varst.RES(jj)) + "T"))
        end 
    end