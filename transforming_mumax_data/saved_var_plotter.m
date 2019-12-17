function [] = saved_var_plotter (SWres,KRV,PM, RES,figno,con,iden)

plotter = zeros(length(KRV), length(PM), length(RES));

for i = 1:length(KRV)
    for j = 1:length(PM)
        for k = 1:length(RES)
               plotter(i,j,k) = nnz(SWres(i,:,j, k));
        end 
    end 
end


nox = 2;
if rem(length(RES),2) ==1
    noy = length(RES)+1/2;
else 
    noy = length(RES)/2;
end

figure(figno); clf;
    for jj = 1:length(RES)
        subplot(nox,noy,jj); imagesc(PM, KRV,plotter(:,:,jj));
        xlabel 'PM size [cm]'; ylabel 'KRV'; title (['Start field = ', num2str(RES(jj)),'T'])
        caxis([min(plotter,[],'all'),max(plotter,[],'all')]); colorbar
        ax = gca; trial = linspace(ax.YLim(1),ax.YLim(2),length(KRV)+1);
        trial = trial - (trial(2)-trial(1))/2;    trial(1) = [];
        yticks(trial);    yticklabels((fliplr(KRV)));
        
        if jj == ceil(noy/2)
            title (compose("Data for " + iden + "\nKRV  = [" + num2str(KRV) + "]\nPM sizes = [" + num2str(PM) + "], Con = " + num2str(con) + "\n \nStart field = " + num2str(RES(jj)) + "T"))
        end 
    end