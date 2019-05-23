for iio = [1,2,3,6]
fig = figure(iio);
print(fig,'-Phpljm477-304b (HP Color LaserJet MFP M477fdw)', '-bestfit','-dwinc')
end

for sl = 1:length(pm_cl)
stravlin = ['C:\Users\Jake\Documents\MATLAB\Git\stray-field-CMC-cloned\Output_plots\Mag_size_param_sweep\2nd run\key ratio for sample mag surface dist = ' num2str(pm_cl(sl)*1000) '_mm.fig'];
fig = figure(sl+74);
set(gcf, 'Position', get(0, 'Screensize')); % Makes the figure full screen
saveas(fig,stravlin);
end