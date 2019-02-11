[Sact,pzact, maxact] = mumax_data_compare(); % acutal magnet
[Srodv,pzrodv, maxrodv] = mumax_data_compare(); % rod constant volume
[Srodd,pzrodd, maxrodd] = mumax_data_compare(); % rod constant dimension
[Scubd,pzcubd, maxcubd] = mumax_data_compare(); % cube constant dimension
[Scubv,pzcubv, maxcubv] = mumax_data_compare(); % cube constant volume

figure(135)
subplot(1,2,1)
plot(pzact,maxact,pzrodv,maxrodv,pzcubv,maxcubv)
xlabel 'Distance from magnet centre (m)'; ylabel 'Field (T)'
title 'Comparison of field characteristics for constant volume shapes'
legend ('Actual magnet' ,'equivalent volume rod' ,'equivalent volume cube')

subplot(1,2,2)
semilogy(pzact,maxact,pzrodv,maxrodv,pzcubv,maxcubv)
xlabel 'Distance from magnet centre (m)'; ylabel 'Field (T)'
title 'Log scale'
legend ({['Actual magnet' newline '20mm OD 6mm ID 20mm L'],['equivalent volume rod' newline '19mm OD, 20mm L'],['equivalent volume cube' newline '18mm L']})

%%

figure(136)
subplot(1,2,1)
plot(pzact,maxact,pzrodd,maxrodd,pzcubd,maxcubd)
title 'Comparison of field characteristics for constant dimension shapes'
xlabel 'Distance from magnet centre (m)'; ylabel 'Field (T)'
legend ('Actual magnet' ,'equivalent dim rod' ,'equivalent dim cube')
 
subplot(1,2,2)
semilogy(pzact,maxact,pzrodd,maxrodd,pzcubd,maxcubd)
xlabel 'Distance from magnet centre (m)'; ylabel 'Field (T)'
title 'Log scale'
xlabel 'Distance from magnet centre (m)'; ylabel 'Field (T)'
title 'Comparison of field characteristics for constant volume shapes'
legend ({['Actual magnet' newline '20mm OD 6mm ID 20mm L'],['equivalent dim rod' newline '20mm OD, 20mm L'],['equivalent dim cube' newline '20mm L']})
