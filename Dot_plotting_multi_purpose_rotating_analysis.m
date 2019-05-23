option = 1;

if option == 1
    
percen_lim = 0.95; 

xmid = round(length(sample.px)/2);
ymid = round(length(sample.py)/2);

for np = 1:length(t)
    
    maxfield = max(max(abs(rotfield(np).fbz)));
    TESTMAT = abs(rotfield(np).fbz) == maxfield;
    [K,L] = find(TESTMAT);
         
    xline = abs(rotfield(np).fbz(K(1),:));
    yline = abs(rotfield(np).fbz(:,L(1)));
    
    fractionx = abs((rotfield(np).fbz(K(1),L(1))) * percen_lim);
    index1x = K(1);
    index2x = find(xline >= fractionx, 1, 'last');
    fwhmx = index2x-index1x + 1; % FWHM in indexes.
    % OR, if you have an x vector
    fwhmxx(np) = sample.px(index2x) - sample.px(index1x);
    
    
    fractiony = abs(rotfield(np).fbz(K(1),L(1)) * percen_lim);
    index1y = L(1);
    index2y = find(yline >= fractiony, 1, 'last');
    fwhmy = index2y-index1y + 1; % FWHM in indexes.
    % OR, if you have an x vector
    fwhmxy(np) = sample.py(index2y) - sample.py(index1y);
    
  end

figure(3)
plot(rad2deg(theta),fwhmxx, rad2deg(theta),fwhmxy)
legend ('x distance to reach 95% max','y distance to reach 95% max')
xlabel 'motion in degrees'; ylabel (['Distance to reach ' num2str(percen_lim*100) '% of max field'])

%% 
elseif option == 2


xlimp = sample.px(length(sample.px));
xlimn = sample.px(1);
ylimp = sample.py(length(sample.py));
ylimn = sample.py(1);

for np = 1:7:length(theta)
    figure(110)
    subplot(2,2,1)
    imagesc(sample.px,sample.py,rotfield(np).fbx)
    axis equal; axis([xlimn,xlimp,ylimn,ylimp]);
    colorbar ; caxis  ([-0.4*Bmax, 0.4*Bmax])
    xlabel 'length(m)', ylabel 'length(m)'; title 'real coord''s (XYZ) - X'
    subplot(2,2,2)
    imagesc(sample.px,sample.py,rotfield(np).fby)
    axis equal; axis([xlimn,xlimp,ylimn,ylimp]);
    colorbar ; caxis  ([-0.1*Bmax, 0.1*Bmax])
    xlabel 'length(m)', ylabel 'length(m)'; title 'Y'
    subplot(2,2,3)
    imagesc(sample.px,sample.py,rotfield(np).fbz)
    axis equal; axis([xlimn,xlimp,ylimn,ylimp]);
    colorbar ; caxis  ([-0.8*Bmax, 0.8*Bmax])
    xlabel 'length(m)', ylabel 'length(m)'; title 'Z'
    subplot(2,2,4)
    imagesc(sample.px,sample.py,rotfield(np).fbtot)
    axis equal; axis([xlimn,xlimp,ylimn,ylimp]);
    colorbar ; caxis  ([-0.8*Bmax, 0.8*Bmax])
    xlabel 'length(m)', ylabel 'length(m)'; title 'Total field'
    
    
    
    %-------------------------------------------------------------------
    
    figure(111)
    subplot(2,2,1)
    imagesc(sample.px,sample.py,Akoun(np).HxAkoun)
    axis equal; axis([xlimn,xlimp,ylimn,ylimp]);
    colorbar ; caxis  ([-0.4*Bmax, 0.4*Bmax])
    xlabel 'length(m)', ylabel 'length(m)'; title 'prime coord''s (X''Y''Z'') - X'''
    subplot(2,2,2)
    imagesc(sample.px,sample.py,Akoun(np).HyAkoun)
    axis equal; axis([xlimn,xlimp,ylimn,ylimp]);
    colorbar ; caxis  ([-0.1*Bmax, 0.1*Bmax])
    xlabel 'length(m)', ylabel 'length(m)'; title 'Y''' 
    subplot(2,2,3)
    imagesc(sample.px,sample.py,Akoun(np).HzAkoun)
    axis equal; axis([xlimn,xlimp,ylimn,ylimp]);
    colorbar ; caxis  ([-0.8*Bmax, 0.8*Bmax])
    xlabel 'length(m)', ylabel 'length(m)'; title 'Z'''
    subplot(2,2,4)
    imagesc(sample.px,sample.py,Akoun(np).modBAkoun)
    axis equal; axis([xlimn,xlimp,ylimn,ylimp]);
    colorbar ; caxis  ([-0.8*Bmax, 0.8*Bmax])
    xlabel 'length(m)', ylabel 'length(m)'; title 'Total field'
    
   
    
    
    
end

%% ----------------------------------------------------------------------
elseif option == 3
    
    for np = 1:length(theta)

        mmaxf(np) = max(max(rotfield(np).fbz));
        mminf(np) = min(min(rotfield(np).fbz));
        mmaxxx(np) = max(max(abs(rotfield(np).fbz)));
        
    end

    thetad = rad2deg(theta);
    figure(82)
    hold on
    plot(thetad,mmaxf)
    plot(thetad,mminf)
    plot(thetad,mmaxxx)
    hline(0,'k','')
    legend ('Max','Min','Max(abs)')    
    xlabel 'Angle (degrees)'; ylabel 'Field Max/Min (T)'
    title 'Maximum and Minimum fields for each plane'
    

end
    


