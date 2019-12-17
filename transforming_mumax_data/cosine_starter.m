
MxB = zeros(length(pm_cl),size(PZ,2));

for pm = 1:length(pm_cl)
    for ww = 1:size(PZ,2)
    
        e = 1e-16;
        idx = find(abs(Mdl_dtl(pm).topmagLinez - PZ(ww))<=e);
        lkpln = Bobj(pm).BZz(:,:,idx);
        MxB(pm,ww) = max(max(lkpln));
        
    end
end
%% ------------------------------------------------------------------------

figure(1) 
semilogy(PZ,MxB)
xlabel 'Magnet surface to sample distance [m]'
ylabel 'Maximum field value [T]'

%% ------------------------------------------------------------------------

 figno = 31; 
    nm =75; % z height location in array;
    sp = 1; % sample space in array 
    pm = 4; % magnet size in array
    
   Cosvar = cos(theta);
   Bcos = MxB(pm, nm).*Cosvar;
    
    figure(figno); clf; subplot(2,1,1);hold on; subplot(2,1,2);hold on;
    for cx = 1:size(swfield,2)
    subplot(2,1,1); plot(Bcos, norm_vol_comp(nm,sp,pm).data(:,cx));
    subplot(2,1,2); plot(Bcos(2:length(Bcos)), diphis(nm,sp,pm).data(:,cx));
    end
    legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
    xlabel 'Field [T]'; ylabel 'Differentiated N particles switched'
 
    subplot(2,1,1)
    legendCell = cellstr(num2str(swfield', 'SF =%-d (T)')); legend(legendCell)
    xlabel 'Field [T]'; ylabel 'Normalised number of particles switched'
  
    clear cx nm figno sp pm