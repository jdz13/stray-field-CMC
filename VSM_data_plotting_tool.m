
% clear

zz = fol_data_ext_function();

%%

figno = 71;


for tr = 1:size(zz,1)
    
    if isempty(zz(tr).data) == 1 
        
        continue
        
    else

    vhd = '.VHD';
    
     newSt = erase(zz(tr).name,vhd);
     newStr = strrep(newSt, '_', ' ');
    
figure(figno); clf; 
subplot (2,2,1); 
plot(zz(tr).data(:,8), zz(tr).data(:,10))
xlabel 'Applied field (Oe)'
ylabel 'Raw moment (X) (memu)'
title (['Raw plot - ', newStr]) 
subplot (2,2,2); 
plot(zz(tr).data(:,8), zz(tr).data(:,12))
xlabel 'Applied field (Oe)'
ylabel 'Raw moment (X) (memu)'
title 'Manipulated plot'
subplot (2,2,3); 
vsmplot(zz(tr).data(:,8), zz(tr).data(:,10))
xlabel 'Applied field (Oe)'
ylabel 'Raw moment (X) (memu)'
title 'Raw plot - slope'
subplot (2,2,4); 
vsmplot(zz(tr).data(:,8), zz(tr).data(:,12))
xlabel 'Applied field (Oe)'
ylabel 'Raw moment (X) (memu)'
title 'Manipulated plot - slope'

figno = figno+1;
clear newStr newSt
    end
end

clear tr