prompt = {'Number of cells - X', 'Number of cells - Y', 'Number of cells - Z'};
definput = {'30','30','30'};
uinput = (inputdlg(prompt,'input',[1 30],definput));
cell_size = zeros(1,length(definput));

for i = 1:length(definput)
    cell_size(i) = str2double(uinput(i));
end

