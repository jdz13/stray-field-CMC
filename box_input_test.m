clear

prompt = {'Number of cells (space separated):','World size (space separated):','MSat','M direction (space separated 3 Vector)'};
input_title = 'Input';
dims = [1 35];
definput = {'50 50 50','1e-3 1e-3 1e-3','1e6','0 0 1'};
answer = inputdlg(prompt,input_title,dims,definput);


grid_size = str2num(cell2mat(answer(1)));
world_range = str2num(cell2mat(answer(2)));
Msat = str2num(cell2mat(answer(3)));
unimagV = str2num(cell2mat(answer(4)));

clear prompt input_title dims definput answer

[space, Akoun] = CMC(grid_size, world_range, unimagV, Msat);

figure(2)
slice(space(1).X,space(1).Y,space(1).Z, Akoun(1).totB, 0,0,0)
caxis([-0.001,0.001])
colorbar

%%
[Bx, By, Bz,range] = Mumax_data_tool();

Checks_on_code_rough
