function [zz] = fol_data_ext_function()
%FOL_DATA_EXT_FUNCTION Summary of this function goes here
%   Detailed explanation goes here


path = uigetdir; % Let the user select the folder that they want to convert.
oldFolder = cd (path); % Change directory and remember the old one to go back to.

tic

zz = dir(path);% This finds all of the files in this folder.

folName = 'Output'; % The output folder name to be placed in the input folder.

% Make the folder.
if ~exist(folName, 'dir')
  mkdir(folName);
end

outpath = [path, '\', folName]; % Make sure that the files are put into the new folder.

for p = 3:length(zz) % 3 as first two are always directories.

    
    
    if zz(p).isdir == 1   % Means that folders and directories are not looked at.
    
    elseif contains(zz(p).name,'.png')==1 % Does not look at the png images.
        
    else
        
    
        
        
filename = [path,'\',zz(p).name];

fid = fopen(filename,'r');
line=fgetl(fid); % fgetl(fileID) returns the next line of the specified file, removing the newline characters.

    clear data
    clear obj
    clear header
    
% This is for skiping lines until we reach the column names.
pattern = 'Time_since_start'; % Tell it what to look for.
while ~contains(line,pattern)==1
    line=fgetl(fid);
end

% This will put the column names into a cell. 
data.columns={};
remain=line;
i=1;
substr = {};

while isempty(remain)==0
    [str, remain] = strtok(remain); % token = strtok(str) parses str from left to right, using whitespace characters as delimiters, and returns part or all of the text in token.
    
    % If statment to ensure column names with spaces are put into one
    % data.column entry.
    if isempty(substr) 
        data.columns{i} = str();
        i=i+1;
    end
    
end

% Because there is a random space at the end colomn.
data.columns(length(data.columns)) = [];

% Add the headers to the object.
obj.headers = data.columns;

% Move down three lines to reach the numerical data.

line=fgetl(fid);line=fgetl(fid);line=fgetl(fid); % Apparently needed four? EDIT. Nope - three. One needed below these (added). 

% Now we have reached the numerical data.
data_line = sscanf(fgetl(fid),'%f'); % A = sscanf(str,formatSpec) reads data from str, converts it according to the format specified by formatSpec, and returns the results in an array.
data.values = zeros(1,length(data.columns));% Preload the array for getting the numerical data.

% This will scan line by line and put it into data.values.
i=1;
while isempty(data_line) ~= 1
    data.values(i,:)=data_line'; % Put the line in the right place.
    
    data_line=sscanf(fgetl(fid),'%f'); % Move to the next line.
    i=i+1; % Move to the next row of data.values.
end

line=fgetl(fid); % new line needed here. Otherwise doesn't pick up the data. 

% Need to get the full loop information - if there is more.

pattern = 'New Section:';
% This is for skiping lines until we reach the column names.
while ~contains(line,pattern)==1
    line=fgetl(fid);
end

data_line = sscanf(fgetl(fid),'%f');
% This will scan line by line and put it into data.values.
while isempty(data_line) ~= 1
    data.values(i,:)=data_line';
    
    data_line=sscanf(fgetl(fid),'%f');
    i=i+1;
end

% Add the data to the object.
obj.data = data.values;

fclose(fid); % Close the file in Matlab.

vhd = '.VHD';
vrd = '.VRD';
prompt = 'What is the file extention? (Case sensitive and please include dot)';

if contains(zz(p).name,vhd)==1
    newStr = erase(zz(p).name,vhd);
    
elseif contains(zz(p).name,vrd)==1
    newStr = erase(zz(p).name,vrd);

else 
    ext = inputdlg(zz(p).name,'input');
    newStr = erase(zz(p).name,ext{1});
    
end

fileout = [outpath,'\', newStr, '.txt']; % Add the .txt extention to make a txt file.
outputsize =  size(data.values,1); 

% When you output data you have to tell it explicitly how you want it
% outputted. As all of our titles are strings, and our data is numeric, we
% can set up a loop which will format the data outputs in this way,
% regardless of their number. 
% Tab delimiter set below - can be changed if needed. 
tiled ='%s'; 
tiled2 = '%f';

% Outputing data with header to outputpath.
    for q = 1:length(data.columns)
        header(q) = data.columns(q);
    end
    
    for q = 1:length(data.columns)-1
        tiled = [tiled, ' \t %s']; % If tab delimiter is not wanted, change here.
        tiled2 = [tiled2, ' \t %f']; % And here.
    end
    
    tiled = [tiled,' \n']; % Finish by defining a new line, to ensure data is in columns.
    tiled2 = [tiled2, ' \n'];
    
    % Export the data out.
    fid = fopen(fileout,'w'); % Open the output file. 
    fprintf(fid,tiled,header{1,:}); % Put in the headers.
    % Output row by row. 
    for j = 1:outputsize
        fprintf(fid,tiled2,data.values(j,:));
    end
    fclose('all'); % Close all open files.
    
    zz(p).data = data.values;
   
    end
    
end
cd (oldFolder) % Change back to the old directory. As not to confuse everyone. 
toc


end

