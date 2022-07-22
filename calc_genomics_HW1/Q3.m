%% Q3 Genomics
clear all; clc; close all;

%% Subsections 1-2

% uploading data and finding YAR068W coordinates
fid = fopen('array.txt');
data = textscan(fid,'%s%c%d%d%s%s%c%c','HeaderLines',3);
fclose(fid);
for i=1:size(data,2)
    for j=1:length(data{:,1})
        if strcmp(data{1,i}(j,1),'YAR068W')
            YAR068W_coordinates=[j,i];
        end
    end
end
%% Subsection 3

% finding coordinates column in the cell array
clear i;clear j;
for i=1:size(data,2)
    if iscell(data{1,i})
        if contains(cell2mat(data{1,i}(1,1)),'..')
            coordinates_column=i;
        end
    end
end
%% Subsection 4

% finding YAL002W ending coordinates
clear i;
for i=1:size(data,2)
    for j=1:length(data{:,1})
        if strcmp(data{1,i}(j,1),'YAL002W')
            YAL002W_coordinates=[j,i];
        end
    end
end

% taking coordinates before and after '..'
boundry_index=strfind(data{1,coordinates_column}...
    {YAL002W_coordinates(1),1},'..');
YAL002W_ending_coordinate=str2double(data{1,coordinates_column}...
    {YAL002W_coordinates(1),1}(boundry_index+2:end));

%% Subsection 5

% creating a seperated array for starting and ending coordinates
clear i; starting_coordinates=[]; ending_coordinates=[];
for i=1:length(data{:,coordinates_column})
    boundry_index=strfind(data{1,coordinates_column}{i,1},'..');
    starting_coordinates(end+1)=str2double(data{1,coordinates_column}...
        {i,1}(1:boundry_index-1));
    ending_coordinates(end+1)=str2double(data{1,coordinates_column}...
        {i,1}(boundry_index+2:end));
end
coordinates_array=cat(1,starting_coordinates,ending_coordinates);
coordinates_array=coordinates_array';