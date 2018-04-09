function [S, VariableNames] = importdata_simplestruct_txt(filenames, Labels, ID)
%
% Script used to import the raw data from .fcs files. Before running the
% name of the file needs to be adjusted. Line 33 and 35. 
%
% Saves a i x 1 struct S with i being the number of measurements and with
% the following fields:
% Data      The raw data with m x n with m being the number of cells
% measured in that measurement and n being the number of surface proteins
% measured
% Labels    Zero, needs to be filled later!
% ID        Zero, needs to be filled later!
%
% Written by G.H. Tinnevelt at Radboud University at 12-6-2015
% Edited by G.H. Tinnevelt at Radboud University at 19-6-2015 (created
% function)
%
%%
% Get the folder with the filenames
folder = uigetdir;
files = dir([folder '\*.txt']);
delimiter = '\t'; %tab delimited text files
if size(files,1) == 0 
    files = dir([folder '\*.csv']); %comma seperated text files
    delimiter = ','; 
end

number_files = size(files,1);
if isempty(filenames)
    filenames = cell(size(files));
    for l1 = 1:number_files
        filenames(l1) = {[folder '\' files(l1).name]};
    end
end
%Initialise struct fields
Data = cell(number_files,1); %Empty data cell
if isempty(Labels)
    Labels = cell(number_files,1); %Empty Labels cell, needs to be filled later
end
if isempty(ID)
    ID = cell(number_files, 1); %Empty ID cell, needs to be filled later
end

% Import data and log transform the data.
for L1 = 1:number_files
    tmp.data = dlmread(filenames{L1},delimiter, 1,0); %load 1 datafile at a time
    Data(L1) = {tmp.data}; %concatenate datafile to matrix
end

%Save struct and VariableNames
S = struct('Data', Data, 'filenames', filenames, 'Labels', Labels, 'ID', ID);
[FileName, pathfile] = uiputfile({'*.mat'},'Save struct as:'); 
save([pathfile FileName ], 'S', '-v7.3') %Change Names 
fid = fopen(filenames{L1});
VariableNames = fgetl(fid); %load the first line of the data containing the VariableNames
[FileName2, pathfile2] = uiputfile({'*.mat'},'Save VariableNames as:'); 
save([pathfile2 FileName2 ], 'VariableNames', '-v7.3')