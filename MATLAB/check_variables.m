%% Check if there are the same variables for all the sample
% Function created to check if, for all the samples in the folder, the
% variables are the same and in the same order

% Output:
% Table with name of the variables measured per sample. 

% Written by R.Folcarelli at Radboud University

function [Var]=check_variables(filenames)

folder = uigetdir;
files = dir([folder '/*.fcs']);
if size(files,1) == 0 
    files = dir([folder '/*.LMD']);
end
number_files = size(files,1);
if isempty(filenames)
    filenames = cell(size(files));
    for l1 = 1:number_files
        filenames(l1) = {[folder '/' files(l1).name]};
    end
end

Var=cell(number_files,1);


for L1 = 1:number_files
    tmp.data = readfcs(filenames{L1})'; %load 1 datafile at a time
    [~, var, ~,~] = readfcs(filenames{L1});
    for L2 = 1:length(var)
        Var{L1,L2}=var(L2);
    end
end

return 