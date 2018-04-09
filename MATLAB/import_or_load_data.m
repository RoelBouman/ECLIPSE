function [S, VariableNames, load_modes] = import_or_load_data(load_modes)
%%
% Function created to either import or load the data prior to MFC data
% analysis.
%
% Input: load_modes (optional) 2x1 cell
% if the load_modes is empty it will asks for questions
% steps:
% load_mode{1} = 'import' or 'load' the data
% load_mode{2} = 0 or 1 (Excel file or not)
%
% Output:
% struct S with i being the number of measurements and with
% the following fields:
% Data      The raw data with m x n with m being the number of cells
% measured in that measurement and n being the number of surface proteins
% measured
% Labels
% ID

% VariableNames (nx1 cell with the variable names)
% load_modes (the options choosen)

% Written by G.H. Tinnevelt at Radboud University at 19-6-2015


%%
%% Question whether to load or import the data
if isempty(load_modes);
    check_input_load_data = 0;
    while check_input_load_data == 0
        load_data = input('Do you need to import the data or load the data from struct (import/load)? ', 's');
        if strcmp(load_data, 'import') || strcmp(load_data, 'load')
            check_input_load_data = 1;
        else
            warning('Wrong input for load_data. Either put import or load.')
        end
    end
else
    load_data = load_modes{1};
end
%% Import data
if strcmp(load_data, 'import')
    data_format = input('What is the extension of the files? \n 1: .FCS \n 2: .LMD \n 3: .TXT \n 4: .CSV \n Extension: ');
    if isempty(load_modes)
        xls_file = input('Do you have a Excel sheet with the information on filenames, labels and ID? \n 0 for no sheet \n 1 for a sheet \n Excel sheet: ');
    else
        xls_file = load_modes{2};
    end
    if  xls_file == 1
        [FileName,PathName] = uigetfile({'*.xls'; '*.xlsx'}, 'Select Excel file with information on filenames, labels and ID');
        xls_file = [PathName FileName];
        [~, ~, xls_sheet] = xlsread(xls_file);
        if data_format == 1 || data_format == 2
            [S, VariableNames] = importdata_simplestruct_fcs(xls_sheet(:,1), xls_sheet(:,2), xls_sheet(:,3));
        else
            [S, VariableNames] = importdata_simplestruct_txt(xls_sheet(:,1), xls_sheet(:,2), xls_sheet(:,3));
        end
    else %No excel file
        if data_format == 1 || data_format == 2
            [S, VariableNames] = importdata_simplestruct_fcs([], [], []);
        else
            [S, VariableNames] = importdata_simplestruct_txt([], [], []);
        end
        warning('Please fill in the Labels and ID and save the struct again!')
    end
    load_modes{2} = xls_file;
end
%% load data
if strcmp(load_data, 'load')
    [FileName,PathName] = uigetfile({'*.mat'}, 'Select struct file with information on filenames, labels and ID');
    struct_name = [PathName FileName];
    Loadtmp = load(struct_name);
    names = fieldnames(Loadtmp);
    S = Loadtmp.(names{:});
    [FileName,PathName] = uigetfile({'*.mat'}, 'Select file with the VariableNames');
    VariableNames_name = [PathName FileName];
    Loadtmp = load(VariableNames_name);
    names = fieldnames(Loadtmp);
    VariableNames = Loadtmp.(names{:});
end

%% Modes
load_modes{1} = load_data;



