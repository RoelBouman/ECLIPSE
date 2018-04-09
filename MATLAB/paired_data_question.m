function paired_data = paired_data_question
%%
% Function that asks if the data is paired or not.
% Output: 
% paired_data is 0: the data is unpaired
% paired_data is 1: the data is paired
%
% If a wrong input is given, the function gives a warning and sets the
% paired_data to 0.
%
% Written by G.H. Tinnevelt at Radboud University at 19-6-2015
%%
paired_data = input('Is the data paired? (yes/no) ', 's');
if strcmp(paired_data, 'yes')
    paired_data = 1;
elseif strcmp(paired_data, 'no')
    paired_data = 0;
else
    warning('Wrong input for paired data. The data is now considered unpaired.')
    paired_data = 0;
end
end