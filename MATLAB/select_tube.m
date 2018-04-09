function [S_selected, tubenr] = select_tube(S, tubenr)
%%
%
% Function used to select a certain tube
%
% Written by G.H. Tinnevelt on 5-11-2015 at Radboud University Nijmegen
%%
Tubes = vertcat(S.Tube);
if isempty(tubenr)
    s = unique(Tubes);
    s2 = 'Tubes: ';
    for l1 = 1:length(s)
        s2 = [s2 ' \n ' num2str(s(l1))];
    end
    s2 = [s2 '\nWhat Tube should be used: '];
    tubenr = input(s2);
end
tube_selected = ismember(Tubes, tubenr);

S_selected = S(tube_selected);

end