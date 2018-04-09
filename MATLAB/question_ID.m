function IDs=question_ID(S_scaled)

disp('IDs of the individuals')

allPIDs = vertcat(S_scaled.ID);
allLabels=vertcat(S_scaled.Labels);
        
for i = 1:length(allPIDs)
 disp(strcat(int2str(i),': ID: ',int2str(allPIDs(i)),' (Label:  ',int2str(allLabels(i)),')'));
 
end

 IDs=input('Which Individuals do you want to visualize in the Control Model? ');% explampe[1,2] for individuals with ID=1 and 2