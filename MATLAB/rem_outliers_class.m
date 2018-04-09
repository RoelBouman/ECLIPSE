function [control_O,S,percentage_retained_idv,percentage_retained_total] =rem_outliers_class(S,Lowerlimit, class)
%% outlier removal for building control model

%Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - data structure (STRUCT) 'S' with dimension 1 x (total_measurements). It must atleast contain:
%   Data input (S.Data)
%   IDs input (S.ID)
%   Labels S.Labels

% - Lowerlimit integer (INT) with dimension 1 x 1
%   Possible values (0 to number_variables)
%   Recommended value 1 to 4
%   Lower values will remove more cells. eg. less conservative

%OUTPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -  data structure (STRUCT) 'S_O' with dimension 1 x (total_measurements)
%   Data (S_O.Data)
%   IDs (S_O.ID)
%   Labels (S_O.Labels)

% - retained percentage of cells per individual
%   'percentage_retained_idv' with dimension (number_individuals x 1)

% - retained percentage average 'percentage_retained_total'
%with dimension (1 x 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate 95% quantiles (5% considered outlier) for each ID,variable.
% If cell of individual has ['Lowerlimit'] or more variables considered outlier, it is
% removed. This is done to correct instrumental errors in FACS data.

Labels=vertcat(S.Labels);
N_Var = length(S(1).Data(1,:));
ID=vertcat(S.ID);

switch class
    case 'control'
        control=vertcat(S(Labels==0));
        ID_contr=vertcat(S(Labels==0).ID);
        N_ID = length(control);
        data = vertcat(control.Data)';
        for i = 1:N_Var
            y(i,:) = quantile(data(i,:),[0.025 0.25 0.50 0.75 0.975]);
        end
        control_temp = control;
        control_temp2 =control;
        control_O = control;
        percentage_retained_idv = zeros(N_ID,1);
        for j = 1:N_ID
            for i = 1:N_Var
                control_temp(j).Data(:,i)=control(j).Data(:,i) <= y(i,1);
                control_temp2(j).Data(:,i)=control(j).Data(:,i) >= y(i,5);
            end
            control_temp(j).Data2(:,1) = ~((sum(control_temp(j).Data,2)+sum(control_temp2(j).Data,2)) > Lowerlimit);
            control_O(j).Data = control(j).Data(control_temp(j).Data2,:);
            percentage_retained_idv(j) = (length(control_O(j).Data) / length(control(j).Data))*100;
        end
        percentage_retained_total = mean(percentage_retained_idv);
        
        S_O=S;
        for l1=1:length(ID_contr)
            S_O(ID==ID_contr(l1)).Data=control_O(ID_contr==ID_contr(l1)).Data;
        end
        
        for i=1:length(S_O)
            S(i).Data=S_O(i).Data;
        end
        
    case 'patient'
        patient=vertcat(S(Labels==1));
        ID_pat=vertcat(S(Labels==1).ID);
        N_ID = length(patient);
        data = vertcat(patient.Data)';
        for i = 1:N_Var
            y(i,:) = quantile(data(i,:),[0.025 0.25 0.50 0.75 0.975]);
        end
        patient_temp = patient;
        patient_temp2 =patient;
        patient_O = patient;
        percentage_retained_idv = zeros(N_ID,1);
        for j = 1:N_ID
            for i = 1:N_Var
                patient_temp(j).Data(:,i)=patient(j).Data(:,i) <= y(i,1);
                patient_temp2(j).Data(:,i)=patient(j).Data(:,i) >= y(i,5);
            end
            patient_temp(j).Data2(:,1) = ~((sum(patient_temp(j).Data,2)+sum(patient_temp2(j).Data,2)) > Lowerlimit);
            patient_O(j).Data = patient(j).Data(patient_temp(j).Data2,:);
            percentage_retained_idv(j) = (length(patient_O(j).Data) / length(patient(j).Data))*100;
        end
        percentage_retained_total = mean(percentage_retained_idv);
        
         S_O=S;
        for l1=length(ID_pat)
            S_O(ID==ID_pat(l1)).Data=patient_O(ID_pat==ID_pat(l1)).Data;
        end
        
    case 'all'
        N_ID = length(S);
        data = vertcat(S.Data)';
        for i = 1:N_Var
            y(i,:) = quantile(data(i,:),[0.025 0.25 0.50 0.75 0.975]);
        end
        S_temp = S;
        S_temp2 =S;
        S_O = S;
        percentage_retained_idv = zeros(N_ID,1);
        for j = 1:N_ID
            for i = 1:N_Var
                S_temp(j).Data(:,i)=S(j).Data(:,i) <= y(i,1);
                S_temp2(j).Data(:,i)=S(j).Data(:,i) >= y(i,5);
            end
            S_temp(j).Data2(:,1) = ~((sum(S_temp(j).Data,2)+sum(S_temp2(j).Data,2)) > Lowerlimit);
            S_O(j).Data = S(j).Data(S_temp(j).Data2,:);
            percentage_retained_idv(j) = (length(S_O(j).Data) / length(S(j).Data))*100;
        end
        percentage_retained_total = mean(percentage_retained_idv);
end