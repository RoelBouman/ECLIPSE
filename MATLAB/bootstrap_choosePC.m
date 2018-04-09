%bootstrap for estabilishing the stability of the PCA space, choose the
%number of the PCs for which the angle between the loadings is the same
%The PCA loadings determined for each resampling are compared for changes.
% Principal components which change significantly from one resampling to the next
% are probably due mostly to noise rather than signal.


% INPUTS
% X:            Data matrix (raw means not centered or scaled)
%
% nb:           number of bootstrap samples (default: no. of rows multiplied by 40)
%               0: no rotation (default)
%               1: rotation based on scores
%               2: rotation based on loadings
%               3: rotation based on a combination of scores and loadings
% PC            number of principal component (choose #variables -1 )
%
% OUTPUTS
% Result_subangle: angle between subspaces (the original one and bootstrapping sample)
% XBtrain:         cells belonging to the train set in the bootstrapping
% XBtest:          cells belonging to the train set in the bootstrapping
% modelPCA:        model, struct containing
% Xp:

% Written by R. Folcarelli

function [Result_subangle, XBtrain, XBtest, modelPCA, Xp]=bootstrap_choosePC(Xp, nb, procrust, PC, pret)

Labels=vertcat(Xp.Labels);
% [Xp, pret] =pret_MFC_Data(X, []);

% switch model
%     case 'control' %  PCA on controls
%         Xp=vertcat(Xp(Labels==0));
%     case 'patient' % PCA on patients;
%         Xp=vertcat(Xp(Labels==1));
%     case 'all'  % PCA on both individuals
%         Xp=Xp;
% end

Xp_control=vertcat(Xp(Labels==0));
% Xp_control_data=vertcat(Xp_control.Data);
Xp_patient=vertcat(Xp(Labels==1));

[Xp_blockscaled]=blockscale(Xp_control);
[T, P, VAR, ~]=PCAecon(Xp_blockscaled);% Normal PCA on tnhe original dataset
Tc=(vertcat(Xp_control.Data))*P;
Tp=(vertcat(Xp_patient.Data))*P;


[I, J]=size(Xp_blockscaled);
TB = nan(I,PC,nb);
PB = nan(J,PC,nb);
Result_subangle=nan(nb,PC);

XBtrain=[];
XBtest=[];
% Q=[];
% T2=[];

modelPCA=struct('Tc',Tc, 'Tp', Tp, 'P',P,'VAR', VAR);

Xbtr=Xp_control;
Xbts=Xp_control;
for i=1:nb
    for k=1:size(Xp_control,1) % bootstrap cells within each individual of the struct X
        ns=size(Xp_control(k).Data,1);
        %         itr=ceil(rand(ns,nb)*ns);
        itr=ceil(rand(ns,nb)*ns);
        its=(1:ns)'; its(itr(:,i))=[];
        Xbtr(k).Data=Xp_control(k).Data(itr(:,i),:);
        Xbtr(k).nb=i;
        Xbts(k).Data=Xp_control(k).Data(its,:);
        Xbts(k).nb=i;
    end
    
    XBtrain=[XBtrain;Xbtr];
    XBtest=[XBtest;Xbts];
    
    [Xb_ptr, ~] =pret_MFC_Data(Xbtr, pret);% for not paired data, log, mean or median center, scaling
    [Xb_pts, ~] =pret_MFC_Data(Xbts, pret);
    Xbtr_data=vertcat(Xb_ptr.Data);
    Xbts_data=vertcat(Xb_pts.Data);
    
    [Xb_ptr_blockscaled]=blockscale(Xb_ptr);
    %     [Xb_pts_blockscaled]=blockscale(Xb_pts);
    [~, Pbtr, ~, ~] = PCAecon(Xb_ptr_blockscaled);%PCA on the current boostrap sample
    Pb=P;
    [~,Qrr] = reorder_reflect(Pb,Pbtr); % reorder and reflect
    
    %     Tbtr = Tbtr*Qrr;
    Pbtr = Pbtr*Qrr;
    Pbtr = Pbtr(:,1:PC);
    Pb= Pb(:,1:PC);
    Tt_proc = T(itr(:,i),1:PC);
    Tbts=Xbts_data*Pbtr;
    
    Tbtr = Xbtr_data*Pbtr;
    Tbtr = Tbtr(:,1:PC);
    
    
    switch procrust
        case 0 %  no rotation
            Qp = eye(PC);
        case 1; % rotation based on scores
            Qp = proc_rot(mnsc(Tt_proc),Tbtr);
        case 2  % rotation based on loadings
            Qp = proc_rot(mnsc(Pb),mnsc(Pbtr));
            %             Qp = proc_rot(mnsc(P),mnsc(Pbtr));
        case 3  % based on a combination of T&P
            x1= P(:,1:PC);
            x2= mnsc(Tt_proc)*sum_squares(Pbtr)/sum_squares(Tbtr);
            y1= Pbtr;
            y2= mnsc(Tbtr)*sum_squares(Pbtr)/sum_squares(Tbtr);
            X = [x1; x2];
            Y = [y1; y2];
            Qp = proc_rot(X,Y);
        case 4 % add by Rita
            P = P(:,1:PC);
            Qp = rotatefactors(Pbtr,'Method','procrustes','Target',P);
    end
    
    Tbtr = Tbtr*Qp;
    Pbtr = Pbtr*Qp;
    Tbts=Xbts_data*Pbtr;    %projection of test on the PCA model
    PB(:,:,i) = Pbtr;
    TB(:,:,i)=Tbtr;
    %     TB_test(:,:,i)=Tbts;
    
    result =angle_subspace(Pb,Pbtr,PC);
    Result_subangle(i,:)=result;
    
end
