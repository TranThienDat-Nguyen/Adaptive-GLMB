function [est, est_N_c]=extract_estimates_partial_smooth(glmb,model,meas,est,prune_flag,prune_thres,K,filter)
% This is an implementation of the partial smoothing technique proposed in:
% @ARTICLE{GLMB_PartialSmooth,
%   author    = {Nguyen, Tran Thien Dat and Kim, Du Yong},
%   title     = {GLMB Tracker with Partial Smoothing},
%   journal   = {Sensors},
%   year      = {2019},
%   volume    = {19},
%   number    = {20},
%   pages     = {4419},
%   publisher = {Multidisciplinary Digital Publishing Institute}}

%extract MAP cardinality and corresponding highest weighted component
[~,mode] = max(glmb.cdn); 
M = mode-1;
T= cell(M,1);
J= zeros(2,M);

[~,idxcmp]= max(glmb.w.*(glmb.n==M));
for m=1:M
    idxptr= glmb.I{idxcmp}(m);
    T{m,1}= glmb.tt{idxptr}.ah;
    J(:,m)= glmb.tt{idxptr}.l;
end
est_N_c = glmb.clutter(idxcmp) ; 

H= cell(M,1);
for m=1:M
   H{m}= [num2str(J(1,m)),'.',num2str(J(2,m))]; 
end

%compute dead & updated & new tracks
[~,~,is]= intersect(est.H,H);
[~,id,in]= setxor(est.H,H);

est.M= M;
est.T= cat(1,est.T(id),T(is),T(in));
est.J= cat(2,est.J(:,id),J(:,is),J(:,in));
est.H= cat(1,est.H(id),H(is),H(in));

%write out estimates in standard format
est.N= zeros(meas.K,1);
est.X= cell(meas.K,1);
est.L= cell(meas.K,1);
for t=1:length(est.T)
    tah= est.T{t};
    traj_length = length(tah) ; 
    ks= est.J(1,t);
    k_end = ks+traj_length-1 ;
    if prune_flag
        if k_end<K && traj_length<=prune_thres
            continue
        end
    end
    kidx = est.J(1,t);
    midx = est.J(2,t); %index of measurement used to initialise this trajectory
    
    m_cell = cell(length(tah),1) ; 
    P_cell = cell(length(tah),1) ; 
    w_cell = cell(length(tah),1) ; 

    % assign birth 
    % (note: birth index in the label is the same as index of measurement 
    % used to initialize the trajectory)
    z_birth = meas.Z{kidx-1}(:, midx) ; %-1 due to birth delay
    m_cell{1} = meas2state(z_birth);
    P_cell{1} = model.P_birth;
    w_cell{1} = 1;
    % forward filtering of trajectory
    for u=2:length(tah)
        [m_cell{u},P_cell{u}] = ukf_predict_multiple(model,m_cell{u-1},P_cell{u-1},filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);
        w_cell{u} = w_cell{u-1} ; 
        k= ks+u-1;
        emm= tah(u);
        if emm > 0
            [qz,m_cell{u},P_cell{u}] = ukf_update_multiple(meas.Z{k}(:,emm),model,m_cell{u},P_cell{u}, ...
                filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);
            w_cell{u}= qz.*w_cell{u}+eps;
            w_cell{u}= w_cell{u}/sum(w_cell{u});
        end

    end
    % backward smoothing of trajectory 
    [~,idxtrk] = max(w_cell{traj_length}) ; 
    est.N(k_end)= est.N(k_end)+1;
    est.X{k_end}= cat(2,est.X{k_end},m_cell{traj_length}(:,idxtrk));
    est.L{k_end}= cat(2,est.L{k_end},est.J(:,t));
    for u= length(tah)-1 : -1 : 1
        m = urts_smooth_estimate(model,w_cell{u},m_cell{u},P_cell{u}, ...
            w_cell{u+1},m_cell{u+1},filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta) ;
        k= ks+u-1;
        est.N(k)= est.N(k)+1;
        est.X{k}= cat(2,est.X{k},m);
        est.L{k}= cat(2,est.L{k},est.J(:,t));
    end
end
end