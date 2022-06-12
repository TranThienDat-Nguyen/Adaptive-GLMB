function est = run_filter_v1(model,meas)
% This is an implementation of the adaptive GLMB filter with non-linear Gaussian 
% models. The filter is proposed in:
% @ARTICLE{adaptive_GLMB,
%   author = {Do, Cong-Thanh and Nguyen, Tran Thien Dat and Moratuwage, Diluka 
%   and Shim, Changbeom and Chung, Yon Dohn},
%   journal = {Signal Processing},
%   title = {Multi-object tracking with an adaptive generalized labeled
%   multi-Bernoulli filter},
%   year = {2022},
%   volume = {196},
%   pages = {108532}}

%=== Setup

%filter parameters
filter.H_upd= 3000;                 %requested number of updated components/hypotheses
filter.H_max= 3000;                 %cap on number of posterior components/hypotheses
filter.hyp_threshold= 1e-10;        %pruning threshold for components/hypotheses

filter.L_max= 100;                  %limit on number of Gaussians in each track - not implemented yet
filter.elim_threshold= 1e-5;        %pruning threshold for Gaussians in each track - not implemented yet
filter.merge_threshold= 4;          %merging threshold for Gaussians in each track - not implemented yet

filter.ukf_alpha= 1;                %scale parameter for UKF - choose alpha=1 ensuring lambda=beta and offset of first cov weight is beta for numerical stability
filter.ukf_beta= 1;                 %scale parameter for UKF
filter.ukf_kappa= 1;                %scale parameter for UKF (alpha=1 preferred for stability, giving lambda=1, offset of beta for first cov weight)

filter.N_max= 100;                  %cap cardinality for cphd

filter.P_G= 0.9999999;                           %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0
filter.beta_factor = 1.1;                        %beta factor

filter.run_flag= 'disp';                   %'disp' or 'silence' for on the fly output
filter.estimator_type = 'partial_smooth' ; %'partial_smooth' or 'standard' estimator
filter.estimator_prune_thres = 5 ;         %prune threshold (trajectory length) for recursive estimator

%initializing structure to store estimates
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);
est.T= {};
est.M= 0;
est.J= []; 
est.H= {};
est.glmb_Nc = zeros(meas.K,1) ; 

%=== Filtering
%initial prior for cphd
cphd.w_update(1)= eps;
cphd.m_update(:,1)= [0.1;0;0.1;0;0];
cphd.P_update(:,:,1)= diag([1 1 1 1 0.1]).^2;
cphd.u_update = 1;
cphd.v_update = 1;
cphd.L_update = 1;

%initialisation of clutter
cphd.PD_init = 0.99;
cphd.Nc_init = round((size(meas.Z{1},2)-cphd.PD_init*0.5*size(meas.Z{1},2))/cphd.PD_init);
cphd.cdn_update = [zeros(cphd.Nc_init-1,1);1;zeros(filter.N_max-cphd.Nc_init+1,1)];

%initialize expected number of clutter generators
cphd.wc_update= cphd.Nc_init;
cphd.uc_update= 1;
cphd.vc_update= 1;
cphd.Lc_update= 1;

%initial prior for glmb
glmb_update.tt= cell(0,1);      %track table for GLMB (cell array of structs for individual tracks)
glmb_update.w= 1;               %vector of GLMB component/hypothesis weights
glmb_update.I= {[]};            %cell of GLMB component/hypothesis labels (labels are indices/entries in track table)
glmb_update.n= 0;               %vector of GLMB component/hypothesis cardinalities
glmb_update.cdn= 1;             %cardinality distribution of GLMB (vector of cardinality distribution probabilities)
glmb_update.clutter = 0 ;
% Initialize birth track for the first frame
tt_birth = cell(size(meas.Z{1},2),1) ; 
for midx = 1 : length(tt_birth)
    tt_birth{midx}.m = meas2state(meas.Z{1}(:,midx)) ; 
    tt_birth{midx}.P = model.P_birth ; 
    tt_birth{midx}.w = 1 ;
    tt_birth{midx}.l = [2;midx] ;
    tt_birth{midx}.r_b = 0.01 ; 
    tt_birth{midx}.ah = [] ; 
    tt_birth{midx}.st = model.init_st ; 
end

est.filter= filter;
%recursive filtering
for k=2:meas.K
    %joint prediction and update
    cphd = est_lambda_pD(tt_birth , cphd, model, filter , meas , k) ; 
    [meas_ru_w , glmb_update] = jointpredictupdate(glmb_update,tt_birth,model,filter,meas,k,cphd.Nc_update,cphd.P_D);      
    H_posterior= length(glmb_update.w);
    tt_birth= gen_meas_driven_birth(model , meas_ru_w, meas.Z{k}, k) ; 
    %pruning and truncation
    glmb_update= prune(glmb_update,filter);
    glmb_update= cap(glmb_update,filter);
    
    %state estimation and display diagnostics
    if strcmp(filter.estimator_type, 'partial_smooth')
        prune_flag = 0 ; 
        prune_thres = filter.estimator_prune_thres ; % threshold (length) to prune short tracjectories
        if k==meas.K % only prune short trajectories at the last time step
            prune_flag = 1 ; 
        end
        [est, est_N_c] = extract_estimates_partial_smooth(glmb_update,model,meas,est,prune_flag,prune_thres,meas.K,filter) ;
        est.glmb_Nc(k) = est_N_c ; 
    elseif strcmp(filter.estimator_type, 'standard')
        [est.X{k},est.N(k),est.L{k},est.glmb_Nc(k)]= extract_estimates(glmb_update,model);
    else
        error('Unknown type of estimator. Check filter.estimator_type')
    end
    est.lambda_c(k)= cphd.Nc_update ; 
    est.num_meas(k) = size(meas.Z{k},2) ;
    display_diaginfo(glmb_update,k,est,filter,H_posterior);
end
% save ukf_est
end

function [meas_ru_w , glmb_nextupdate] = jointpredictupdate(glmb_update,tt_birth,model,filter,meas,k,N_c,P_D)
%---generate next update

%create surviving tracks - via time prediction (single target CK)
tt_survive= cell(length(glmb_update.tt),1);                                                                                 %initialize cell array
for tabsidx=1:length(glmb_update.tt)
    [mtemp_predict,Ptemp_predict]= ukf_predict_multiple(model,glmb_update.tt{tabsidx}.m,glmb_update.tt{tabsidx}.P,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);     %kalman prediction for GM
    tt_survive{tabsidx}.m= mtemp_predict;                                                                                   %means of Gaussians for surviving track
    tt_survive{tabsidx}.P= Ptemp_predict;                                                                                   %covs of Gaussians for surviving track
    tt_survive{tabsidx}.w= glmb_update.tt{tabsidx}.w;                                                                       %weights of Gaussians for surviving track
    tt_survive{tabsidx}.l= glmb_update.tt{tabsidx}.l;                                                                       %track label
    tt_survive{tabsidx}.ah= glmb_update.tt{tabsidx}.ah;                                                                     %track association history (no change at prediction) 
end

%create predicted tracks - concatenation of birth and survival
glmb_predict.tt= cat(1,tt_birth,tt_survive);                                                                                %copy track table back to GLMB struct

%gating by tracks
if filter.gate_flag
    for tabidx=1:length(glmb_predict.tt)
        glmb_predict.tt{tabidx}.gatemeas= gate_meas_ukf_idx(meas.Z{k},filter.gamma,model,glmb_predict.tt{tabidx}.m,glmb_predict.tt{tabidx}.P,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);
    end
else
    for tabidx=1:length(glmb_predict.tt)
        glmb_predict.tt{tabidx}.gatemeas= 1:size(meas.Z{k},2);
    end
end

%precalculation loop for average survival/death probabilities
NB = length(tt_birth) ; 
r_birth = zeros(NB , 1) ; 
for bidx = 1 : NB
    r_birth(bidx) = tt_birth{bidx}.r_b ; 
end
avps= [r_birth; zeros(length(glmb_update.tt),1)];
for tabidx=1:length(glmb_update.tt)
    avps(NB+tabidx)= model.P_S;
end
avqs= 1-avps;

%create updated tracks (single target Bayes update)
m= size(meas.Z{k},2);                                   %number of measurements
avpd= zeros(length(glmb_predict.tt),1);
tt_update= cell((1+m)*length(glmb_predict.tt),1);       %initialize cell array
%missed detection tracks (legacy tracks)
for tabidx= 1:length(glmb_predict.tt)
    tt_update{tabidx}= glmb_predict.tt{tabidx};         %same track table
    tt_update{tabidx}.ah= [tt_update{tabidx}.ah; 0];    %track association history (updated for missed detection)
    avpd(tabidx) = P_D ; 
end
avqd= 1-avpd;
%measurement updated tracks (all pairs)
allcostm= zeros(length(glmb_predict.tt),m);
for tabidx= 1:length(glmb_predict.tt)
    for emm= glmb_predict.tt{tabidx}.gatemeas
        stoidx= length(glmb_predict.tt)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted_tracks*j + i)
        [qz_temp,m_temp,P_temp] = ukf_update_multiple(meas.Z{k}(:,emm),model,glmb_predict.tt{tabidx}.m,glmb_predict.tt{tabidx}.P,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);    %kalman update for this track and this measurement
        w_temp= qz_temp.*glmb_predict.tt{tabidx}.w+eps;                                                                                 %unnormalized updated weights
        tt_update{stoidx}.m= m_temp;                                                                                                    %means of Gaussians for updated track
        tt_update{stoidx}.P= P_temp;                                                                                                    %covs of Gaussians for updated track
        tt_update{stoidx}.w= w_temp/sum(w_temp);                                                                                        %weights of Gaussians for updated track
        tt_update{stoidx}.l = glmb_predict.tt{tabidx}.l;                                                                                %track label
        tt_update{stoidx}.ah= [glmb_predict.tt{tabidx}.ah; emm];                                                                        %track association history (updated with new measurement)
        allcostm(tabidx,emm)= sum(w_temp);                                                                                              %predictive likelihood
    end
end
glmb_nextupdate.tt= tt_update;  



%copy track table back to GLMB struct
%joint cost matrix
jointcostm= [diag(avqs) ...
             diag(avps.*avqd) ...
             repmat(avps.*avpd,[1 m]).*allcostm/(N_c*model.pdf_c)];
%gated measurement index matrix
gatemeasidxs= zeros(length(glmb_predict.tt),m);
for tabidx= 1:length(glmb_predict.tt)
    gatemeasidxs(tabidx,1:length(glmb_predict.tt{tabidx}.gatemeas))= glmb_predict.tt{tabidx}.gatemeas;
end
gatemeasindc= gatemeasidxs>0;
         

%component updates
runidx= 1;
meas_ru = [] ;
for pidx=1:length(glmb_update.w)
    %calculate best updated hypotheses/components
    cpreds= length(glmb_predict.tt);
    nbirths= NB;
    nexists= length(glmb_update.I{pidx});
    ntracks= nbirths + nexists;
    tindices= [1:nbirths nbirths+glmb_update.I{pidx}'];                                                                                 %indices of all births and existing tracks  for current component
    lselmask= false(length(glmb_predict.tt),m); lselmask(tindices,:)= gatemeasindc(tindices,:);                                         %logical selection mask to index gating matrices
    mindices= unique_faster(gatemeasidxs(lselmask));                                                                                    %union indices of gated measurements for corresponding tracks
    costm= jointcostm(tindices,[tindices cpreds+tindices 2*cpreds+mindices]);                                                           %cost matrix - [no_birth/is_death | born/survived+missed | born/survived+detected]
    neglogcostm= -log(costm);                                                                                                           %negative log cost
    [uasses,nlcost]= gibbswrap_jointpredupdt_custom(neglogcostm,round(filter.H_upd*sqrt(glmb_update.w(pidx))/sum(sqrt(glmb_update.w))));%murty's algo/gibbs sampling to calculate m-best assignment hypotheses/components
    uasses(uasses<=ntracks)= -inf;                                                                                                      %set not born/track deaths to -inf assignment
    uasses(uasses>ntracks & uasses<= 2*ntracks)= 0;                                                                                     %set survived+missed to 0 assignment
    uasses(uasses>2*ntracks)= uasses(uasses>2*ntracks)-2*ntracks;                                                                       %set survived+detected to assignment of measurement index from 1:|Z|    
    uasses(uasses>0)= mindices(uasses(uasses>0));                                                                                       %restore original indices of gated measurements
    
    %generate corrresponding jointly predicted/updated hypotheses/components
    for hidx=1:length(nlcost)
        meas_ru_temp = false (size(meas.Z{k},2) , 1) ;
        update_hypcmp_tmp= uasses(hidx,:)'; 
        meas_ru_temp(update_hypcmp_tmp(update_hypcmp_tmp>0)) = true ; 
        meas_ru = [meas_ru , meas_ru_temp] ; 
        update_hypcmp_idx= cpreds.*update_hypcmp_tmp+[(1:nbirths)'; nbirths+glmb_update.I{pidx}];
        glmb_nextupdate.w(runidx)= -N_c+m*log(N_c*model.pdf_c)+log(glmb_update.w(pidx))-nlcost(hidx);                                             %hypothesis/component weight
        glmb_nextupdate.I{runidx}= update_hypcmp_idx(update_hypcmp_idx>0);                                                                                              %hypothesis/component tracks (via indices to track table)
        glmb_nextupdate.n(runidx)= sum(update_hypcmp_idx>0);   
        glmb_nextupdate.clutter(runidx) = m-sum(update_hypcmp_tmp>0) ; %hypothesis/component cardinality
        runidx= runidx+1;
    end
end

glmb_nextupdate.w= exp(glmb_nextupdate.w-logsumexp(glmb_nextupdate.w));                                                                                                                 %normalize weights
meas_ru_w = meas_ru * glmb_nextupdate.w' ; 
%extract cardinality distribution
for card=0:max(glmb_nextupdate.n)
    glmb_nextupdate.cdn(card+1)= sum(glmb_nextupdate.w(glmb_nextupdate.n==card));                                                                                                       %extract probability of n targets
end

%remove duplicate entries and clean track table
glmb_nextupdate= clean_update(clean_predict(glmb_nextupdate));
end

function glmb_out= cap(glmb_in,filter)
%cap total number of components to specified maximum
if length(glmb_in.w) > filter.H_max
    [~,idxsort]= sort(glmb_in.w,'descend');
    idxkeep=idxsort(1:filter.H_max);
    glmb_out.tt= glmb_in.tt;
    glmb_out.w= glmb_in.w(idxkeep);
    glmb_out.I= glmb_in.I(idxkeep);
    glmb_out.n= glmb_in.n(idxkeep);
    glmb_out.clutter= glmb_in.clutter(idxkeep) ; 
    glmb_out.w= glmb_out.w/sum(glmb_out.w);
    for card=0:max(glmb_out.n)
        glmb_out.cdn(card+1)= sum(glmb_out.w(glmb_out.n==card));
    end
else
    glmb_out= glmb_in;
end
end

function glmb_temp= clean_predict(glmb_raw)
%hash label sets, find unique ones, merge all duplicates
for hidx= 1:length(glmb_raw.w)
    glmb_raw.hash{hidx}= sprintf('%i*',sort(glmb_raw.I{hidx}(:)'));
end

[cu,~,ic]= unique(glmb_raw.hash);

glmb_temp.tt= glmb_raw.tt;
glmb_temp.w= zeros(length(cu),1);
glmb_temp.I= cell(length(cu),1);
glmb_temp.n= zeros(length(cu),1);
glmb_temp.clutter = zeros(length(cu),1);
for hidx= 1:length(ic)
        glmb_temp.w(ic(hidx))= glmb_temp.w(ic(hidx))+glmb_raw.w(hidx);
        glmb_temp.I{ic(hidx)}= glmb_raw.I{hidx};
        glmb_temp.n(ic(hidx))= glmb_raw.n(hidx);
        glmb_temp.clutter(ic(hidx))= glmb_raw.clutter(hidx);
end
glmb_temp.cdn= glmb_raw.cdn;
end

function glmb_clean= clean_update(glmb_temp)
%flag used tracks
usedindicator= zeros(length(glmb_temp.tt),1);
for hidx= 1:length(glmb_temp.w)
    usedindicator(glmb_temp.I{hidx})= usedindicator(glmb_temp.I{hidx})+1;
end
trackcount= sum(usedindicator>0);

%remove unused tracks and reindex existing hypotheses/components
newindices= zeros(length(glmb_temp.tt),1); newindices(usedindicator>0)= 1:trackcount;
glmb_clean.tt= glmb_temp.tt(usedindicator>0);
glmb_clean.w= glmb_temp.w;
glmb_clean.clutter= glmb_temp.clutter ; 
for hidx= 1:length(glmb_temp.w)
    glmb_clean.I{hidx}= newindices(glmb_temp.I{hidx});
end
glmb_clean.n= glmb_temp.n;
glmb_clean.cdn= glmb_temp.cdn;
end

function glmb_out= prune(glmb_in,filter)
%prune components with weights lower than specified threshold
idxkeep= find(glmb_in.w > filter.hyp_threshold);
glmb_out.tt= glmb_in.tt;
glmb_out.w= glmb_in.w(idxkeep);
glmb_out.I= glmb_in.I(idxkeep);
glmb_out.n= glmb_in.n(idxkeep);
glmb_out.clutter = glmb_in.clutter(idxkeep) ; 
glmb_out.w= glmb_out.w/sum(glmb_out.w);
for card=0:max(glmb_out.n)
    glmb_out.cdn(card+1)= sum(glmb_out.w(glmb_out.n==card));
end
end

function display_diaginfo(glmb,k,est,filter,H_posterior)
if ~strcmp(filter.run_flag,'silence')
    disp([' time= ',num2str(k),...
        ' #eap cdn=' num2str((0:(length(glmb.cdn)-1))*glmb.cdn(:)),...
        ' #var cdn=' num2str((0:(length(glmb.cdn)-1)).^2*glmb.cdn(:)-((0:(length(glmb.cdn)-1))*glmb.cdn(:))^2,4),...
        ' #est card=' num2str(est.N(k),4),...
        ' #est clutter=' num2str(est.glmb_Nc(k)),...
        ' #comp post=',num2str(H_posterior)]);
end
end
