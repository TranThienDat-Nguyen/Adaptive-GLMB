% For more information on the robust CPHD filter, please refer to the
% following paper:
% @ARTICLE{robust_CPHD,
%   author = {R. P. S. Mahler and B.-T. Vo and B.-N. Vo},
%   journal = {IEEE Transactions on Signal Processing},
%   title = {CPHD Filtering With Unknown Clutter Rate and Detection Profile},
%   year = {2011},
%   month = {August},
%   volume = {59},
%   number = {8},
%   pages = {3497-3513}} 
function [cphd_out, qz_temp, m_temp, P_temp] = cphd_update(cphd_in, cdn_in, model, filter, phd, meas, w_in, w_birth, k)
    
    w_predict = phd.w_predict ; 
    m_predict = phd.m_predict ; 
    P_predict = phd.P_predict ; 
    u_predict = phd.u_predict ; 
    v_predict = phd.v_predict ; 
    
    Lc_predict= model.Lc_birth + cphd_in.Lc_update;
    uc_predict= [model.u_cb; cphd_in.uc_update];
    vc_predict= [model.v_cb; cphd_in.vc_update];
    wc_predict= [model.w_cb; model.clutter_P_S*cphd_in.wc_update];
    
    %---cardinality prediction 
    %surviving cardinality distribution
    survive_cdn_predict = zeros(filter.N_max+1,1);
    survival_factor= (sum(w_in)*model.P_S + sum(cphd_in.wc_update)*model.clutter_P_S)/(sum(w_in)+sum(cphd_in.wc_update));
    if isnan(survival_factor), survival_factor=0; end %catch the degernerate zero case
    for j=0:filter.N_max
        idxj=j+1;
        terms= zeros(filter.N_max+1,1);
        for ell=j:filter.N_max
            idxl= ell+1;
            terms(idxl) = exp(sum(log(1:ell))-sum(log(1:j))-sum(log(1:ell-j))+j*log(survival_factor)+(ell-j)*log(1-survival_factor))*cdn_in(idxl);
        end
        survive_cdn_predict(idxj) = sum(terms);
    end

    %predicted cardinality= convolution of birth and surviving cardinality distribution
    cdn_predict = zeros(filter.N_max+1,1);
    for n=0:filter.N_max
        idxn=n+1;
        terms= zeros(filter.N_max+1,1);
        for j=0:n
            idxj= j+1;
            terms(idxj)= exp(-(sum(w_birth)+model.lambda_cb)+(n-j)*log(sum(w_birth)+model.lambda_cb)-sum(log(1:n-j)))*survive_cdn_predict(idxj);
        end
        cdn_predict(idxn) = sum(terms);
    end

    %normalize predicted cardinality distribution
    cdn_predict = cdn_predict/sum(cdn_predict);
    
    %---intensity update
    %number of measurements
    m= size(meas.Z{k},2);
    
    %pre calculation for Kalman update parameters
    L_predict = size(m_predict,2) ;  
    if m~=0
        [qz_temp,m_temp,P_temp] = ukf_update_multiple(meas.Z{k},model,m_predict,P_predict, ...
                filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);
    else
        qz_temp = zeros(L_predict,0) ; 
    end
    
    %There is no elementary symmetric functions in lambda-CPHD filter
    %pre calculation for Upsilon0 and Upsilon1
    
    missed_factor= 1- (sum(w_predict.*(u_predict./(u_predict+v_predict)))+ sum(wc_predict.*(uc_predict./(uc_predict+vc_predict))))/(sum(w_predict)+sum(wc_predict));
    if isnan(missed_factor), missed_factor=0; end %catch the degernerate zero case
    terms0 = zeros(filter.N_max+1,1);
    for n=0:filter.N_max
        idxn= n+1;
        if n < m
            terms0(idxn) = eps(0);
        else
            terms0(idxn) =  exp(sum(log(1:n))-sum(log(1:n-m))+(n-m)*log(missed_factor));
        end   
    end
    Upsilon0 = terms0;
    
    terms1 = zeros(filter.N_max,1);
    for n=0:filter.N_max
        idxn= n+1;
        if n < m+1
            terms1(idxn) = eps(0);
        else
            terms1(idxn) =  exp(sum(log(1:n))-sum(log(1:n-(m+1)))+(n-(m+1))*log(missed_factor));
        end   
    end
    Upsilon1 = terms1;
    
    % misdetection term
    m_update = m_predict;
    P_update = P_predict;   
    w_update= 1/(sum(phd.w_predict)+sum(wc_predict))*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict)*(beta(u_predict,v_predict+1)./beta(u_predict,v_predict)).*w_predict;
    u_update= u_predict; v_update= v_predict+1;
    
    % detection terms
    clut_w_term= zeros(m,1); w_cumsum= sum(w_update);
    start_pt= length(w_predict)+1;
    for ell=1:m
        
       %calculate weights first 
       w_t= beta(u_predict(:)+1,v_predict(:))./beta(u_predict(:),v_predict(:)).*w_predict(:).*qz_temp(:,ell);
       K_t= sum(wc_predict.*(uc_predict./(uc_predict+vc_predict)))*model.pdf_c;
       C= K_t + sum((u_predict(:)'./(u_predict(:)'+v_predict(:)')).*w_predict(:)'.*qz_temp(:,ell)');
       w_t= w_t/C; w_cumsum= w_cumsum + sum(w_t);
       clut_w_term(ell)= K_t/C;
        
       %find terms with weights that won't be truncated at the end
       idx= find( w_t > filter.elim_threshold ); %idx= find( w_t > 0 ); %if you want to do them all
       end_pt= start_pt-1 + length(idx);
       w_update(start_pt:end_pt)= w_t(idx);
        
       %update for these terms
       for j=1:length(idx)
           m_update(:,start_pt-1+j) = m_temp(:,idx(j),ell);
       end
       P_update(:,:,start_pt:end_pt)= P_temp(:,:,idx);
       u_update(start_pt:end_pt,1)= u_predict(idx)+1; v_update(start_pt:end_pt,1)= v_predict(idx);
       start_pt= end_pt+ 1;
    end
    w_update= w_cumsum*w_update/sum(w_update);    
    
    
    %---update for clutter intensity
    wc_update= 1/(sum(w_predict)+sum(wc_predict))*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict)*(beta(uc_predict,vc_predict+1)./beta(uc_predict,vc_predict)).*wc_predict;
    uc_update= uc_predict; vc_update= vc_predict+ 1;
    
    %detection terms (m of them)
    wc_cumsum= sum(wc_update);
    start_pt= Lc_predict+1; 
    for ell=1:m
        
        %calculate weights first
        wc_t= beta(uc_predict(:)+1,vc_predict(:))./beta(uc_predict(:),vc_predict(:)).*wc_predict(:).*model.pdf_c; 
        K_t= sum(wc_predict.*(uc_predict./(uc_predict+vc_predict)))*model.pdf_c;
        C= K_t + sum((u_predict(:)'./(u_predict(:)'+v_predict(:)')).*w_predict(:)'.*qz_temp(:,ell)');
        wc_t= wc_t/C; wc_cumsum= wc_cumsum + sum(wc_t);
        
        %find terms with weights that won't be truncated at the end
        idx= find( wc_t > filter.elim_threshold ); %idx= find( w_t > 0 ); %if you want to do them all
        end_pt= start_pt-1 + length(idx);
        wc_update(start_pt:end_pt,1)= wc_t(idx);
        
        %update for these terms
        uc_update(start_pt:end_pt,1)= uc_predict(idx)+1;  vc_update(start_pt:end_pt,1)= vc_predict(idx);
        start_pt= end_pt+ 1;
    end
    wc_update= wc_cumsum*wc_update/sum(wc_update);
        
     
    % cardinality update

    %pruning, merging, capping
    [~,~,~,~,~,wc_update,uc_update,vc_update]= gaus_prune_3(w_update,m_update,P_update,u_update,v_update,wc_update,uc_update,vc_update,filter.elim_threshold);    

    Nc_update = sum((uc_update(:)'./(uc_update(:)'+vc_update(:)')).*wc_update(:)');
    
    %capping
    [~,~,~,~,~,wc_update,uc_update,vc_update]= gaus_cap_3(w_update,m_update,P_update,u_update,v_update,wc_update,uc_update,vc_update,filter.L_max);               
    Lc_cap = length(wc_update);
    
    Lc_update = Lc_cap;
    
    % Output assignment (only care about clutter)
    cphd_out.Nc_update = Nc_update  ;
    cphd_out.Lc_update = Lc_update ; 
    cphd_out.uc_update = uc_update ; 
    cphd_out.vc_update = vc_update ; 
    cphd_out.wc_update = wc_update ; 