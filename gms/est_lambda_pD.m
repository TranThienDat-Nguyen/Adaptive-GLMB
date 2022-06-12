% % This is the MATLAB code for the lambda-pD-CPHD filter proposed in
% % R. Mahler, B.-T. Vo and B.-N. Vo, "CPHD Filtering with unknown clutter rate and detection profile," IEEE Transactions on Signal Processing, Vol. 59, No. 8, pp. 3497-3513, 2011.
% % http://ba-ngu.vo-au.com/vo/MVV_PHDrobust.pdf
% 
% % ---BibTeX entry
% % @ARTICLE{RobustCPHD,
% % author={R. P. S. Mahler and B.-T. Vo and B.-N. Vo},
% % journal={IEEE Transactions on Signal Processing},
% % title={CPHD Filtering With Unknown Clutter Rate and Detection Profile},
% % year={2011},
% % month={August},
% % volume={59},
% % number={8},
% % pages={3497-3513}} 
% %---

function cphd_out = est_lambda_pD(tt_birth , cphd_in , model , filter ,meas , k)
    %---preprocess the birhts
    m_birth = zeros(model.x_dim , length(tt_birth)) ; 
    P_birth = zeros(model.x_dim , model.x_dim , length(tt_birth)) ; 
    w_birth = zeros(1,length(tt_birth)) ; 
    u_b = zeros(length(tt_birth),1) ; 
    v_b = zeros(length(tt_birth),1) ; 
    for bidx =1 : length(tt_birth)
        m_birth(:,bidx) = tt_birth{bidx}.m ; 
        P_birth(:,:,bidx) = tt_birth{bidx}.P ;
        w_birth(bidx) = tt_birth{bidx}.r_b ;
        u_b(bidx) = tt_birth{bidx}.st(1) ;
        v_b(bidx) = tt_birth{bidx}.st(2) ;
    end
    %---intensity prediction 
    [m_predict,P_predict] = kalman_predict_multiple(model,cphd_in.m_update,cphd_in.P_update);                       %surviving components
    w_predict= model.P_S*cphd_in.w_update;                                                                  %surviving weights

    %--- prediction of beta parameters
    u_predict = zeros(cphd_in.L_update,1); v_predict = zeros(cphd_in.L_update,1);
    for j=1:cphd_in.L_update
        beta_avg_tmp= cphd_in.u_update(j)/(cphd_in.u_update(j)+cphd_in.v_update(j)); 
        beta_var_tmp= (cphd_in.u_update(j)*cphd_in.v_update(j))/((cphd_in.u_update(j)+cphd_in.v_update(j))^2*(cphd_in.u_update(j)+cphd_in.v_update(j)+1)); 
        beta_tht_tmp= (beta_avg_tmp*(1-beta_avg_tmp))/min(filter.beta_factor*beta_var_tmp,beta_avg_tmp*(1-beta_avg_tmp))- 1;
        u_predict(j)= beta_tht_tmp*beta_avg_tmp;  v_predict(j)= beta_tht_tmp*(1-beta_avg_tmp);
    end   
    
    m_predict= cat(2,m_birth,m_predict); P_predict=cat(3,P_birth,P_predict);            %append birth components
    w_predict= cat(1,w_birth',w_predict);                                                      %append birth weights
    L_predict= length(tt_birth)+cphd_in.L_update;                                                              %number of predicted components    
    
    u_predict = cat(1,u_b,u_predict);
    v_predict = cat(1,v_b,v_predict);
    
    %---prediction for clutter intensity
    Lc_predict= model.Lc_birth + cphd_in.Lc_update;
    uc_predict= [model.u_cb; cphd_in.uc_update];
    vc_predict= [model.v_cb; cphd_in.vc_update];
    wc_predict= [model.w_cb; model.clutter_P_S*cphd_in.wc_update];
    
    %---cardinality prediction 
    %surviving cardinality distribution
    survive_cdn_predict = zeros(filter.N_max+1,1);
    survival_factor= (sum(cphd_in.w_update)*model.P_S + sum(cphd_in.wc_update)*model.clutter_P_S)/(sum(cphd_in.w_update)+sum(cphd_in.wc_update));
    if isnan(survival_factor), survival_factor=0; end %catch the degernerate zero case
    for j=0:filter.N_max
        idxj=j+1;
        terms= zeros(filter.N_max+1,1);
        for ell=j:filter.N_max
            idxl= ell+1;
            terms(idxl) = exp(sum(log(1:ell))-sum(log(1:j))-sum(log(1:ell-j))+j*log(survival_factor)+(ell-j)*log(1-survival_factor))*cphd_in.cdn_update(idxl);
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
    if m~=0
        [qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k},model,m_predict,P_predict);
    end
    
    %There is no elementary symmetrci functions in lambda-CPHD filter
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
    w_update= 1/(sum(w_predict)+sum(wc_predict))*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict)*(beta(u_predict,v_predict+1)./beta(u_predict,v_predict)).*w_predict;
    u_update= u_predict; v_update= v_predict+1;
    
    % detection terms
    clut_w_term= zeros(m,1); w_cumsum= sum(w_update);
    start_pt= L_predict+1;
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
        wc_update(start_pt:end_pt)= wc_t(idx);
        
        %update for these terms
        uc_update(start_pt:end_pt,1)= uc_predict(idx)+1;  vc_update(start_pt:end_pt,1)= vc_predict(idx);
        start_pt= end_pt+ 1;
    end
    wc_update= wc_cumsum*wc_update/sum(wc_update);
        
     
    % cardinality update
    for n=0:filter.N_max
        idxn= n+1;
        if n < m
            cdn_update(idxn) = 0;
        else
            cdn_update(idxn) =  cdn_predict(idxn)*Upsilon0(idxn)/(Upsilon0'*cdn_predict);
        end   
    end
    cdn_update = cdn_update/sum(cdn_update);
    
    
    %pruning, merging, capping
    [w_update,m_update,P_update,u_update,v_update,wc_update,uc_update,vc_update]= gaus_prune_3(w_update,m_update,P_update,u_update,v_update,wc_update,uc_update,vc_update,filter.elim_threshold);    

    pD_tmp = (u_update(:)'./(u_update(:)'+v_update(:)'));
    pD = w_update(:)'*pD_tmp'/sum(w_update);
    Nc_update = sum((uc_update(:)'./(uc_update(:)'+vc_update(:)')).*wc_update(:)');
    
    %capping
    [w_update,m_update,P_update,u_update,v_update,wc_update,uc_update,vc_update]= gaus_cap_3(w_update,m_update,P_update,u_update,v_update,wc_update,uc_update,vc_update,filter.L_max);               
    L_cap  = length(w_update); Lc_cap = length(wc_update);
    
    L_update= L_cap; Lc_update = Lc_cap;
    
    % Output assignment
    %---main
    cphd_out.m_update = m_update ;
    cphd_out.w_update = w_update ;
    cphd_out.P_update = P_update ;
    cphd_out.L_update = L_update ;
    cphd_out.cdn_update = cdn_update ;
    %---clutter
    cphd_out.u_update = u_update ; 
    cphd_out.v_update = v_update ; 
    cphd_out.Nc_update = Nc_update  ;
    cphd_out.Lc_update = Lc_update ; 
    cphd_out.uc_update = uc_update ; 
    cphd_out.vc_update = vc_update ; 
    cphd_out.wc_update = wc_update ; 
    %---P_D
    cphd_out.P_D = pD ; 
    
end

            