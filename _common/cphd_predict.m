function [phd, w_in, w_birth, cdn_in] = cphd_predict(glmb_in, tt_survive, tt_birth, model, filter)
% note:
% number cphd_in is the number of sensors...
% prediction is common for all cphd
% clutter rate is unique for each cphd
%% common CPHD prediction
% assume only 1 component of gaussian for each tt
% grab the cardinality
cdn_in = zeros(1, filter.N_max+1) ; 
cdn_in(1:length(glmb_in.cdn)) = glmb_in.cdn ; 
cdn_in = cdn_in/sum(cdn_in) ; 

m_predict = zeros(model.x_dim, length(tt_survive)) ; 
P_predict = zeros(model.x_dim, model.x_dim, length(tt_survive)) ; 
w_predict = zeros(1, length(tt_survive)) ; 
u_predict = zeros(1, length(tt_survive)) ;
v_predict = zeros(1, length(tt_survive)) ;
for tidx = 1 : length(tt_survive)
    m_predict(:, tidx) = tt_survive{tidx}.m ; 
    P_predict(:, :, tidx) = tt_survive{tidx}.P ;
    u_predict(:,tidx) = tt_survive{tidx}.st(1) ; 
    v_predict(:,tidx) = tt_survive{tidx}.st(2) ; 
end

for cmp = 1 : length(glmb_in.w)
    w_predict(glmb_in.I{cmp}) = w_predict(glmb_in.I{cmp}) + exp(glmb_in.w(cmp)) ; 
end
% then normalize the mixture weight
% w_predict = w_predict/sum(w_predict) ; 
w_in = w_predict ; 
w_predict = w_predict * model.P_S ; 
%  births stuff (births go first in the array)
m_birth = zeros(model.x_dim, length(tt_birth)) ;
P_birth = zeros(model.x_dim, model.x_dim, length(tt_birth)) ;
u_birth = zeros(1, length(tt_birth)) ;
v_birth = zeros(1, length(tt_birth)) ;
w_birth = zeros(1, length(tt_birth)) ;
for bidx = 1 : length(tt_birth)
    m_birth(:,bidx) = tt_birth{bidx}.m ; 
    P_birth(:,:,bidx) = tt_birth{bidx}.P ;
    u_birth(:,bidx) = tt_birth{bidx}.st(1) ;
    v_birth(:,bidx) = tt_birth{bidx}.st(2) ;
    w_birth(bidx) = tt_birth{bidx}.r_b ;
end
% now convolve births
phd.m_predict= cat(2,m_birth,m_predict); 
phd.P_predict=cat(3,P_birth,P_predict);            %append birth components
phd.w_predict= cat(2,w_birth,w_predict);           %append birth weights                                                          %number of predicted components    
phd.u_predict = cat(2,u_birth,u_predict);
phd.v_predict = cat(2,v_birth,v_predict);


    
