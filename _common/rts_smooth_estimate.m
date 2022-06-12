function m_s = rts_smooth_estimate(model,w,m,P,w_s_p,m_s_p)
% Only pick the maximum components to smooth (since I only smooth the estimate)
[~,ii] = max(w) ; 
[~,jj] = max(w_s_p) ;
% Prediction
m_predict = model.F * m(:,ii) ; 
P_predict = model.F * P(:,:,ii) * model.F' + model.Q ;
% Linear smoothing
D = P(:,:,ii) * model.F' / P_predict ;  
m_s = m(:,ii) + D*(m_s_p(:,jj) - m_predict) ; 
