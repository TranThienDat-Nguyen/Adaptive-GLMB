function m_s = urts_smooth_estimate(model,w,m,P,w_s_p,m_s_p,alpha,kappa,beta)
% Only pick the maximum components to smooth (since I only smooth the estimate)
[~,ii] = max(w) ; 
[~,jj] = max(w_s_p) ;
% Prediction
[X_ukf,u]= ut( [m(:,ii); zeros(model.v_dim,1) ], blkdiag(P(:,:,ii),model.Q),alpha,kappa);
X_pred= gen_newstate_fn(model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.v_dim),:) ); 
m_predict = X_pred*u(:);
X_temp_covariance_p= X_pred- repmat(m_predict,[1 length(u)]); % predicted X_temp
u(1)= u(1)+(1-alpha^2+beta);
P_n_p= X_temp_covariance_p*diag(u)*X_temp_covariance_p'; % not smoothed predicted covriance
X_temp_covariance = X_ukf(1:model.x_dim,:) - repmat(m(:,ii),[1 length(u)]) ; % Current time X_temp
C_n_p= X_temp_covariance*diag(u)*X_temp_covariance_p' ; % cross-covariance between 2 time steps
D = C_n_p /P_n_p ; % smooth gain
m_s = m(:,ii) + D*(m_s_p(:,jj) - m_predict) ; % smoothed mean