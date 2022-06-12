function model= gen_model

% basic parameters
model.x_dim= 4;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector

% dynamical model parameters (CV model)
model.T= 1;                                           %sampling period
model.A0= [ 1 model.T; 0 1 ];                         %transition matrix                     
model.F= [ model.A0 zeros(2,2); zeros(2,2) model.A0 ];
model.B0= [ (model.T^2)/2; model.T ];
model.B= [ model.B0 zeros(2,1); zeros(2,1) model.B0 ];
model.sigma_v = 2;
model.Q= (model.sigma_v)^2* model.B*model.B';         %process noise covariance

% survival/death parameters
model.P_S= .99;
model.Q_S= 1-model.P_S;

% birth parameters
B_birth= diag([ 10; 10; 10; 10 ]);     % std of Gaussians
model.P_birth = B_birth*B_birth';      % cov of Gaussians


% observation model parameters (noisy x/y only)
model.H= [ 1 0 0 0 ; 0 0 1 0 ];    %observation matrix
model.D= diag([ 10; 10 ]); 
model.R= model.D*model.D';         %observation noise covariance

% detection parameters
model.P_D= .95;             %to generate scenario, not known to the filter
model.Q_D= 1-model.P_D;     %probability of missed detection in measurements
model.init_st = [90 , 10] ; %beta parameters of birth tracks
% clutter parameters
model.range_c = [ -1000 1000; -1000 1000 ];                  %uniform clutter region
model.pdf_c = 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
model.lambda_cb = 1;                                         %birth rate for clutter targets
model.w_cb = 1;
model.clutter_P_S = 0.9;                                     %survival probability for clutter targets
model.clutter_P_D = 0.95;                                    %detection probability for clutter targets
model.Lc_birth = 1; model.u_cb = 1; model.v_cb = 1;          %parameters of birth clutter
model.clutter_Nt = 300;                                      %number of clutter generators
model.lambda_c = 30 ;                                        %clutter rate to generate scenario, not known to the filter
