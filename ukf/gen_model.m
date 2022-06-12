function model= gen_model
% basic parameters
model.x_dim= 5;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector
model.v_dim= 3;   %dimension of process noise
model.w_dim= 2;   %dimension of observation noise

% dynamical model parameters (CT model)
% state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
model.T= 1;                         %sampling period
model.sigma_vel= 2;
model.sigma_turn= (pi/180);   %std. of turn rate variation (rad/s)
model.bt= model.sigma_vel*[ (model.T^2)/2; model.T ];
model.B2= [ model.bt zeros(2,2); zeros(2,1) model.bt zeros(2,1); zeros(1,2) model.T*model.sigma_turn ];
model.B= eye(model.v_dim);
model.Q= model.B*model.B';

% survival/death parameters
model.P_S= .99;
model.Q_S= 1-model.P_S;

% birth parameters (LMB birth model, single component only)
B_birth = diag([ 10; 10; 10; 10; 3*(pi/180) ]); 
model.P_birth = B_birth * B_birth' ; 

% observation model parameters (noisy r/theta only)
% measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
model.D= diag([ 1*(pi/180); 5 ]);       %std for angle and range noise
model.R= model.D*model.D';              %covariance for observation noise

% detection parameters
model.P_D= .95;             %probability of detection in measurements
model.Q_D= 1-model.P_D;     %probability of missed detection in measurements
model.init_st = [90 , 10] ; %beta parameters of birth tracks

% clutter parameters
% model.lambda_c= 15;                                       %poisson average rate of uniform clutter (per scan)
model.range_c= [ -pi/2 pi/2; 0 2000 ];                      %uniform clutter on r/theta
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density

model.lambda_cb = 1;                                %birth rate for clutter targets
model.w_cb = 1;
model.clutter_P_S = 0.9;                            %survival probability for clutter targets
model.clutter_P_D = 0.95;                           %detection probability for clutter targets
model.Lc_birth = 1; model.u_cb = 1; model.v_cb = 1; %parameters of birth clutter
model.clutter_Nt = 300;                             %number of clutter generators
model.lambda_c= 30 ;                                %clutter rate to generate scenario, not known to the filter




