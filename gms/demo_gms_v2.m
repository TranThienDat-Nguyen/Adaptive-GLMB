% This is a demonstration of the adaptive GLMB filter with linear Gaussian 
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

disp('Running adaptive GLMB filter with linear Gaussian models (v2) ...')
model = gen_model;
truth = gen_truth(model);
meas = gen_meas(model,truth);
est = run_filter_v2(model,meas);
handles = plot_results(model,truth,meas,est);