%% Info:
% This is a demonstration of the adaptive GLMB filter proposed in:
% @ARTICLE{adaptive_GLMB,
%   author = {Do, Cong-Thanh and Nguyen, Tran Thien Dat and Moratuwage, Diluka 
%   and Shim, Changbeom and Chung, Yon Dohn},
%   journal = {Signal Processing},
%   title = {Multi-object tracking with an adaptive generalized labeled
%   multi-Bernoulli filter},
%   year = {2022},
%   volume = {196},
%   pages = {108532}}
%
% The 'v2' implementation shares similar structure to the robust multi-sensor
% GLMB filter proposed in:
% @ARTICLE{robust_MSGLMB,
%   author = {Do, Cong-Thanh and Nguyen, Tran Thien Dat and Nguyen, Hoa Van},
%   journal = {Signal Processing},
%   title = {Robust multi-sensor generalized labeled multi-Bernoulli filter},
%   year = {2022},
%   volume = {192},
%   pages = {108368}}
%
% The original joint predict-update GLMB filter is proposed in:
% @ARTICLE{GLMB,
%   author = {B.-N. Vo and B.-T. Vo and H. Hung},
%   journal ={IEEE Transactions on Signal Processing},
%   title = {An Efficient Implementation of the Generalized Labeled Multi-Bernoulli Filter},
%   year = {2017},
%   month = {Apr}
%   volume = {65},
%   number = {8},
%   pages = {1975-1987}}
%
% The original robust CPHD filter is proposed in:
% @ARTICLE{robust_CPHD,
%   author = {R. P. S. Mahler and B.-T. Vo and B.-N. Vo},
%   journal = {IEEE Transactions on Signal Processing},
%   title = {CPHD Filtering With Unknown Clutter Rate and Detection Profile},
%   year = {2011},
%   month = {August},
%   volume = {59},
%   number = {8},
%   pages = {3497-3513}} 
%
% The partial smoothing estimator is proposed in:
% @ARTICLE{GLMB_PartialSmooth,
%   author = {Nguyen, Tran Thien Dat and Kim, Du Yong},
%   title = {GLMB Tracker with Partial Smoothing},
%   journal = {Sensors},
%   year = {2019},
%   volume = {19},
%   number = {20},
%   pages = {4419},
%   publisher = {Multidisciplinary Digital Publishing Institute}}
%
% Version date: 12/06/2022.
% Contact: tranthiendat.nguyen@gmail.com
%% Prepare workspace
clc
close all
addpath(genpath(pwd))
disp('Preparing workspace ...')
clear all ;
dbstop if error ;
warning('off','all') ;
%% Run demonstrations
demo_gms_v1 % linear dynamic scenario version 1
demo_gms_v2 % linear dynamic scenario version 2
demo_ukf_v1 % non-linear dynamic scenario version 1
demo_ukf_v2 % non-linear dynamic scenario version 2