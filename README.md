# Multi-Object Tracking with an Adaptive GLMB Filter
This is an implementation of the adaptive GLMB filter proposed in:\
@ARTICLE{adaptive_GLMB, \
  &nbsp;&nbsp;&nbsp;&nbsp;author={Do, Cong-Thanh and Nguyen, Tran Thien Dat and Moratuwage, Diluka and Shim, Changbeom and Chung, Yon Dohn}, \
  &nbsp;&nbsp;&nbsp;&nbsp;journal={Signal Processing},  \
  &nbsp;&nbsp;&nbsp;&nbsp;title={Multi-object tracking with an adaptive generalized labeled multi-Bernoulli filter}, \
  &nbsp;&nbsp;&nbsp;&nbsp;year={2022},\
  &nbsp;&nbsp;&nbsp;&nbsp;volume={196},\
  &nbsp;&nbsp;&nbsp;&nbsp;pages={108532}}\
The paper is available at https://www.sciencedirect.com/science/article/abs/pii/S0165168422000792 \
A pre-print version is available at https://arxiv.org/abs/2008.00413
# Implementation Notes
'v1' is the original filter proposed in the paper. \
'v2' is based on the filtering algorithm proposed in the following paper: \
@ARTICLE{robust_MSGLMB, \
  &nbsp;&nbsp;&nbsp;&nbsp;author={Do, Cong-Thanh and Nguyen, Tran Thien Dat and Nguyen, Hoa Van}, \
  &nbsp;&nbsp;&nbsp;&nbsp;journal={Signal Processing},  \
  &nbsp;&nbsp;&nbsp;&nbsp;title={Robust multi-sensor generalized labeled multi-Bernoulli filter}, \
  &nbsp;&nbsp;&nbsp;&nbsp;year={2022},\
  &nbsp;&nbsp;&nbsp;&nbsp;volume={192},\
  &nbsp;&nbsp;&nbsp;&nbsp;pages={108368}}
# Run Instructions
Use the file 'demo.m' to run the simulation. In that file:
- Line 67 runs demonstration with linear Gaussian models and filter implementation 'v1' 
- Line 68 runs demonstration with linear Gaussian models and filter implementation 'v2' 
- Line 69 runs demonstration with non-linear Gaussian models and filter implementation 'v1' 
- Line 70 runs demonstration with non-linear Gaussian models and filter implementation 'v2' 
# Performance Notes
The algorithms occasionally overestimate the cardinality due to the nature of the measurement-driven birth model.
# Acknowledgments
This implementation is based on MATLAB RFS tracking toolbox provided by Prof. Ba-Tuong Vo at http://ba-tuong.vo-au.com/codes.html.
# Contact
For any queries please contact me at tranthiendat.nguyen@gmail.com.\
Copyright (C) 2022, Tran Thien Dat Nguyen.
