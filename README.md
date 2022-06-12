# Multi-Object Tracking with an Adaptive GLMB Filter
This is an implementation of the adaptive GLMB filter proposed in:
```
@ARTICLE{adaptive_GLMB, 
  author={Do, Cong-Thanh and Nguyen, Tran Thien Dat and Moratuwage, Diluka and Shim, Changbeom and Chung, Yon Dohn}, 
  journal={Signal Processing}, 
  title={Multi-object tracking with an adaptive generalized labeled multi-Bernoulli filter}, 
  year={2022},
  volume={196},
  pages={108532}}
```
The paper is available at https://www.sciencedirect.com/science/article/abs/pii/S0165168422000792 \
A pre-print version is available at https://arxiv.org/abs/2008.00413
# Implementation Notes
'v1' is the original filter proposed in the paper. \
'v2' is based on the filtering algorithm proposed in [1] \
Schematic of the 'v2' filter \
<img src="https://github.com/TranThienDat-Nguyen/adaptive-GLMB/blob/main/v2_filter_schematic.png" width="600" height="272"> \
In 'v2' implementation, the detection probability for each track is processed by GLMB filter.
# Run Instructions
Use the file 'demo.m' to run the demonstrations. In that file:
- Line 67 runs demonstration with linear Gaussian models and filter implementation 'v1' 
- Line 68 runs demonstration with linear Gaussian models and filter implementation 'v2' 
- Line 69 runs demonstration with non-linear Gaussian models and filter implementation 'v1' 
- Line 70 runs demonstration with non-linear Gaussian models and filter implementation 'v2' 

You can choose to use either the standard estimator or partial smooth estimator proposed in [2] by setting the parameter 'filter.estimator_type' in files 'gms/run_filter_v1.m', 'gms/run_filter_v2.m', 'ukf/run_filter_v1.m' and 'ukf/run_filter_v2.m'. Partial smooth estimator is used as default.
# Performance Notes
The algorithms occasionally overestimate the cardinality due to the nature of the measurement-driven birth model.
# Acknowledgment
This implementation is based on MATLAB RFS tracking toolbox provided by Prof. Ba-Tuong Vo at http://ba-tuong.vo-au.com/codes.html.
# References
[1] C.-T. Do, T.T.D. Nguyen, and H.V. Nguyen. 2022. "Robust multi-sensor generalized labeled multi-Bernoulli filter" Signal Processing 192, pp. 108368. \
https://www.sciencedirect.com/science/article/pii/S0165168421004059 \
[2] T.T.D. Nguyen, and D.Y. Kim. 2019. "GLMB Tracker with Partial Smoothing" Sensors 19, pp. 4419.\
https://www.mdpi.com/1424-8220/19/20/4419
# Contact
For any queries please contact me at tranthiendat.nguyen@gmail.com.\
Copyright (C) 2022, Tran Thien Dat Nguyen.
