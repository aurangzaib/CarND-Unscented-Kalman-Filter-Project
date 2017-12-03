# Unscented Kalman Filter Project


| Notes    | |
|:-----------|:-------------|
| **Source Code**  | [https://github.com/aurangzaib/CarND-Unscented-Kalman-Filter-Project](https://github.com/aurangzaib/CarND-Unscented-Kalman-Filter-Project)  |
| **Overview**  | `README.md`  |
| **Setup**  | `setup.md`  |
| **NIS calculations**  | `NIS_Laser.csv`|
|							 |`NIS_Radar.csv`| 
| **How to run**  | `mkdir build && cd build` | 
| |`cmake .. && make`     	|
| |`./UnscentedKF`     		|

## Introduction:

Extended Kalman Filter algorithm is used to predict is the position (px, py) and velocity (vx, vy) based on Laser and Radar data provided by car (simulator).

The accuracy of the prediction is calculated using RMSE method.

The steps of the project are the following:

- Read the sensors data and ground truth values from Laser and Radar.

- Initialize Unscented Kalman Filter.

- Predict State Vector (`x`) and State Covariance (`P`).

- Update `x` and `P` for Laser.

- Update `x` and `P` for Radar.

- Calculate Root Mean Square Error (`RMSE`) to get the accuracy of estimation against the ground truth values.

## Sensor data and Predicted data:

The data provided by Laser sensor is in following form:

| px | py | timestamp | gt_px | gt_px | gt_vx | gt_vy | gt_yho | gt_yho _dot | 
|:-----------|:-------------|:-----------|:-------------|:-----------|:-------------|:-----------|:-------------|:-----------|

The data provided by Radar sensor is in following form:

| tho | psi | rho_dot | timestamp | gt_px | gt_px | gt_vx | gt_vy | gt_yho | gt_yho _dot | 
|:-----------|:-------------|:-----------|:-------------|:-----------|:-------------|:-----------|:-------------|:-----------|:-----------|

The data provided by server using Extended Kalman Filter is in following form:

| est_px | est_py | est_vx | est_vy | mse | 
|:-----------|:-------------|:-----------|:-------------|:-----------|


## Explanation of the code:

The implementation of UKF is divided into 4 files:

`main.cpp`

`ukf.cpp`

`tools.cpp`

Following table summarizes the purpose of each file:

| File | Explanation |
|:-----------|:-------------|
|**main.cpp**| |
|				| Get measurement data of Laser and Radar from simulator |
| 				| Call `ProcessMeasurement` method of class `FusionEKF` | 
|				| Get estimated position (px, py) and velocity (vx, vy) |
|				| Call `CalculateRMSE` method of class `Tools` to get RMSE value |
|**ukf.cpp**||
|`Constructor` | |
| 				| Initialize sensor flags `use_laser_` and `use_radar_` | 
|				| Initialize process noises `std_a_` and `std_yawwd_` | 
|				| Initialize measurement noise std `std_laspx_`, `std_laspy_`, `std_radr_`, `std_radphi_`, `std_radrd_`|
|				| Initialize process noises `std_a_` and `std_yawwd_` | 
|				| Initialize measurement noises `R_laser_` and `R_radar_`|
|`Init`| |
|				|	Initialize `x` and `P` for first measurement |
|`Prediction`| |
|				|	Calculate augmented sigma points `Xsig_aug` and corviarance `P_aug` |
|				|	Calculated predicted sigma points `Xsig_pred_` considering process noises |
|				|	Calculate state mean and covariance `x_` and `P_` |
|`UpdateLidar` & `UpdateRadar`||
|| Transform sigma points in measurement space `Z_sig`|
||Calculate mean predicted measurement `z_pred`|
||Calculate measurement covariance matrix `S`|
||Calculate cross correlation matrix `Tc` and kalman gain `Kc`|
||Calculate update state mean and covariance`x_` and `P_`|
||Caculate NIS for Laser `NIS_laser` and Radar `NIS_radar`|
|**tools.cpp**|  |
| 				| Implement `CalculateRMSE` to find RMSE values | 