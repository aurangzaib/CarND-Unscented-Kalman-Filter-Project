#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // should be set to true in UKF initialization
  is_initialized_ = false;

  /**************
   * SENSOR FLAGS
   **************/

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  /**************
   * MATRICES DIMENSIONS
   **************/

  // number of prediction states
  n_x_ = 5;
  // number of prediction augmented states
  n_aug_ = 7;
  // radar measurement noise dimension
  n_z_radar_ = 3;
  // laser measurement noise dimension
  n_z_laser_ = 2;
  // sigma points
  n_sigma_point_ = 2 * n_aug_ + 1;
  // spread parameter, actual value will be assigned when used
  lambda_ = 3 - n_aug_;

  /**************
   * STATE VECTOR AND MATRIX
   **************/

  // initial state vector
  x_ = VectorXd(n_x_);
  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  /**************
   * PROCESS NOISES
   **************/

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2; // was initially 30
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2; // was initially 30

  /**************
   * MEASUREMENT NOISES STD
   **************/

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;
  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;
  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;
  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**************
   * MEASUREMENT NOISES
   **************/

  // measurement noise - laser
  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_.fill(0.0);
  R_laser_(0, 0) = pow(std_laspx_, 2);
  R_laser_(1, 1) = pow(std_laspy_, 2);
  // measurement noise - radar
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_.fill(0.0);
  R_radar_(0, 0) = pow(std_radr_, 2);
  R_radar_(1, 1) = pow(std_radphi_, 2);
  R_radar_(2, 2) = pow(std_radrd_, 2);

  /**************
   * WEIGHT
   **************/

  weights_ = VectorXd(n_sigma_point_);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // reset time
  time_us_ = 0;
  EPS_ = 0.0001;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  // if not initialized
  if (!is_initialized_) {
    Init(meas_pack);
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // compute the time elapsed between the current and previous measurements
  // difference of current and previous timestamps
  // dt - expressed in seconds
  const double dt = (meas_pack.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_pack.timestamp_;
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // Radar updates
  if (meas_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_pack);
  }
    // Laser updates
  else if (meas_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_pack);
  }
}

/**
 * initialize P_ state covariance matrix
 * initialize x_ state prediction vector
 */
void UKF::Init(MeasurementPackage meas_pack) {

  // Laser measurements
  // when sensor type is laser and use_laser is true
  // laser -> px, py
  if (meas_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    const auto px = fabs(meas_pack.raw_measurements_[0]) > EPS_ ? meas_pack.raw_measurements_[0] : EPS_,
        py = fabs(meas_pack.raw_measurements_[1]) > EPS_ ? meas_pack.raw_measurements_[1] : EPS_;

    x_ << px,                   // px
          py,                   // py
          0,                    // v
          0,                    // yaw
          0;                    // yaw_d

    P_.fill(0.0);
    // state covariance matrix
    P_(0, 0) = pow(std_laspx_, 2);
    P_(1, 1) = pow(std_laspy_, 2);
    P_(2, 2) = 1; // std_v
    P_(3, 3) = 1; // std_yho
    P_(4, 4) = 1; // std_yho_dot
  }

    // Radar measurements
    // when sensor type is radar and use_radar is true
    // radar -> rho, psi, rho_dot
  else if (meas_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Convert radar from polar to cartesian coordinates and initialize state.
    const auto rho = meas_pack.raw_measurements_[0],
        phi = meas_pack.raw_measurements_[1],
        rho_dot = meas_pack.raw_measurements_[2];

    const auto px = rho * cos(phi),
        py = rho * sin(phi),
        vx = rho_dot * cos(phi),
        vy = rho_dot * sin(phi);

    x_ << px,                                     // px
        py,                                     // py
        sqrt(pow(vx, 2) + pow(vy, 2)),          // v
        0,                                      // yaw
        0;                                      // yaw_d

    P_.fill(0.0);
    // state covariance matrix
    P_(0, 0) = pow(std_radr_, 2);
    P_(1, 1) = pow(std_radr_, 2);
    P_(2, 2) = 1;                   // std_v
    P_(3, 3) = 1; // std_yho
    P_(4, 4) = 1;  // std_yho_dot
  }

  // Save the initial timestamp for dt calculation
  time_us_ = meas_pack.timestamp_;

  // done initializing, no need to predict or update
  is_initialized_ = true;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
   1- sigma points generation
   2- sigma points augmentation
   3- sigma points prediction
   4- state mean and covariance prediction
  */

  // create vectors and matrices for sigma generation, augmentation and prediction

  // create matrix with predicted sigma points as columns (5x15)
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_point_);

  /***********************************
   *  1- sigma points generation
   ************************************/

  //create augmented mean vector (7x1)
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance (7x15)
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented sigma point matrix (7x15)
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_point_);

  // spread factor
  lambda_ = 3 - n_aug_;

  // initialize X augmented
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  // initialize P augmented
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = pow(std_a_, 2);
  P_aug(6, 6) = pow(std_yawdd_, 2);

  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  double design_factor = sqrt(lambda_ + n_aug_);
  for (int col = 1; col <= n_aug_; col++) {
    VectorXd root_factor = design_factor * A.col(col - 1);
    // n_aug+1 columns
    Xsig_aug.col(col) << x_aug + root_factor;
    // 2n_aug+1 columns
    Xsig_aug.col(col + n_aug_) << x_aug - root_factor;
  }
  // 1st column
  Xsig_aug.col(0) = x_aug;

  /***********************************
   *  2- sigma points prediction
   ************************************/

  // process noise vector Vk
  VectorXd Vk(5);
  // predict sigma points
  // avoid division by zero
  // write predicted sigma points into right column
  for (int col = 0; col < (n_sigma_point_); col++) {

    // common process noise terms
    /*
     * Px
     * Py
     * V
     * yaw
     * yaw_dot
     * Va
     * Vyaw_dot
     */
    double v        = Xsig_aug(2, col),
           yaw      = Xsig_aug(3, col),
           yaw_dot  = Xsig_aug(4, col),
           v_a      = Xsig_aug(5, col),
           v_yaw    = Xsig_aug(6, col);

    double cos_yaw = cos(yaw),
           sin_yaw = sin(yaw);

    // process noise equations:
    Vk <<
          0.5 * pow(dt, 2) * v_a * cos_yaw,
          0.5 * pow(dt, 2) * v_a * sin_yaw,
          dt * v_a,
          (yaw_dot * dt) + (0.5 * pow(dt, 2) * v_yaw),
          (dt * v_yaw);

    // handling division by zero
    if (fabs(yaw_dot) > EPS_) {
      Vk(0) += (v / yaw_dot) * (+sin(yaw + yaw_dot * dt) - sin_yaw);
      Vk(1) += (v / yaw_dot) * (-cos(yaw + yaw_dot * dt) + cos_yaw);
    } else {
      Vk(0) += v * cos_yaw * dt;
      Vk(1) += v * sin_yaw * dt;
    }

    // sigma predicted points
    Xsig_pred_.col(col) << Xsig_aug.col(col).head(5) + Vk;
  }

  /***********************************
   *  3- state mean and covariance prediction
   ************************************/

  // predict state mean
  x_.fill(0.0);
  for (int loop = 0; loop < n_sigma_point_; loop++) {
    x_ += weights_(loop) * Xsig_pred_.col(loop);
  }

  // predict state covariance matrix
  P_.fill(0.0);
  for (int loop = 0; loop < n_sigma_point_; loop++) {
    VectorXd x_diff = Xsig_pred_.col(loop) - x_;
    // angle normalization
    while (x_diff(3) > M_PI)  x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    // state covariance matrix
    P_ += weights_(loop) * x_diff * x_diff.transpose();
  }
}

/**
 * Use lidar data to update the belief about the object's
 * position. Modify the state vector, x_, and covariance, P_.
 * calculate the lidar NIS
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /*
  1- measurement prediction
  2- measurement update
  */

  /***********************************
   *  1- measurement prediction
   ************************************/

  // create matrix for sigma points in measurement space
  MatrixXd Z_sig = MatrixXd(n_z_laser_, n_sigma_point_); // (2x15)
  Z_sig.fill(0.0);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_laser_);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_laser_, n_z_laser_);

  // transform sigma points into measurement space
  z_pred.fill(0.0);
  for (int col = 0; col < n_sigma_point_; col++) {
    float px = Xsig_pred_(0, col);
    float py = Xsig_pred_(1, col);
    Z_sig(0, col) = px;
    Z_sig(1, col) = py;

    // calculate mean predicted measurement
    z_pred += weights_(col) * Z_sig.col(col);
  }

  // calculate measurement covariance matrix S
  S.fill(0.0);
  for (int col = 0; col < n_sigma_point_; col++) {
    VectorXd z_diff = Z_sig.col(col) - z_pred;
    S += weights_(col) * z_diff * z_diff.transpose();
  }
  S += R_laser_;

  /***********************************
   *  2- measurement update
   ************************************/

  // new measurement values
  VectorXd new_z = VectorXd(n_z_laser_);
  new_z << meas_package.raw_measurements_[0], // px
           meas_package.raw_measurements_[1]; // py

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_laser_);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int col = 0; col < n_sigma_point_; col++) {
    // state and measurements difference
    VectorXd x_diff = Xsig_pred_.col(col) - x_;
    VectorXd z_diff = Z_sig.col(col) - z_pred;

    // angle normalization
    while (x_diff(3) > M_PI)  x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    // cross correlation
    Tc += weights_(col) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // measurements difference (residual)
  VectorXd z_diff = new_z - z_pred;

  // update state vector and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // laser normalized innovation squared (NIS) calculation
  double NIS_laser = z_diff.transpose() * S.inverse() * z_diff;

  // save NIS results in a file
  std::ofstream file;
  file.open("../NIS_Laser.csv", ios::app);
  if (file.is_open()) {
    file << NIS_laser << ",\n";
    file.close();
  }
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   1- measurement prediction
   2- measurement update. calculate the radar NIS.
  */

  /***********************************
   *  1- radar measurement prediction
   ************************************/

  // create matrix for sigma points in measurement space
  MatrixXd Z_sig = MatrixXd(n_z_radar_, n_sigma_point_);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);

  // transform sigma points into measurement space
  z_pred.fill(0.0);
  for (int col = 0; col < n_sigma_point_; col++) {

    double px   = Xsig_pred_(0, col),
           py   = Xsig_pred_(1, col),
           v    = Xsig_pred_(2, col),
           yaw  = Xsig_pred_(3, col),
           vx   = v * cos(yaw),
           vy   = v * sin(yaw);

    if (px < EPS_) px = EPS_;
    if (py < EPS_) py = EPS_;

    double rho = sqrt(pow(px, 2) + pow(py, 2));
    // avoid division by zero
    rho = max(fabs(rho), EPS_);
    double phi = atan2(py, px);
    double rho_dot = (px * vx + py * vy) / rho;

    // rho
    Z_sig(0, col) = rho;
    Z_sig(1, col) = phi;
    Z_sig(2, col) = rho_dot;

    // calculate mean predicted measurement
    z_pred += weights_(col) * Z_sig.col(col);
  }

  // normalize rho in z_pred
//  while (z_pred(1) > M_PI)  z_pred(1) -= 2. * M_PI;
//  while (z_pred(1) < -M_PI) z_pred(1) += 2. * M_PI;

  // calculate measurement covariance matrix S
  S.fill(0.0);
  for (int col = 0; col < n_sigma_point_; col++) {
    VectorXd z_diff = Z_sig.col(col) - z_pred;

    //angle normalization
    while (z_diff(1) > M_PI)  z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    // update measurement noise
    S += weights_(col) * z_diff * z_diff.transpose();
  }
  S += R_radar_;

  /***********************************
   *  2- radar measurement update
   ************************************/

  // create matrix for cross correlation Tc - (5x3)
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int col = 0; col < n_sigma_point_; col++) {
    // state and measurements difference
    VectorXd x_diff = Xsig_pred_.col(col) - x_;
    VectorXd z_diff = Z_sig.col(col)      - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI)  z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    while (x_diff(3) > M_PI)  x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    // cross correlation
    Tc += weights_(col) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // new measurement values
  VectorXd z_ = VectorXd(n_z_radar_);
  z_ << meas_package.raw_measurements_[0], // rho
      meas_package.raw_measurements_[1], // phi
      meas_package.raw_measurements_[2]; // rho_dot

  // measurements difference
  VectorXd z_diff = z_ - z_pred;

  // angle normalization
  while (z_diff(1) > M_PI)  z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  // update state vector and covariance
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // radar normalized innovation squared (NIS) calculation
  double NIS_radar = z_diff.transpose() * S.inverse() * z_diff;
  ofstream radar_file;
  radar_file.open("../NIS_Radar.csv", ios::app);
  // save NIS results in a file
  if (radar_file.is_open()) {
    radar_file << NIS_radar << ",\n";
    radar_file.close();
  }
}
