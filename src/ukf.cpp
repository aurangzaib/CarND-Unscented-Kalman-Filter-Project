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

  is_initialized_ = false;

  /**************
   * SENSOR FLAGS
   **************/

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  /**************
   * STATE VECTOR AND MATRIX
   **************/

  // initial state vector
  x_ = VectorXd(5);
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0.0);
  
  /**************
   * PROCESS NOISES
   **************/

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2; // was initially 30
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2; // was initially 30

  /**************
   * LASER MEASUREMENT NOISES
   **************/

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;
  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  /**************
   * RADAR MEASUREMENT NOISES
   **************/

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;
  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**************
   * MATRICES DIMENSIONS
   **************/

  // number of prediction states
  n_x_ = 5;
  // number of prediction augmented states
  n_aug_ = 7;
  // spread parameter
  lambda_ = 0;
  // sigma points weights (15x1)
  int n_sigma_ = 2 * n_aug_ + 1;
  // weight
  weights_ = VectorXd(n_sigma_);

  // Initialize NIS for radar and lidar
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;

  // reset time
  time_us_ = 0;
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
  if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_pack);
  }
  // Laser updates
  else {
    UpdateLidar(meas_pack);
  }
}
/**
 * initialize P_ state covariance matrix
 * initialize x_ state prediction vector
 */
void UKF::Init(MeasurementPackage meas_pack) {
  const double EPS = 0.0001; 
 
  // Laser measurements
  // when sensor type is laser and use_laser is true
  // laser -> px, py
  if (meas_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    auto px = fabs(meas_pack.raw_measurements_[0]) > EPS ? meas_pack.raw_measurements_[0] : EPS,
         py = fabs(meas_pack.raw_measurements_[1]) > EPS ? meas_pack.raw_measurements_[1] : EPS;

    x_ << px,                   // px
          py,                   // py
          0,                    // v
          0,                    // yaw
          0;                    // yaw_d

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
    auto rho = meas_pack.raw_measurements_[0],
         phi = meas_pack.raw_measurements_[1],
         rho_dot = meas_pack.raw_measurements_[2];

    auto px = rho * cos (phi),
         py = rho * sin (phi),
         vx = rho_dot * cos(phi),
         vy = rho_dot * sin(phi);

    x_ << px,                                     // px
          py,                                     // py
          sqrt(pow(vx, 2) + pow(vy, 2)),          // v
          0,                                      // yaw
          0;                                      // yaw_d

      // state covariance matrix
      P_(0, 0) = pow(std_radr_, 2);
      P_(1, 1) = pow(std_radr_, 2);
      P_(2, 2) = 1;                   // std_v
      P_(3, 3) = pow(std_radphi_, 2); // std_yho 
      P_(4, 4) = pow(std_radphi_, 2);  // std_yho_dot
  }

  // set weights
  int n_sigma_ = 2 * n_aug_ + 1;
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int loop = 1; loop < n_sigma_; loop++) {
    weights_(loop) = 0.5 / (lambda_ + n_aug_);
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

  const double EPS = 0.0001;
  // create vectors and matrices for sigma generation, augmentation and prediction

  // spread factor
  lambda_ = 3 - n_x_;

  //create sigma point matrix (5x11)
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //create augmented mean vector (7x1)
  VectorXd x_aug = VectorXd(n_aug_);
  
  //create augmented state covariance (7x15)
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented sigma point matrix (7x15)
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create matrix with predicted sigma points as columns (5x15)
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  /***********************************
   *  1- sigma points generation
   ************************************/

  //calculate square root of P
  MatrixXd Proot = P_.llt().matrixL();

  double design_factor = sqrt(lambda_ + n_x_);

  // set sigma points
  for (int col = 1; col <= n_x_; col++) {
    VectorXd root_factor = design_factor * Proot.col(col - 1);
    // nx+1 columns
    Xsig.col(col)        << x_ + root_factor;
    // 2nx+1 columns
    Xsig.col(col + n_x_) << x_ - root_factor;
  }
  // 1st column
  Xsig.col(0) << x_;

  /***********************************
   *  2- sigma points augmentation
   ************************************/

   // spread factor
   lambda_ = 3 - n_aug_;

  // initialize X augmented
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  // initilizae P augmented
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = pow(std_a_, 2);
  P_aug(6, 6) = pow(std_yawdd_, 2);

  // create square root matrix
  Proot = P_aug.llt().matrixL();

  // create augmented sigma points
  design_factor = sqrt(lambda_ + n_aug_);
  for (int col = 1; col <= n_aug_; col++) {
    VectorXd root_factor = design_factor * Proot.col(col - 1);
    // n_aug+1 columns
    Xsig_aug.col(col)          << x_aug + root_factor;
    // 2n_aug+1 columns
    Xsig_aug.col(col + n_aug_) << x_aug - root_factor;
  }
  // 1st column
  Xsig_aug.col(0) = x_aug;

  /***********************************
   *  3- sigma points prediction
   ************************************/

  // process noise vector Vk
  VectorXd Vk(5);
  // predict sigma points
  // avoid division by zero
  // write predicted sigma points into right column
  for (int col = 0; col < (2 * n_aug_ + 1); col++) {

    // common process noise terms
    /* example of a column in Xsig_aug
     * Px
     * Py
     * V
     * yaw
     * yaw_dot
     * Va
     * Vyaw_dot
     */
    double v     = Xsig_aug(2, col),
        yaw      = Xsig_aug(3, col),
        yaw_dot  = Xsig_aug(4, col),
        v_a      = Xsig_aug(5, col),
        v_yaw    = Xsig_aug(6, col),
        cos_yaw  = cos(yaw),
        sin_yaw  = sin(yaw);

    // process noise equations:
    Vk <<
        (0.5) * pow(dt, 2) * v_a * cos_yaw ,
        (0.5) * pow(dt, 2) * v_a * sin_yaw,
        (dt * v_a),
        (yaw_dot * dt) + (0.5 * pow(dt, 2) * v_yaw),
        (dt * v_yaw);

    // handling division by zero
    if (fabs(yaw_dot) > EPS) {
      Vk(0) += (v / yaw_dot) * (sin(yaw + yaw_dot * dt) - sin_yaw);
      Vk(1) += (v / yaw_dot) *-(cos(yaw + yaw_dot * dt) + cos_yaw);
    } else {
      Vk(0) += v * cos_yaw * dt;
      Vk(1) += v * sin_yaw * dt;
    }

    // sigma predicted poitns
    Xsig_pred_.col(col) << Xsig_aug.col(col) + Vk;
  }

  /***********************************
   *  4- state mean and covariance prediction
   ************************************/

  // predict state mean
  for (int loop = 0; loop < 2 * n_aug_ + 1; loop++) {
    x_ += weights_(loop) * Xsig_pred_.col(loop);
  }

  // predict state covariance matrix
  P_.fill(0.0);
  for (int loop = 0; loop < 2 * n_aug_ + 1; loop++) {
    
    VectorXd x_diff = Xsig_pred_.col(loop) - x_;
    
    // angle normalization
    if (x_diff(3) > +M_PI) {
      x_diff(3) -= 2. * M_PI;
    }
    else if (x_diff(3) < -M_PI) {
      x_diff(3) += 2. * M_PI;
    }
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

  // set measurement dimension, lidar can measure px and py
  int n_z = 2;

  // create matrix for sigma points in measurement space
  MatrixXd Z_sig = MatrixXd(n_z, 2 * n_aug_ + 1); // (2x15)

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  z_pred.fill(0.0);

  // transform sigma points into measurement space
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    float px  = Xsig_pred_(0, col);
    float py  = Xsig_pred_(1, col);
    Z_sig(0, col) = px;
    Z_sig(1, col) = py;

    // calculate mean predicted measurement
    z_pred += weights_(col) * Z_sig.col(col);
  }

  // calculate measurement covariance matrix S
  MatrixXd R(n_z, n_z); // (2x2)
  
  R.fill(0.0);
  R(0, 0) = pow(std_laspx_, 2);
  R(1, 1) = pow(std_laspy_, 2);
  
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    VectorXd z_diff = Z_sig.col(col) - z_pred;
    S += weights_(col) * z_diff * z_diff.transpose();
  }
  S += R;

  /***********************************
   *  2- measurement update
   ************************************/

  // new measurement values
  VectorXd new_z = VectorXd(n_z);
  new_z << meas_package.raw_measurements_[0], // px
           meas_package.raw_measurements_[1]; // py

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // calculate cross correlation matrix
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    // state and measurements difference
    VectorXd x_diff = Xsig_pred_.col(col) - x_;
    VectorXd z_diff = Z_sig.col(col) - z_pred;
    // angle normalization
    if (x_diff(3) > +M_PI) x_diff(3) -= 2. * M_PI;
    else if (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    
    // cross correlation
    Tc += weights_(col) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // measurements difference (residual)
  VectorXd z_diff = new_z - z_pred;

  // laser NIS calculation
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  // update state vector and covariance
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
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

  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Z_sig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  z_pred.fill(0.0);

  // transform sigma points into measurement space
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    float px  = Xsig_pred_(0, col),
          py  = Xsig_pred_(1, col),
          v   = Xsig_pred_(2, col),
          yaw = Xsig_pred_(3, col),
          vx  = v * cos(yaw),
          vy  = v * sin(yaw);

    // rho      
    Z_sig(0, col) = sqrt(pow(px, 2) + pow(py, 2));
    // phi
    Z_sig(1, col) = atan2(py, px);
    // rho dot
    Z_sig(2, col) = (px * vx + py * vy) / Z_sig(0, col);

    // calculate mean predicted measurement
    z_pred += weights_(col) * Z_sig.col(col);
  }

  // calculate measurement covariance matrix S
  MatrixXd R(n_z, n_z);
  R.fill(0.0);
  R(0, 0) = pow(std_radr_, 2);
  R(1, 1) = pow(std_radphi_, 2);
  R(2, 2) = pow(std_radrd_, 2);

  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    VectorXd z_diff = Z_sig.col(col) - z_pred;
    
    //angle normalization
    if (z_diff(1) > +M_PI) z_diff(1) -= 2. * M_PI;
    else if (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    
    S += weights_(col) * z_diff * z_diff.transpose();
  }
  S += R;

  /***********************************
   *  2- radar measurement update
   ************************************/

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    // state and measurements difference
    VectorXd x_diff = Xsig_pred_.col(col) - x_;
    VectorXd z_diff = Z_sig.col(col) - z_pred;
    
    // angle normalization
    if (z_diff(1) > +M_PI) z_diff(1) -= 2. * M_PI;
    else if (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    if (x_diff(3) > +M_PI) x_diff(3) -= 2. * M_PI;
    else if (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    
    // cross correlation
    Tc += weights_(col) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // new measurement values
  VectorXd z_ = VectorXd(n_z);
  z_ << meas_package.raw_measurements_[0], // rho
        meas_package.raw_measurements_[1], // phi
        meas_package.raw_measurements_[2]; // rho_dot

  // measurements difference
  VectorXd z_diff = z_ - z_pred;

  // angle normalization
  if (z_diff(1) > +M_PI) z_diff(1) -= 2. * M_PI;
  else if (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  // radar NIS calculation
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  // update state vector and covariance
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}
