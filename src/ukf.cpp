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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2; // was initially 30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2; // was initially 30

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  // reset time
  time_us_ = 0;
  // number of prediction states
  n_x_ = 5;
  // number of prediction augmented states
  n_aug_ = 7;

  int n_sigma_ = 2 * n_aug_ + 1;
  // spread parameter
  lambda_ = 3 - n_aug_;

  // sigma points weights (15x1)
  weights_ = VectorXd(n_sigma_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  // if not initialized
  // initialize P_ state covariance matrix
  // initialize x_ state prediction vector
  // based on radar or laser
  if (!is_initialized_) {

    P_(0, 0) = 0.85; // as std_laspx_ is 0.15
    P_(1, 1) = 0.85; // as std_laspy_ is 0.15
    P_(2, 2) = 1;
    P_(3, 3) = 1;
    P_(4, 4) = 1;

    // Radar measurements
    if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      auto  rho = meas_pack.raw_measurements_[0],
            phi = meas_pack.raw_measurements_[1],
            rho_dot = meas_pack.raw_measurements_[2];
      auto  vx = rho_dot*cos(phi),
            vy = rho_dot*sin(phi);
      x_ << rho * cos(phi),                                     // px
            rho * sin(phi),                                     // py
            sqrt(pow(vx, 2) + pow(vy, 2)),                      // v
            0,                                                  // yaw
            0;                                                  // yaw_d
    }

    // LASER measurements
    else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
      auto  px = fabs(meas_pack.raw_measurements_[0]) > 0.001 ? meas_pack.raw_measurements_[0]: 0.001,
            py = fabs(meas_pack.raw_measurements_[1]) > 0.001 ? meas_pack.raw_measurements_[1]: 0.001;
      x_ << px,                   // px
            py,                   // py
            0,                    // v
            0,                    // yaw
            0;                    // yaw_d
    }

    //set weights
    int n_sigma_ = 2 * n_aug_ + 1;
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int loop = 1; loop < n_sigma_; loop++) {
      weights_(loop) = 0.5 / (lambda_ + n_aug_);
    }

    // Save the initial timestamp for dt calculation
    time_us_ = meas_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // compute the time elapsed between the current and previous measurements
  // difference of current and previous timestamps
  const double delta_t = (meas_pack.timestamp_ - time_us_) / 1000000.0;    //dt - expressed in seconds
  time_us_ = meas_pack.timestamp_;

  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  cout << "sensor type: " << meas_pack.sensor_type_ << endl;
  if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_pack);
  } else {
    // Laser updates
    UpdateLidar(meas_pack);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Estimate the object's location. Modify the state vector x_.
  Predict sigma points, the state, and the state covariance matrix.

  Steps:
   1- sigma points generation
   2- sigma points augmentation
   3- sigma points prediction
   4- state mean and covariance prediction
  */

  cout << "inside prediction" << endl;

  // 1- sigma points generation

  // number of sigma points
  int n_sigma_ = 2 * n_aug_ + 1;
  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, n_sigma_);
  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  double design_factor = sqrt(lambda_ + n_x_);
  for (int col = 1; col <= n_x_; col++) {
    VectorXd root_factor = design_factor * A.col(col - 1);
    Xsig.col(col) << x_+ root_factor;             // nx+1 columns
    Xsig.col(col + n_x_) << x_ - root_factor;     // 2nx+1 columns
  }
  Xsig.col(0) << x_;                              // 1st column

  cout << "Xsigma: \n" << Xsig << endl;

  // 2- sigma points augmentation

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);
  //create augmented mean state
  //x_aug is vector of 7 elements
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  //create augmented covariance matrix
  //P_aug is matrix of 7x15 elements
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = pow(std_a_, 2);
  P_aug(6, 6) = pow(std_yawdd_, 2);
  //create square root matrix
  A = P_aug.llt().matrixL();
  //create augmented sigma points
  design_factor = sqrt(lambda_ + n_aug_);
  for (int col = 1; col <= n_aug_; col++) {
    VectorXd root_factor = design_factor * A.col(col - 1);
    Xsig_aug.col(col) << x_aug + root_factor;              // n_aug+1 columns
    Xsig_aug.col(col + n_aug_) << x_aug - root_factor;     // 2n_aug+1 columns
  }
  Xsig_aug.col(0) = x_aug;                                 // 1st column

  cout << "X sigma aug: \n" << Xsig_aug << endl;

  // 3- sigma points prediction

  // create matrix with predicted sigma points as columns (5x15)
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  // process noise vector
  VectorXd processNoise(5);
  // predict sigma points
  // avoid division by zero
  // write predicted sigma points into right column
  for (int matrix_col = 0; matrix_col < (2 * n_aug_ + 1); matrix_col++) {
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
    double v     = Xsig_aug(2, matrix_col),
        yaw      = Xsig_aug(3, matrix_col),
        yaw_dot  = Xsig_aug(4, matrix_col),
        v_a      = Xsig_aug(5, matrix_col),
        v_yaw    = Xsig_aug(6, matrix_col),
        cos_yaw  = cos(yaw),
        sin_yaw  = sin(yaw);
    // process noise equations:
    processNoise <<
        (0.5) * pow(delta_t, 2) * cos_yaw * v_a,
        (0.5) * pow(delta_t, 2) * sin_yaw * v_a,
        (delta_t * v_a),
        (yaw_dot * delta_t) + (0.5 * pow(delta_t, 2) * v_yaw),
        (delta_t * v_yaw);
    // handling zero cases
    if (fabs(yaw_dot) > 0.001) {
      processNoise(0) += (v / yaw_dot) * (sin(yaw + yaw_dot * delta_t) - sin_yaw);
      processNoise(1) += (v / yaw_dot) * (-cos(yaw + yaw_dot * delta_t) + cos_yaw);
    } else {
      processNoise(0) += v * cos_yaw * delta_t;
      processNoise(1) += v * sin_yaw * delta_t;
    }
    // loop through rows. assign state and process values
    for (int matrix_row = 0; matrix_row < n_x_; matrix_row++) {
      Xsig_pred_(matrix_row, matrix_col) = Xsig_aug(matrix_row, matrix_col);
      Xsig_pred_(matrix_row, matrix_col) += processNoise(matrix_row);
    }
  }
  cout << "X sigma predicted: \n" << Xsig_pred_ << endl;

  // 4- state mean and covariance prediction

  // predict state mean
  for (int loop = 0; loop < n_sigma_; loop++) {
    x_ += weights_(loop) * Xsig_pred_.col(loop);
  }

  cout << "x: \n" << x_ << endl;

  //predict state covariance matrix
  P_.fill(0.0);
  for (int loop = 0; loop < n_sigma_; loop++) {
    VectorXd x_diff = Xsig_pred_.col(loop) - x_;
    cout << "x diff: " << x_diff << endl;
    //angle normalization
    while (x_diff(3) > +M_PI) {
      x_diff(3) -= 2. * M_PI;
    }
    while (x_diff(3) < -M_PI) {
      x_diff(3) += 2. * M_PI;
    }
    // state covariance matrix
    P_ += weights_(loop) * x_diff * x_diff.transpose();
  }

  cout << "P: \n" << P_ << endl;

  cout << "end prediction" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  /*
  Steps:
  1- measurement prediction
  2- measurement update
  */

  // 1- measurement prediction

  //set measurement dimension, lidar can measure px and py
  int n_z = 2;
  //create matrix for sigma points in measurement space
  MatrixXd Z_sig = MatrixXd(n_z, 2 * n_aug_ + 1); // (2x15)
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  //transform sigma points into measurement space
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    float px  = Xsig_pred_(0, col),
          py  = Xsig_pred_(1, col);
    Z_sig(0, col) = px;
    Z_sig(1, col) = py;
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    z_pred += weights_(col) * Z_sig.col(col);
  }
  //calculate measurement covariance matrix S
  MatrixXd R(n_z, n_z);
  R.fill(0.0);
  R(0, 0) = pow(std_laspx_, 2);
  R(1, 1) = pow(std_laspy_, 2);
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    VectorXd z_diff = Z_sig.col(col) - z_pred;
    S += weights_(col) * z_diff * z_diff.transpose();
  }
  S += R;

  // 2- measurement update

  // new measurement values
  VectorXd z_ = VectorXd(n_z);
  z_ << meas_package.raw_measurements_[0], // px
        meas_package.raw_measurements_[1]; // py
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    // state and measurements difference
    VectorXd x_diff = Xsig_pred_.col(col) - x_;
    VectorXd z_diff = Z_sig.col(col) - z_pred;
    // angle normalization
    while (x_diff(3) > +M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    // cross correlation
    Tc += weights_(col) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // measurements difference
  VectorXd z_diff = z_ - z_pred;
  // update state vector and covariance
  x_ += K * z_diff;
  P_ += -K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.

  Steps:
   1- measurement prediction
   2- measurement update
  */

  // 1- measurement prediction

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Z_sig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  //transform sigma points into measurement space
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    float px  = Xsig_pred_(0, col),
          py  = Xsig_pred_(1, col),
          v   = Xsig_pred_(2, col),
          yaw = Xsig_pred_(3, col);
    Z_sig(0, col) = sqrt(pow(px, 2) + pow(py, 2));
    Z_sig(1, col) = atan(py / px);
    Z_sig(2, col) = (px * cos(yaw) * v + py * sin(yaw) * v) / Z_sig(0, col);
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    z_pred += weights_(col) * Z_sig.col(col);
  }
  //calculate measurement covariance matrix S
  MatrixXd R(n_z, n_z);
  R.fill(0.0);
  R(0, 0) = pow(std_radr_, 2);
  R(1, 1) = pow(std_radphi_, 2);
  R(2, 2) = pow(std_radrd_, 2);
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    VectorXd z_diff = Z_sig.col(col) - z_pred;
    //angle normalization
    while (z_diff(1) > +M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    S += weights_(col) * z_diff * z_diff.transpose();
  }
  S += R;

  // 2- measurement update

  // new measurement values
  VectorXd z_ = VectorXd(n_z);
  z_ << meas_package.raw_measurements_[0], // rho
       meas_package.raw_measurements_[1], // phi
       meas_package.raw_measurements_[2]; // rho_dot
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
    // state and measurements difference
    VectorXd x_diff = Xsig_pred_.col(col) - x_;
    VectorXd z_diff = Z_sig.col(col) - z_pred;
    // angle normalization
    while (z_diff(1) > +M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    while (x_diff(3) > +M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    // cross correlation
    Tc += weights_(col) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // measurements difference
  VectorXd z_diff = z_ - z_pred;
  // angle normalization
  while (z_diff(1) > +M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
  // update state vector and covariance
  x_ += K * z_diff;
  P_ += -K * S * K.transpose();

  cout << "updated X radar\n: " << x_ << endl;
  cout << "updated P radar\n: " << P_ << endl;
}
