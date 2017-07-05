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
  use_radar_ = true;

  // initial state vector - what's in this vector?
  // px [m]
  // py [m]
  // v [m/s]
  // rho [rad]
  // rho_dot [rad/s]
  x_ = VectorXd(5);
  x_.fill(0.0);  // TODO: is it ok to fill with zeros?

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0.0); // TODO: is it ok to fill with zeros?

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30; // TODO: too much??

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30; // TODO: too much??

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
  n_x_ = 5;
  n_aug_ = 7;
  
  // Sigma point weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); i++) {
    weights_(i) = 1 / (2 * (lambda_ + n_aug_));
  }
  
  // Predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_aug_ * 2 + 1);
  
  // Sigma point spreading param
  lambda_ = 0.1;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  
  */
  
  if (!is_initialized_) {
// TODO: init with first measurement?
//    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
//      auto rho = measurement_pack.raw_measurements_[0];
//      auto phi = measurement_pack.raw_measurements_[1];
//      auto ro_dot = measurement_pack.raw_measurements_[2];
//      ekf_.x_ << cos(phi) * rho, sin(phi) * rho, ro_dot * cos(phi), ro_dot * sin(phi);
//    }
//    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
//      ekf_.x_ << measurement_pack.raw_measurements_[0],
//      measurement_pack.raw_measurements_[1],
//      0,
//      0;
//    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  long dt = meas_package.timestamp_ - time_us_;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Generate augmented sigma points
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i + 1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+ 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  // Predict new Xsig_pred
  for (unsigned int i = 0; i < Xsig_aug.cols(); i++) {
    VectorXd col = Xsig_aug.col(i);
    float v = col(2);
    float yaw = col(3);
    float yaw_dot = col(4);
    float nu_a = col(5);
    float nu_yaw_dd = col(6);
    
    VectorXd change = VectorXd(5);
    if (yaw_dot == 0) {
      change << v * cos(yaw) * delta_t,
      v * sin(yaw) * delta_t,
      0,
      yaw_dot * delta_t,
      0;
    } else {
      change << v * (sin(yaw + yaw_dot * delta_t) - sin(yaw)) / yaw_dot,
      v * (-cos(yaw + yaw_dot * delta_t) + cos(yaw)) / yaw_dot,
      0,
      yaw_dot * delta_t,
      0;
    }
    
    VectorXd noise = VectorXd(5);
    noise << delta_t * delta_t * cos(yaw) * nu_a / 2,
    delta_t * delta_t * sin(yaw) * nu_a / 2,
    delta_t * nu_a,
    delta_t * delta_t * nu_yaw_dd / 2,
    delta_t * nu_yaw_dd;
    
    Xsig_pred_.col(i) = col.head(5) + change + noise;
  }
  
  // reset and predict state mean
  x_.fill(0.0);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
  
  // reset and predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    P_ += weights_(i) * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  
  */
  
  //transform sigma points into measurement space
  //calculate mean predicted measurement
  //calculate measurement covariance matrix S
  
  // TODO: should these two belongs to "this" and be saved in object state?
  VectorXd z_pred = VectorXd(3); // mean predicted measurement, works for radar only (3 params)
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1); // create matrix for sigma points in measurement space
  
  z_pred.fill(0.0);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd sp = Xsig_pred_.col(i);
    float px = sp(0);
    float py = sp(1);
    float v = sp(2);
    float yaw = sp(3);
    //    float yawd = sp(4);
    
    VectorXd col = VectorXd(3);
    float sqrt_px_py = sqrt(px * px + py * py);
    col << sqrt_px_py, // rho
    atan2(py, px), // range
    (px * cos(yaw) * v + py * sin(yaw) * v) / sqrt_px_py; // rho dot
    
    // normalize angle to (0, 2 * PI)
    while (col(1)> M_PI) col(1)-=2.*M_PI;
    while (col(1)<-M_PI) col(1)+=2.*M_PI;
    
    Zsig.col(i) = col;
    z_pred += weights_(i) * col;
  }
  
  MatrixXd R = MatrixXd(3,3);
  R << std_radr_ * std_radr_, 0, 0,
  0, std_radphi_ * std_radphi_, 0,
  0, 0, std_radrd_ * std_radrd_;
  
  // initialize predicted measurement covariance matrix S
  MatrixXd S = MatrixXd(3,3);
  
  S.fill(0.0);
  for (int i = 0; i < Zsig.cols(); i++) {
    MatrixXd tmp = Zsig.col(i) - z_pred;
    S += weights_(i) * tmp * tmp.transpose();
  }
  S += R;
  
  MatrixXd Tc = MatrixXd(n_x_, 3); // 3 for radar
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }
  MatrixXd K = Tc * S.inverse();
  VectorXd z = meas_package.raw_measurements_;
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
}
