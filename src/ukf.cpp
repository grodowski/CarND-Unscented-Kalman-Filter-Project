#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

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
  x_ = VectorXd::Zero(5);
  
  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 10;

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

  lambda_ = 3 - n_aug_;
  
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); i++) {
    weights_(i) = 1 / (2 * (lambda_ + n_aug_));
  }

  // Predicted sigma points
  Xsig_pred_ = MatrixXd::Zero(n_x_, n_aug_ * 2 + 1);
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
  cout << "Entering ProcessMeasurement" << endl;
  
  if (!is_initialized_) {
    time_us_ = meas_package.timestamp_;
    
    x_ << 1, 1, 1, 1, 0.1;
    P_ << 0.15,    0, 0, 0, 0,
          0, 0.15, 0, 0, 0,
          0,    0, 1, 0, 0,
          0,    0, 0, 1, 0,
          0,    0, 0, 0, 1;
    
    is_initialized_ = true;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       Convert radar from polar to cartesian coordinates and initialize state.
       */
      auto rho = meas_package.raw_measurements_[0];
      auto phi = meas_package.raw_measurements_[1];
      x_(0) = cos(phi) * rho;
      x_(1) = sin(phi) * rho;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }
    return;
  }
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) {
    return;
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) {
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; // s from us
  cout << "Delta t: " << delta_t << endl;
  
  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  time_us_ = meas_package.timestamp_;
  cout << "Exiting ProcessMeasurement" << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  cout << "Entering Prediction" << endl;
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Generate augmented sigma points
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
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
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
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
    if (abs(yaw_dot) < 0.0001) { // avoid zero division
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
  
  // reset and predict state mean and state covariance matrix
  x_.fill(0.0);
  P_.fill(0.0);
  
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
  
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Tools::ConstrainAngle(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
  
  cout << "x" << endl << x_ << endl;
  cout << "P" << endl << P_ << endl;
  
  cout << "Exiting Prediction" << endl;
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
  cout << "Entering UpdateLidar" << endl;
  
  MatrixXd Zsig = MatrixXd::Zero(2, 2 * n_aug_ + 1);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }
  
  MatrixXd R = MatrixXd(2, 2);
  R << std_laspx_ * std_laspx_, 0.0,
       0.0, std_laspy_ * std_laspy_;
  
  Update(Zsig, R, meas_package.raw_measurements_);
  cout << "Exiting UpdateLidar" << endl;
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

  cout << "Entering UpdateRadar" << endl;
  
  //transform sigma points into measurement space
  //calculate mean predicted measurement
  //calculate measurement covariance matrix S
  
  MatrixXd Zsig = MatrixXd::Zero(3, 2 * n_aug_ + 1); // create matrix for sigma points in measurement space
  
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd sp = Xsig_pred_.col(i);
    float px = sp(0);
    float py = sp(1);
    float v = sp(2);
    float yaw = sp(3);
    //    float yawd = sp(4);
    
    VectorXd col = VectorXd(3);
    double sqrt_px_py = sqrt(px * px + py * py);
    double constr = Tools::ConstrainAngle(atan2(py, px));
    col << sqrt_px_py,
           constr,
           (px * cos(yaw) * v + py * sin(yaw) * v) / sqrt_px_py;
    Zsig.col(i) = col;
  }
  
  MatrixXd R = MatrixXd(3, 3);
  R << std_radr_ * std_radr_, 0.0, 0.0,
       0.0, std_radphi_ * std_radphi_, 0.0,
       0.0, 0.0, std_radrd_ * std_radrd_;
  
  // initialize predicted measurement covariance matrix S
  Update(Zsig, R, meas_package.raw_measurements_);
  cout << "Exiting UpdateRadar" << endl;
}

void UKF::Update(MatrixXd Zsig, MatrixXd R, VectorXd raw_measurements) {
  int n_z = Zsig.rows();
  
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i = 0; i < Zsig.cols(); i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  
  for (int i = 0; i < Zsig.cols(); i++) {
    MatrixXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = Tools::ConstrainAngle(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S += R;
  
  cout << "Zsig" << endl << Zsig << endl;
  cout << "z_pred" << endl << z_pred << endl;
  cout << "S measurement covariance" << endl << S << endl;
  
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    z_diff(1) = Tools::ConstrainAngle(z_diff(1));
    x_diff(3) = Tools::ConstrainAngle(x_diff(3));
    
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  cout << "S measurement covariance inverse" << endl << S.inverse() << endl;
  cout << "Tc" << endl << Tc << endl;
  
  MatrixXd K = Tc * S.inverse();
  cout << "K" << endl << K << endl;
  
  VectorXd z_diff = raw_measurements - z_pred;
  z_diff(1) = Tools::ConstrainAngle(z_diff(1));
  
  cout << "z_diff" << endl << z_diff << endl;
  
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  
  cout << "x" << endl << x_ << endl;
  cout << "P" << endl << P_ << endl;
}
