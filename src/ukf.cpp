#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define E_P_S 0.001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // state vector dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  lambda_ = 3 - n_x_;

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  ///* Sigma points dimension
  n_sig_ = 2 * n_aug_ + 1;

  // Initialize weights.
  weights_ = VectorXd(n_sig_);
  weights_.fill(0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // Initialize measurement noice covarieance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;

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
  if ( !is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0]; // Range
      double phi = meas_package.raw_measurements_[1]; // Bearing
      double rho_dot = meas_package.raw_measurements_[2]; // velocity of rho
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
  	  double vy = rho_dot * sin(phi);
      double v_f = sqrt(vx * vx + vy * vy);
      x_ << x, y, v_f, 0, 0;
    } else {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // Saving first timestamp in seconds
    time_us_ = meas_package.timestamp_ ;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // Calculate dt
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  // Prediction step
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
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
  
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // Creating sigma points.
  MatrixXd Xsig_aug = GenerateSigmaPoints(x_aug, P_aug, lambda_, n_sig_);
  // 2. Predict Sigma Points.
  Xsig_pred_ = PredictSigmaPoints(Xsig_aug, delta_t, n_x_, n_sig_, std_a_, std_yawdd_);
  // 3. Predict Mean and Covariance
  //predicted state mean
  x_ = Xsig_pred_ * weights_;

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    NormalizeAngleOnComponent(x_diff, 3);

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
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

  // 1. Predit measurement
  int n_zz = 2;
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_zz, n_sig_);

  //Mean predicted measure
  VectorXd z_pred = VectorXd(n_zz);
  z_pred.fill(0.0);
  for (int k=0; k < n_sig_; k++) {
      z_pred = z_pred + weights_(k) * Zsig.col(k);
  }

  //Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_zz,n_zz);
  S.fill(0.0);
  for (int k = 0; k < n_sig_; k++) {  //2n+1 simga points
    //residual
    VectorXd z_dif = Zsig.col(k) - z_pred;

    S = S + weights_(k) * z_dif * z_dif.transpose();
  }

  
  //Add measure noise covariance matrix
  S = S + R_lidar_;

  // 2. Update state
  // Incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_zz);

  Tc.fill(0.0);
  for (int k = 0; k < n_sig_; k++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(k) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(k) - x_;

    Tc = Tc + weights_(k) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K_M = Tc * S.inverse();

  //residual
  VectorXd z_dif = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K_M * z_dif;
  P_ = P_ - K_M*S*K_M.transpose();

  //NIS Lidar Update
  NIS_laser_ = z_dif.transpose() * S.inverse() * z_dif;
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
  // Radar measument dimension
  int n_zz = 3;
  // 1. Predict measurement
  MatrixXd Zsig = MatrixXd(n_zz, n_sig_);
  //transform sigma points into measurement space
  for (int k = 0; k < n_sig_; k++) {  //2n+1 simga points

    double p_x = Xsig_pred_(0,k);
    double p_y = Xsig_pred_(1,k);
    double v  = Xsig_pred_(2,k);
    double yaw = Xsig_pred_(3,k);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // Measurement mod
    Zsig(0,k) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,k) = atan2(p_y,p_x);                                 //phi
    Zsig(2,k) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //Mean predicted measurement
  VectorXd z_pred = VectorXd(n_zz);
  z_pred.fill(0.0);
  for (int k=0; k < n_sig_; k++) {
      z_pred = z_pred + weights_(k) * Zsig.col(k);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_zz,n_zz);
  S.fill(0.0);
 
 
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
    //Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //Angle normalization
    NormalizeAngleOnComponent(z_diff, 1);

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  
  //Add measurement noise covariance matrix
  S = S + R_radar_;

  // Incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;

  //Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_zz);

  Tc.fill(0.0);
  for (int k = 0; k < n_sig_; k++) {  //2n+1 simga points

   
    VectorXd z_diff = Zsig.col(k) - z_pred;
    //Angle normalization
    NormalizeAngleOnComponent(z_diff, 1);

    // State difference
    VectorXd x_dif = Xsig_pred_.col(k) - x_;
    
	//Angle normalization
    NormalizeAngleOnComponent(x_dif, 3);

    Tc = Tc + weights_(k) * x_dif * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K_M = Tc * S.inverse();

  VectorXd z_dif = z - z_pred;

  //Angle normalization
  NormalizeAngleOnComponent(z_dif, 1);

  //update state mean and covariance matrix
  x_ = x_ + K_M * z_dif;
  P_ = P_ - K_M*S*K_M.transpose();

  //NIS Update
  NIS_radar_ = z_dif.transpose() * S.inverse() * z_dif;
}

/**
 *  Normalized the component `index` of the vector `vector` to be inside [-M_PI, M_PI] interval.
 */
void UKF::NormalizeAngleOnComponent(VectorXd vector, int index) {
  while (vector(index)> M_PI) vector(index)-=2.*M_PI;
  while (vector(index)<-M_PI) vector(index)+=2.*M_PI;
}

/**
 * Predits sigma points.
 **/
MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig, double delta_t, int n_x, int n_sig, double nu_am, double nu_yawdd) {
  MatrixXd Xsig_pred = MatrixXd(n_x, n_sig);
  //Predict sigma points
  for (int k = 0; k< n_sig; k++)
  {
    //extract values for better readability
    double px = Xsig(0,k);
    double py = Xsig(1,k);
    double v = Xsig(2,k);
    double yaw = Xsig(3,k);
    double yawd = Xsig(4,k);
    double nu_a = Xsig(5,k);
    double nu_yawdd = Xsig(6,k);

    //Predicted state values
    double px_p, py_p;

    //Avoid division by zero
    if (fabs(yawd) > E_P_S) {
        px_p = px + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = py + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = px + v*delta_t*cos(yaw);
        py_p = py + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //Add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //Write predicted sigma point into right column
    Xsig_pred(0,k) = px_p;
    Xsig_pred(1,k) = py_p;
    Xsig_pred(2,k) = v_p;
    Xsig_pred(3,k) = yaw_p;
    Xsig_pred(4,k) = yawd_p;
  }

  return Xsig_pred;
}

/**
 *   Generate sigma points:
 */
MatrixXd UKF::GenerateSigmaPoints(VectorXd x, MatrixXd P, double lambda, int n_sig) {
  int n = x.size();
  //Create sigma point matrix
  MatrixXd Xsig = MatrixXd( n, n_sig );

  //Calculate square root of P
  MatrixXd A = P.llt().matrixL();

  Xsig.col(0) = x;

  double lambda_plue_n_x_sqrt = sqrt(lambda + n);
  for (int k = 0; k < n; k++){
      Xsig.col( k + 1 ) = x + lambda_plue_n_x_sqrt * A.col(k);
      Xsig.col( k + 1 + n ) = x - lambda_plue_n_x_sqrt * A.col(k);
  }
  return Xsig;
}
