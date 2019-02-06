#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <unistd.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

  // Weights of sigma points
  //Eigen::VectorXd weights_;

  //state dimension
  n_x_ = 5;

  //augmented state dimension
  n_aug_ = 7;

  //Number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  //sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

     // Matrix to hold sigma points
     Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

     // Vector for weights
     weights_ = VectorXd(n_sig_);
     weights_.fill(0.5/(lambda_ + n_aug_));
     weights_(0) = lambda_/(lambda_ + n_aug_);
     cout << "weights_" << weights_ << endl;


     time_us_ = 0;

     NIS_radar_ = 0;
     NIS_laser_ = 0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  //cout << "start ProcessMeasurement" << endl;
   if(!is_initialized_){
     time_us_ = meas_package.timestamp_;

     if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
       P_ << 0.3, 0, 0, 0, 0,
             0, 0.3, 0, 0, 0,
             0, 0, 0.3, 0, 0,
             0, 0, 0, 0.3, 0,
             0, 0, 0, 0, 0.3;

       double range = meas_package.raw_measurements_(0);
       double bearing = meas_package.raw_measurements_(1);
       x_ << range * cos(bearing),
             range * sin(bearing),
             0,
             0,
             0;

     }
     else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
       P_ << 0.3, 0, 0, 0, 0,
             0, 0.3, 0, 0, 0,
             0, 0, 0.3, 0, 0,
             0, 0, 0, 0.3, 0,
             0, 0, 0, 0, 0.3;
       x_ << meas_package.raw_measurements_[0],
             meas_package.raw_measurements_[1],
             0,
             0,
             0;
     }
     is_initialized_ = true;
   }
   else{
     double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
     time_us_ = meas_package.timestamp_; // put outside the loop

     Prediction(dt);

     if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
       UpdateRadar(meas_package);
     }
     else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
       UpdateLidar(meas_package);
     }
   }
cout << "finished ProcessMeasurement" << endl;

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
  //cout << "start Prediction" << endl;
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0);
  x_aug.head(n_x_) = x_;
  //cout << "start Prediction1" << endl;
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  //cout << "start Prediction2" << endl;
  MatrixXd A = P_aug.llt().matrixL();
  //A = A * ;
  //cout << "start Prediction3" << endl;
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  //cout << "start Prediction3.1" << endl;
  Xsig_aug.fill(0);
  //cout << "start Prediction3.2" << endl;
  Xsig_aug.col(0) = x_aug;
  //cout << "start Prediction3.3" << endl;
  for(int i = 0; i < n_aug_; i++){
    //cout << "start Prediction3.4" << endl;
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    //cout << "start Prediction3.5" << endl;
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
    //cout << "start Prediction3.6" << endl;
    }
  //cout << "start Prediction4" << endl;
  for (int i = 0; i < n_sig_; i++){
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yaw_dot = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yaw_doubledot = Xsig_aug(6, i);
    double cos_yaw = cos(yaw);
    double sin_yaw = sin(yaw);

    double p_x_pred = p_x;
    double p_y_pred = p_y;
    double v_pred = v;
    double yaw_pred = yaw;
    double yaw_dot_pred = yaw_dot;
    //cout << "start Prediction5" << endl;
    if(fabs(yaw_dot) > 0.0001){
      double angle = yaw + yaw_dot * delta_t;
      double v_yaw_dot = v/yaw_dot;

      p_x_pred += v_yaw_dot * (sin(angle) - sin_yaw);
      p_y_pred += v_yaw_dot * (-cos(angle) + cos_yaw);
      yaw_pred += yaw_dot * delta_t;
      //cout << "start Prediction6" << endl;
      }
    else{
      //cout << "start Prediction6.5" << endl;
      double v_dt = v*delta_t;
      p_x_pred += v_dt * cos_yaw;
      p_y_pred += v_dt * sin_yaw;
    }
    //cout << "start Prediction7" << endl;
    double half_dt2 = 0.5 * delta_t * delta_t;

    p_x_pred += half_dt2 * nu_a * cos_yaw;
    p_y_pred += half_dt2 * nu_a * sin_yaw;
    v_pred += delta_t * nu_a;
    yaw_pred +=half_dt2 * nu_yaw_doubledot;
    yaw_dot_pred += delta_t * nu_yaw_doubledot;

    //cout << "start Prediction8" << endl;
    Xsig_pred_(0,i) = p_x_pred;
    Xsig_pred_(1,i) = p_y_pred;
    Xsig_pred_(2, i) = v_pred;
    Xsig_pred_(3, i) = yaw_pred;
    Xsig_pred_(4, i) = yaw_dot_pred;
  }
  //cout << "start Prediction9" << endl;
  x_.fill(0);
  //cout << "start Prediction9.1" << endl;
  for(int i = 0; i < n_sig_; i++){
    //cout << "start Prediction9.2" << endl;
    //cout << "weights_" << weights_ << endl;
    //cout << "Xsig_pred_" << Xsig_pred_ << endl;
    x_ += weights_(i) * Xsig_pred_.col(i);
    //cout << "start Prediction9.3" << endl;
  }
  cout << "Prediction x_" << endl << x_ << endl;
  //cout << "start Prediction10" << endl;
  P_.fill(0);
  for(int i=0; i < n_sig_; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    x_diff(3) = NormalizeAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
  cout << "Prediction P_" << endl << P_ << endl;
cout << "finished Prediction" << endl;
}



void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  * TODO: Complete this function! Use lidar data to update the belief
  * about the object's position. Modify the state vector, x_, and
  * covariance, P_.
  * You can also calculate the lidar NIS, if desired.
  */
  //usleep(500000);
  cout << "start UpdateLidar" << endl;

  MatrixXd Zsig = MatrixXd(2, n_sig_);
  Zsig.row(0) = Xsig_pred_.row(0);
  Zsig.row(1) = Xsig_pred_.row(1);

  VectorXd z_pred = VectorXd(2);
  z_pred.fill(0);
  for (int i = 0; i < n_sig_; i++){
    z_pred += weights_(i) * Zsig.col(i);
  }

  MatrixXd S = MatrixXd(2, 2);
  S.fill(0);
  for (int i = 0; i < n_sig_; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  S(0, 0) += std_laspx_*std_laspx_;
  S(1, 1) += std_laspy_*std_laspy_;

  MatrixXd T_c = MatrixXd(n_x_, 2);
  T_c.fill(0);
  for (int i = 0; i< n_sig_; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;

    T_c += weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd S_i = S.inverse();
  MatrixXd K = T_c*S_i;
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
  cout << "UpdateLidar x_" << endl << x_ << endl;
  cout << "UpdateLidar P_" << endl << P_ << endl;
  NIS_laser_ = z_diff.transpose() * S_i * z_diff;
  cout << "finished UpdateLidar" << endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  * TODO: Complete this function! Use radar data to update the belief
  * about the object's position. Modify the state vector, x_, and
  * covariance, P_.
  * You can also calculate the radar NIS, if desired.
  */
  //usleep(500000);
  cout << "start UpdateRadar" << endl;


  MatrixXd Zsig = MatrixXd(3, n_sig_);
  Zsig.fill(0);

  for(int i = 0; i < n_sig_; i++){
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double range = sqrt(p_x * p_x + p_y * p_y);
    Zsig(0,i) = range;
    Zsig(1, i) = atan2(p_y,p_x);
    Zsig(2,i) = (p_x * cos(yaw)*v + p_y * sin(yaw) * v) / range;
  }

  VectorXd z_pred = VectorXd(3);
  z_pred.fill(0);
  for(int i = 0; i < n_sig_; i++){
    z_pred += weights_(i) * Zsig.col(i);
  }

  MatrixXd S = MatrixXd(3, 3);
  S.fill(0);
  for (int i = 0; i < n_sig_; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  S(0,0) += std_radr_*std_radr_;
  S(1,1) += std_radphi_*std_radphi_;
  S(2,2) += std_radrd_*std_radrd_;

  MatrixXd Tc = MatrixXd(n_x_, 3);
  Tc.fill(0);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd Si = S.inverse();

  MatrixXd K = Tc * Si;

  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  z_diff(1) = NormalizeAngle(z_diff(1));

  x_ += K* z_diff;
  P_ -= K * S * K.transpose();

  cout << "UpdateRadar x_" << endl << x_ << endl;
  cout << "UpdateRadar P_" << endl << P_ << endl;
  NIS_radar_ = z_diff.transpose() * Si * z_diff;
  cout << "finished UpdateRadar" << endl;
}

double UKF::NormalizeAngle(double angle){
  //cout << "start NormalizeAngle" << endl;
  while(angle > M_PI || angle < -M_PI){
    if(angle > M_PI){
      angle -= 2 * M_PI;
    }
    else if(angle < -M_PI){
      angle += 2 * M_PI;
    }
  }
  return angle;
  //cout << "finished NormalizeAngle" << endl;
}
