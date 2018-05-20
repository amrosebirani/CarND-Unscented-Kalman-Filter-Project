#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

static double normalize(double angle)
{
  while (angle > M_PI) {
    angle -= 2.*M_PI;
  };
  while (angle < -M_PI) {
    angle += 2.*M_PI;
  };
  
  return angle;
}

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .4;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .4;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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

  is_initialized_ = false;

  previous_timestamp_ = 0;

  n_aug_ = 7;

  n_x_ = 5;

  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2*n_aug_+1);
  weights_[0] = lambda_/(lambda_ + n_aug_);
  for (int i=1; i<=n_aug_; i++) {
      weights_[i] = .5 / (lambda_ + n_aug_);
      weights_[i+n_aug_] = .5 / (lambda_ + n_aug_);
  }

  R_Radar_ = MatrixXd(3, 3);
  R_Radar_ << pow(std_radr_, 2), 0, 0,
       0, pow(std_radphi_, 2), 0,
       0, 0, pow(std_radrd_, 2);

  R_Laser_ = MatrixXd(2, 2);
  R_Laser_ << pow(std_laspx_, 2), 0,
       0, pow(std_laspy_, 2);
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

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double x = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
      double y = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);
      double vx = meas_package.raw_measurements_[2]*cos(meas_package.raw_measurements_[1]);
      double vy = meas_package.raw_measurements_[2]*sin(meas_package.raw_measurements_[1]);
      double v = sqrt(pow(vx, 2) + pow(vy, 2));
      x_ << x, y, v, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
     x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  cout << "prediction done" << endl;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_) {
      cout << "perform radar update" << endl;
      UpdateRadar(meas_package);
    }
  } else {
    if (use_laser_) {
      cout << "perform lidar update" << endl;
      UpdateLidar(meas_package);
    }
  }

  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  cout << "NIS_ = " << NIS_ << endl;
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
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  x_aug.head(5) = x_;
  x_aug[5] = 0;
  x_aug[6] = 0;
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.row(n_aug_-1) << 0,0,0,0,0,0,pow(std_yawdd_, 2);
  P_aug.row(n_aug_-2) << 0,0,0,0,0,pow(std_a_, 2),0;
  MatrixXd A = P_aug.llt().matrixL();
  Xsig_aug.col(0) = x_aug;
  MatrixXd sqm = MatrixXd(n_aug_, n_aug_);
  sqm = sqrt(3) * A;
  for (int i=1; i<=n_aug_;i++) {
    Xsig_aug.col(i) = x_aug + sqm.col(i-1);
    Xsig_aug.col(i+n_aug_) = x_aug - sqm.col(i-1);
  }
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  for (int i=0; i<=2*n_aug_; i++) {
      VectorXd det(n_x_);
      VectorXd stoc(n_x_);
      VectorXd coli = Xsig_aug.col(i);
      stoc << .5*pow(delta_t, 2)*cos(coli[3])*coli[5],
              .5*pow(delta_t, 2)*sin(coli[3])*coli[5],
              delta_t*coli[5],
              .5*pow(delta_t, 2)*coli[6],
              delta_t*coli[6];
      if ( fabs(coli[4]) > 0.001 ) {
          det << coli[2]/coli[4]*(sin(coli[3] + coli[4]*delta_t) - sin(coli[3])),
                 coli[2]/coli[4]*(-cos(coli[3] + coli[4]*delta_t) + cos(coli[3])),
                 0,
                 coli[4]*delta_t,
                 0;
      } else {
          det << coli[2]*cos(coli[3])*delta_t,
                 coli[2]*sin(coli[3])*delta_t,
                 0,
                 coli[4]*delta_t,
                 0;
      }
      VectorXd xk(5);
      xk << coli[0], coli[1], coli[2], coli[3], coli[4];
      Xsig_pred.col(i) = det + stoc + xk;
  }
  Xsig_pred_ = Xsig_pred;
  std::cout << "Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  
  for (int i=0; i<2*n_aug_+1; i++) {
      x = x + weights_[i]*Xsig_pred.col(i);
  }
  x_ = x;
  std::cout << "x_ prediction done" << std::endl;
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);
  
  for (int i=0; i<2*n_aug_+1; i++) {
    VectorXd xd = Xsig_pred.col(i) - x;
    xd[3] = normalize(xd[3]);
    P = P + weights_[i]*xd*xd.transpose();
  }

  P_ = P;
  std::cout << "P_ prediction done" << std::endl;
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
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i=0; i<=2*n_aug_; i++) {
    VectorXd coli = Xsig_pred_.col(i);
    VectorXd zk(n_z);
    zk << coli[0], coli[1];
    Zsig.col(i) = zk;
  }
  for (int i=0; i<=2*n_aug_; i++) {
    z_pred = z_pred + weights_[i]*Zsig.col(i);
  }
  for (int i=0; i<=2*n_aug_; i++) {
    VectorXd zdd = Zsig.col(i) - z_pred;
    S = S + weights_[i]*zdd * zdd.transpose();
  }
  S = S + R_Laser_;

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  VectorXd zd = z - z_pred;
  Tc.fill(0.0);
  for (int i=0; i<=2*n_aug_; i++) {
      Tc = Tc + weights_[i]*(Xsig_pred_.col(i) - x_)*(Zsig.col(i) - z_pred).transpose();
  }
  MatrixXd K = Tc*S.inverse();
  x_ = x_+ K * zd;
  
  P_ = P_ - K*S*K.transpose();
  NIS_ = zd.transpose() * S.inverse() * zd;
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
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i=0; i<=2*n_aug_; i++) {
    VectorXd coli = Xsig_pred_.col(i);
    VectorXd zk(n_z);
    if (fabs(coli[0]) > .0001 || fabs(coli[1]) > .0001) {
      zk << sqrt(pow(coli[0], 2) + pow(coli[1], 2)),
            atan2(coli[1], coli[0]),
            (coli[0]*cos(coli[3])*coli[2] + coli[1]*sin(coli[3])*coli[2])/sqrt(pow(coli[0], 2) + pow(coli[1], 2));
    }
    else {
      zk << 0.0, 0.0, 0.0;
    }
    Zsig.col(i) = zk;
  }
  for (int i=0; i<=2*n_aug_; i++) {
    z_pred = z_pred + weights_[i]*Zsig.col(i);
  }
  for (int i=0; i<=2*n_aug_; i++) {
    VectorXd zdd = Zsig.col(i) - z_pred;
    zdd[1] = normalize(zdd[1]);
    S = S + weights_[i]*zdd * zdd.transpose();
  }
  S = S + R_Radar_;

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  VectorXd zd = z - z_pred;
  zd[1] = normalize(zd[1]);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i=0; i<=2*n_aug_; i++) {
      VectorXd zdk = Zsig.col(i) - z_pred;
      VectorXd xd = Xsig_pred_.col(i) - x_;
      xd[3] = normalize(xd[3]);
      zdk[1] = normalize(zdk[1]);
      Tc = Tc + weights_[i]*xd* zdk.transpose();
  }
  MatrixXd K = Tc*S.inverse();
  x_ = x_+ K * zd;
  
  P_ = P_ - K*S*K.transpose();
  NIS_ = zd.transpose() * S.inverse() * zd;
}
