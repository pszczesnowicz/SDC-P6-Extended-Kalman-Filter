#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


// Class constructor.
FusionEKF::FusionEKF() {
  
  is_initialized_ = false;

  previous_timestamp_ = 0;
  
  noise_ax_ = 9;
  noise_ay_ = 9;
  
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  // Measurement noise covariance matrix - laser.
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // Measurement noise covariance matrix - radar.
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  // Measurement matrix - laser.
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
}


// Class destructor.
FusionEKF::~FusionEKF() {}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  
  if (!is_initialized_) {
    
    cout << "EKF: " << endl;
    
    ekf_.x_ = VectorXd(4);

    // State vector x & y position initialization using the first RADAR measurements of rho and phi.
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      ekf_.x_(0) = measurement_pack.raw_measurements_(0) * cos(measurement_pack.raw_measurements_(1));
      ekf_.x_(1) = measurement_pack.raw_measurements_(0) * sin(measurement_pack.raw_measurements_(1));
      ekf_.x_(2) = 0;
      ekf_.x_(3) = 0;
    }
    
    // State vector x & y position initialization using the first LASER measurements of x and y.
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
      ekf_.x_(2) = 0;
      ekf_.x_(3) = 0;
    }
    
    // Uncertainty covariance matrix initialization.
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ <<  0.1, 0, 0, 0,
                0, 0.1, 0, 0,
                0, 0, 1000, 0,
                0, 0, 0, 1000;
    
    // State transition matrix initialization.
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ <<  1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0,
                0, 0, 0, 1;
    
    // Process noise covariance matrix initialization.
    ekf_.Q_ = MatrixXd(4, 4);
    
    previous_timestamp_ = measurement_pack.timestamp_;

    is_initialized_ = true;
    
    return;
  }

  // Prediction step.
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  if (dt > 0) {
    // State transition matrix update.
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    
    // Process noise covariance matrix update.
    float dt_squared = pow(dt, 2.0);
    float dt_cubed = pow(dt, 3.0) / 2.0;
    float dt_quad = pow(dt, 4.0) / 4.0;
    
    ekf_.Q_ <<  (dt_quad * noise_ax_), 0, (dt_cubed * noise_ax_), 0,
    0, (dt_quad * noise_ay_), 0, (dt_cubed * noise_ay_),
    (dt_cubed * noise_ax_), 0, (dt_squared * noise_ax_), 0,
    0, (dt_cubed * noise_ay_), 0, (dt_squared * noise_ay_);
    
    ekf_.Predict();
  }

  // Update step
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
      
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
