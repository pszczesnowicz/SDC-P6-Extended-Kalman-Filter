#include "kalman_filter.h"
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

KalmanFilter::KalmanFilter() {}


KalmanFilter::~KalmanFilter() {}


void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  
}


void KalmanFilter::Predict() {
  
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
  
}


void KalmanFilter::Update(const VectorXd &z) {
  
  VectorXd y = z - H_ * x_;
  MatrixXd PHt = P_ * H_.transpose();
  MatrixXd S = H_ *  PHt + R_;
  MatrixXd K = PHt * S.inverse();
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
  
}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  VectorXd h(3);
  
  // Rho
  h(0) = sqrt(pow(x_(0), 2.0) + pow(x_(1), 2.0));
  
  // Phi
  if (x_(1) == 0 && x_(0) == 0) {
    h(1) = 0;
  }
  else {
    h(1) = atan2(x_(1), x_(0));
  }

  // Rho dot
  if (fabs(h(0)) < 0.0001) {
    h(2) = 0;
  }
  else {
    h(2) = (x_(0) * x_(2) + x_(1) * x_(3)) / h(0);
  }

  VectorXd y = z - h;
  
  if (y(1) > M_PI) {
    y(1) -= 2.0 * M_PI;
  }
  else if (y(1) < -M_PI) {
    y(1) += 2.0 * M_PI;
  }
  
  MatrixXd PHt = P_ * H_.transpose();
  MatrixXd S = H_ *  PHt + R_;
  MatrixXd K = PHt * S.inverse();
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
  
}
