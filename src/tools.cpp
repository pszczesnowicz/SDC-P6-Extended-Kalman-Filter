#include <iostream>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

// Constructor.
Tools::Tools() {}

// Destructor.
Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  VectorXd root_mean_squared_error(4);
  root_mean_squared_error << 0, 0, 0, 0;
  
  if(estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    cout << "CalculateRMSE() - Error: Estimation and ground truth vectors not equal in size.";
    return root_mean_squared_error;
  }
  
  for(int i = 0; i < estimations.size(); ++i) {
    VectorXd squared_residuals = pow((estimations[i] - ground_truth[i]).array(), 2.0);
    root_mean_squared_error += squared_residuals;
  }
  
  // Calculates the mean.
  root_mean_squared_error = root_mean_squared_error / estimations.size();
  
  // Calculates the squared root.
  root_mean_squared_error = root_mean_squared_error.cwiseSqrt();
  
  return root_mean_squared_error;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
  
  MatrixXd jacobian(3, 4);
  jacobian << 0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0;
  
  float position_x = x_state(0);
  float position_y = x_state(1);
  float velocity_x = x_state(2);
  float velocity_y = x_state(3);
  
  float divisor = pow(position_x, 2.0) + pow(position_y, 2.0);
  
  if(fabs(divisor) < 0.0001) {
    cout << "CalculateJacobian() - Error: Division by zero." << endl;
    return jacobian;
  }
  
  float divisor_root = sqrt(divisor);
  float divisor_three_halves = divisor * divisor_root;
  
  jacobian << (position_x / divisor_root), (position_y / divisor_root), 0, 0,
              (- position_y / divisor), (position_x / divisor), 0, 0,
              ((position_y * (velocity_x * position_y - velocity_y * position_x)) / divisor_three_halves), ((position_x * (velocity_y * position_x - velocity_x * position_y)) / divisor_three_halves), (position_x / divisor_root), (position_y / divisor_root);
  
  return jacobian;
}
