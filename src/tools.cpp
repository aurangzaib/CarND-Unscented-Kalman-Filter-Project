#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  // Check validity
  if (estimations.size() <= 0) {
    std::cout << "CalculateRMSE - Division by zero ERROR" << std::endl;
    return rmse;
  } else if (estimations.size() != ground_truth.size()) {
    std::cout << "CalculateRMSE - Estimation size and ground_truth size differ ERROR" << std::endl;
    return rmse;
  }

  // Accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}