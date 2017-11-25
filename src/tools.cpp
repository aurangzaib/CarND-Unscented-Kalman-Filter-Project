#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.empty() || estimations.size() != ground_truth.size()) return rmse;

  // accumulate squared residuals
  // todo: see if this loop works are expected
  // its implemented differently compared to udacity
  for (unsigned int i = 0; i < estimations.size(); ++i) {
    for (unsigned int j = 0; j < rmse.size(); j++) {
      auto residual = estimations[i][j] - ground_truth[i][j];
      rmse[j] += pow(residual, 2);
    }
  }

  // root mean residuals
  for (int i = 0; i < rmse.size(); i++) {

    // calculate the mean
    rmse[i] /= estimations.size();

    // calculate the squared root
    rmse[i] = sqrt(rmse[i]);
  }

  //return the result
  return rmse;
}