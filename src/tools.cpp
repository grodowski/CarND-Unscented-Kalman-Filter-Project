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
  rmse << 0,0,0,0;
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() == 0) {
    throw invalid_argument("received 0 length vector");
  }
  if (estimations.size() != ground_truth.size()) {
    throw invalid_argument("estimations and ground_truth does not match");
  }
  
  
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }
  //calculate the mean
  rmse = rmse / estimations.size();
  return rmse.array().sqrt();
}
