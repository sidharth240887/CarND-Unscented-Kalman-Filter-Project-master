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
  rmse << 0,0,0,0;

  if(estimations.size() == 0){
    cout << "ERROR CalculateRMSE () - The estimations vec is empty" << endl;
    return rmse;
  }

  if(ground_truth.size() == 0){
    cout << "ERROR CalculateRMSE () - The ground truth vec is empty" << endl;
    return rmse;
  }

  unsigned int cap = estimations.size();
  if(cap != ground_truth.size()){
    cout << "ERROR CalculateRMSE () - The ground-truth and estimations vectors must have same size." << endl;
    return rmse;
  }

  for(unsigned int i=0; i < estimations.size(); ++i){
    VectorXd dif = estimations[i] - ground_truth[i];
    dif = dif.array()*dif.array();
    rmse += dif;
  }

  rmse = rmse / cap;
  rmse = rmse.array().sqrt();
  return rmse;
}