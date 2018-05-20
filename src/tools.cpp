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

  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
	    cout << "Bad Input" << endl;
	    return rmse;
	}

  for(int i=0; i < estimations.size(); ++i){
    VectorXd diff = estimations[i] - ground_truth[i];
    VectorXd sqd = diff.array()*diff.array();
		rmse[0] = rmse[0] + sqd[0];
		rmse[1] = rmse[1] + sqd[1];
		rmse[2] = rmse[2] + sqd[2];
		rmse[3] = rmse[3] + sqd[3];
	}

  rmse = rmse/estimations.size();

  rmse = rmse.array().sqrt();

  return rmse;
}