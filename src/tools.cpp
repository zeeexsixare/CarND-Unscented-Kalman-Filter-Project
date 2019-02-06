#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
//using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse.fill(0);

  int est_size = estimations.size();

  if (est_size != ground_truth.size() || est_size == 0) {
    //cout << "Invalid vectors passed to CalculateRMSE." << std::endl;
    return rmse;
  }
  else{
    for (int i = 0; i < est_size; ++i){
      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array() * residual.array();
      rmse += residual;
    }
    rmse = rmse / est_size;
    rmse = rmse.array().sqrt();
    return rmse;
    }
  }
