#pragma once
#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // For converting back and forth between radians and degrees.
  const double pi();
  double deg2rad(double x);
  double rad2deg(double x);
  const double angel = 19.0;

  double steer_value;
  double throttle_value;
  vector<double> x_pred;
  vector<double> y_pred;

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  void Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
  void TransformCoordinates(vector<double> &x_vals, vector<double> y_vals, double p_x, double p_y, double psi);
};
