#pragma once
#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();
  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuation
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
  void TransformMapCoordinates(vector<double> &x_vals, vector<double> &y_vals,
                               double px, double py, double psi);

  // For converting back and forth between radians and degrees.
  const double pi();
  double deg2rad(double x);
  double rad2deg(double x);

  vector<double> x_vals;
  vector<double> y_vals;

  double steer_value;
  double throttle_value;

  const double angel = 18.9; // 25; 21; 22; 21.5; 19; 15;
  const double limit = std::numeric_limits<double>::max(); //1.0e19;
};
