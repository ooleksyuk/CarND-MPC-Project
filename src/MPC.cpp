#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Set the time step length and duration
size_t N = 10;//25;
double dt = 0.1;//0.05;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

const double angel = 19.0;

// Set up reference values for cte, epsi and velocity
double ref_v = 60.0; //mph
double ref_cte = 0.0;
double ref_epsi = 0.0;

// Define and initialize vector X variable and actuator
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // Implement MPC
    fg[0] = 0;
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    // Reference 'state cost'. Define cost related reference state and
    // anything related to it;
    for (int t = 0; t < N; t++) {
      fg[0] += 3000 * CppAD::pow(vars[cte_start + t] - ref_cte, 2);
      fg[0] += 3000 * CppAD::pow(vars[epsi_start + t] - ref_epsi, 2);
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // Punish actuations, minimize it's use.
    for (int t = 0; t < N - 1; t++) {
      fg[0] += 500 * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += 5 * CppAD::pow(vars[a_start], 2);
      fg[0] += 5 * CppAD::pow(vars[delta_start], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < N - 2; t++) {
      fg[0] += 200 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += 10 * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    // Set up model constrains.
    // Initial constrains. Add 1 to each of the starting indices due to cost being located at
    // index 0 of 'fg'.
    // This bumps up the position of all other values.
    fg[x_start + 1]    = vars[x_start];
    fg[y_start + 1]    = vars[y_start];
    fg[psi_start + 1]  = vars[psi_start];
    fg[v_start + 1]    = vars[v_start];
    fg[cte_start + 1]  = vars[cte_start];
    fg[epsi_start + 1] = vars[epsi_start];

    for (int t = 0; t < N - 1; t++) {
      // Set up x1, y1, psi1, v1, cte1, epsi1 at time 't + 1'.
      AD<double> x1    = vars[x_start + t + 1];
      AD<double> y1    = vars[y_start + t + 1];
      AD<double> psi1  = vars[psi_start + t + 1];
      AD<double> v1    = vars[v_start + t + 1];
      AD<double> cte1  = vars[cte_start + t + 1];
      AD<double> epsi1 = vars[epsi_start + t + 1];

      // Set up x0, y0, psi0, v0, cte0, epsi0 at time 't'.
      AD<double> x0    = vars[x_start + t];
      AD<double> y0    = vars[y_start + t];
      AD<double> psi0  = vars[psi_start + t];
      AD<double> v0    = vars[v_start + t];
      AD<double> cte0  = vars[cte_start + t];
      AD<double> epsi0 = vars[epsi_start + t];

      // Only consider the actuation at time 't'.
      AD<double> delta0 = vars[delta_start + t];
      AD<double> a0 = vars[a_start + t];

      if (t > 1) {
        delta0 = vars[delta_start + t - 1];
        a0 = vars[a_start + t - 1];
      }

      // Define and initialize 3rd order derivative
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 * coeffs[3] * CppAD::pow(x0, 2));

      // Constraint this value to be 0.
      // Equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[x_start + t + 2]    = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[y_start + t + 2]    = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[psi_start + t + 2]  = psi1 - (psi0 + v0 * delta0/ Lf * dt);
      fg[v_start + t + 2]    = v1 - (v0 + a0 * dt);
      fg[cte_start + t + 2]  = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[epsi_start + t + 2] = epsi1 - (psi0 - psides0 + v0  * delta0/ Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

const double MPC::pi() { return M_PI; }
double MPC::deg2rad(double x) { return x * pi() / 180; }
double MPC::rad2deg(double x) { return x * 180 / pi(); }

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  x_vals.clear();
  y_vals.clear();
//  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];
  size_t n_vars = N * 6 + (N - 1) * 2;
  // Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // Set lower and upper limits for variables.
  // Set all non-actuators upper and lower limits
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19; // -std::numeric_limits<double>::min();
    vars_upperbound[i] = 1.0e19; // std::numeric_limits<double>::max();
  }

  // The upper and lower limits of delta are set to -25 and 25 degrees (values in radians).
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = - deg2rad(angel);
    vars_upperbound[i] = deg2rad(angel);
  }

  // Acceleration/decceleration upper and lower limits.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  std::cout << "Delta " << solution.x[delta_start + 2] << std::endl;
  std::cout << "Error psi " << solution.x[epsi_start + 1] << std::endl;
  std::cout << "Accuracy " << solution.x[a_start + 2] << std:: endl;

  for(int i = 0; i < N - 1; i++) {
    x_vals.push_back(solution.x[x_start + i + 1]);
    y_vals.push_back(solution.x[y_start + i + 1]);
  }

  return {
      solution.x[x_start + 1],
      solution.x[y_start + 1],
      solution.x[psi_start + 1],
      solution.x[v_start + 1],
      solution.x[cte_start + 1],
      solution.x[epsi_start + 1],
      (solution.x[delta_start + 1] + solution.x[delta_start]) / deg2rad(angel),
      solution.x[a_start + 1] + solution.x[a_start]
  };
}

void MPC::TransformCoordinates(vector<double> &x_vals, vector<double> y_vals,
                               double p_x, double p_y, double psi) {
  vector<double> t_x;
  vector<double> t_y;

  for(int i = 0; i < x_vals.size(); i ++) {
    double x_diff = x_vals[i] - p_x;
    double y_diff = y_vals[i] - p_y;
    t_x.push_back(x_diff * cos(-psi) + y_diff * sin(-psi));
    t_y.push_back(x_diff * sin(-psi) + y_diff * cos(-psi));
  }

  x_vals = t_x;
  y_vals = t_y;

  return;
}