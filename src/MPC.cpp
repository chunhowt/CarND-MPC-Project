#include <limits>

#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

size_t N = 10;
double dt = 0.1;  // 100 ms.

// FG_eval takes in a single vector of all the states and actuator values over
// all the time period. Thus, we flatten these values into a 1-D vector and
// here, we use compute the start index of each of the group of value.
// Note, in general, each state value happens N times and actuator values
// happens for N - 1 time period.
const int X_START = 0;
const int Y_START = X_START + N;
const int PHI_START = Y_START + N;
const int V_START = PHI_START + N;
const int CTE_START = V_START + N;
const int EPSI_START = CTE_START + N;
const int STEERING_START = EPSI_START + N;
const int THROTTLE_START = STEERING_START + N - 1;

// Corresponding to 25 degrees, as mentioned for this project.
const double MAX_STEERING_RADS = 0.44;

const double IDEAL_SPEED = 50;

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

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  // fg a vector of cost + constraints, vars is the vectors containing states and actuators.
  void operator()(ADvector& fg, const ADvector& vars) {
    // First, set up cost, fg[0].
    fg[0] = 0;

    // error compared to ideal trajectory.
    for (int i = 0; i < N; ++i) {
      // Minimize cross track error.
      fg[0] += CppAD::pow(vars[CTE_START + i], 2);
      // Minimize orientation error.
      fg[0] += CppAD::pow(vars[EPSI_START + i], 2);
      // Make sure the vehicle keep moving even if there is no error left.
      fg[0] += CppAD::pow(vars[V_START + i] - IDEAL_SPEED, 2); 
    }

    // Minimize the change in actuator value to smooth the driving.
    for (int i = 0; i < N - 2; ++i) {
      fg[0] += CppAD::pow(vars[STEERING_START + i + 1] - vars[STEERING_START + i], 2);
      fg[0] += CppAD::pow(vars[THROTTLE_START + i + 1] - vars[THROTTLE_START + i], 2);
    }

    // Set up constraints.
    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + X_START] = vars[X_START];
    fg[1 + Y_START] = vars[Y_START];
    fg[1 + PHI_START] = vars[PHI_START];
    fg[1 + V_START] = vars[V_START];
    fg[1 + CTE_START] = vars[CTE_START];
    fg[1 + EPSI_START] = vars[EPSI_START];

    // The rest of the constraints
    for (int i = 0; i < N - 1; i++) {
      // The state at time t+1 .
      AD<double> x1 = vars[X_START + i + 1];
      AD<double> y1 = vars[Y_START + i + 1];
      AD<double> psi1 = vars[PHI_START + i + 1];
      AD<double> v1 = vars[V_START + i + 1];
      AD<double> cte1 = vars[CTE_START + i + 1];
      AD<double> epsi1 = vars[EPSI_START + i + 1];

      // The state at time t.
      AD<double> x0 = vars[X_START + i];
      AD<double> y0 = vars[Y_START + i];
      AD<double> psi0 = vars[PHI_START + i];
      AD<double> v0 = vars[V_START + i];
      AD<double> cte0 = vars[CTE_START + i];
      AD<double> epsi0 = vars[EPSI_START + i];

      // Get actuator at time t.
      AD<double> curr_steering = vars[STEERING_START + i];
      AD<double> curr_throttle = vars[THROTTLE_START + i];

      // We use polynomial of degree 3.
      AD<double> cte_f = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      AD<double> epsi_f = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0 * x0);

      // Now, we use the kinematic vehicle model.
      fg[1 + X_START + i + 1] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + Y_START + i + 1] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + PHI_START + i + 1] = psi1 - (psi0 + v0 * curr_steering / Lf * dt);
      fg[1 + V_START + i + 1] = v1 - (v0 + curr_throttle * dt);

      // and the cross track / orientation error.
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[1 + CTE_START + i + 1] = cte1 - ((cte_f - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + EPSI_START + i + 1] = epsi1 - ((psi0 - epsi_f) + v0 * curr_steering / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs,
                          vector<double>* mpc_xs, vector<double>* mpc_ys) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // N timesteps of states, N - 1 timesteps of actuators.
  size_t n_vars = state.size() * N + 2 * (N - 1);
  // N timesteps of states.
  size_t n_constraints = N * state.size();

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  vars[X_START] = state[0];
  vars[Y_START] = state[1];
  vars[PHI_START] = state[2];
  vars[V_START] = state[3];
  vars[CTE_START] = state[4];
  vars[EPSI_START] = state[5];

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // Set lower and upper limits for variables.
  // For state values, it can be +/- infinity.
  for (int i = 0; i < STEERING_START; ++i) {
    vars_lowerbound[i] = std::numeric_limits<double>::min();
    vars_upperbound[i] = std::numeric_limits<double>::max();
  }
  // Set the limit for steering actuator.
  for (int i = STEERING_START; i < THROTTLE_START; ++i) {
    vars_lowerbound[i] = - MAX_STEERING_RADS;
    vars_upperbound[i] = MAX_STEERING_RADS;
  }
  // Set the limit for throttle actuator.
  for (int i = THROTTLE_START; i <  n_vars; ++i) {
    vars_lowerbound[i] = -1;
    vars_upperbound[i] = 1;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  // Set initial state as constraints.
  constraints_lowerbound[X_START] = constraints_upperbound[X_START] = state[0];
  constraints_lowerbound[Y_START] = constraints_upperbound[Y_START] = state[1];
  constraints_lowerbound[PHI_START] = constraints_upperbound[PHI_START] = state[2];
  constraints_lowerbound[V_START] = constraints_upperbound[V_START] = state[3];
  constraints_lowerbound[CTE_START] = constraints_upperbound[CTE_START] = state[4];
  constraints_lowerbound[EPSI_START] = constraints_upperbound[EPSI_START] = state[5];

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

  // Save the predicted trajectory and pass them out.
  for (int i = 0; i < N; ++i) {
    mpc_xs->push_back(solution.x[X_START + i]);
    mpc_ys->push_back(solution.x[Y_START + i]);
  }
  return {solution.x[STEERING_START + 1], solution.x[THROTTLE_START + 1]};
}
