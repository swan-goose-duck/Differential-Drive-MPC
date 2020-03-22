#include "MPC.h"

using CppAD::AD;
using namespace std;

class FG_eval {

  public:

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    Eigen::VectorXd coeffs; // Fitted road curve polynomial coefficients

    FG_eval(Eigen::VectorXd Kin) {
      this->coeffs = coeffs; 
    }
  
    AD<double> cost_cte, cost_etheta, cost_vel; // costs

    void operator()(ADvector& fg, const ADvector& vars) {
      // fg a vector containing the cost and all constraints
      // vars is a vector containing all states and actuations for _mpc_steps "lookahead" states and actuations.

      //*********************************************************
      //* COST DEFINED HERE
      //*********************************************************
      fg[0] = 0;
      cost_cte =  0;
      cost_etheta = 0;
      cost_vel = 0;
      fg[0] = 0.0;

      for (int i = 0; i < _mpc_steps; i++) {

        fg[0] += _w_cte    * CppAD::pow(vars[_cte_start + i]    - _ref_cte, 2);    // cross deviation error
        fg[0] += _w_etheta * CppAD::pow(vars[_etheta_start + i] - _ref_etheta, 2); // heading error
        fg[0] += _w_vel    * CppAD::pow(vars[_v_start + i]      - _ref_vel, 2);    // speed error

        // for debugging
        cost_cte    +=  _w_cte    * CppAD::pow(vars[_cte_start + i]    - _ref_cte, 2);
        cost_etheta +=  _w_etheta * CppAD::pow(vars[_etheta_start + i] - _ref_etheta, 2); 
        cost_vel    +=  _w_vel    * CppAD::pow(vars[_v_start + i]      - _ref_vel, 2); 
      }

      for (int i = 0; i < _mpc_steps - 1; i++) {
        fg[0] += _w_angvel * CppAD::pow(vars[_angvel_start + i], 2);
        fg[0] += _w_accel  * CppAD::pow(vars[_a_start + i], 2);
      }

      for (int i = 0; i < _mpc_steps - 2; i++) {
        fg[0] += _w_angvel_d * CppAD::pow(vars[_angvel_start + i + 1] - vars[_angvel_start + i], 2);
        fg[0] += _w_accel_d  * CppAD::pow(vars[_a_start + i + 1]      - vars[_a_start + i],      2);
      }

      //*********************************************************
      //* CONSTRAINTS DEFINED HERE
      //*********************************************************

      // given state does not vary
      fg[1 + _x_start]      = vars[_x_start];
      fg[1 + _y_start]      = vars[_y_start];
      fg[1 + _theta_start]  = vars[_theta_start];
      fg[1 + _v_start]      = vars[_v_start];
      fg[1 + _cte_start]    = vars[_cte_start];
      fg[1 + _etheta_start] = vars[_etheta_start];

      // constraints based on our kinematic model
      for (int i = 0; i < _mpc_steps - 1; ++i) {
        
        // The state at time t+1 .
        AD<double> x1      = vars[_x_start      + i + 1];
        AD<double> y1      = vars[_y_start      + i + 1];
        AD<double> theta1  = vars[_theta_start  + i + 1];
        AD<double> v1      = vars[_v_start      + i + 1];
        AD<double> cte1    = vars[_cte_start    + i + 1];
        AD<double> etheta1 = vars[_etheta_start + i + 1];
        
        // The state at time t.
        AD<double> x0      = vars[_x_start      + i];
        AD<double> y0      = vars[_y_start      + i];
        AD<double> theta0  = vars[_theta_start  + i];
        AD<double> v0      = vars[_v_start      + i];
        AD<double> cte0    = vars[_cte_start    + i];
        AD<double> etheta0 = vars[_etheta_start + i];
        
        // Only consider the actuation at time t.
        AD<double> w0 = vars[_angvel_start + i]; // angular velocity
        AD<double> a0 = vars[_a_start + i];      // acceleration

        // desired py and psi
        AD<double> f0 = 0.0;
        for (int i = 0; i < coeffs.size(); i++) {
            f0 += coeffs[i] * CppAD::pow(x0, i); // f(x0)
        }
        
        AD<double> trj_grad0 = 0.0;
        for (int i = 1; i < coeffs.size(); i++) {
            trj_grad0 += i*coeffs[i] * CppAD::pow(x0, i-1); // f'(x0)
        }
        trj_grad0 = CppAD::atan(trj_grad0);

        // relationship of current state + actuations and next state
        // based on our kinematic model
        const auto xf      = x0 + v0 * CppAD::cos(theta0) * _dt;
        const auto yf      = y0 + v0 * CppAD::sin(theta0) * _dt;
        const auto thetaf  = theta0 +  w0 * _dt;
        const auto vf      = v0 + a0 * _dt;
        const auto ctef    = (f0 - y0) + (v0 * CppAD::sin(etheta0) * _dt);
        const auto ethetaf = (theta0 - trj_grad0) + w0 * _dt;

        // store the constraint expression of two consecutive states
        fg[2 + _x_start      + i] = x1 - xf;
        fg[2 + _y_start      + i] = y1 - yf;
        fg[2 + _theta_start  + i] = theta1 - thetaf;
        fg[2 + _v_start      + i] = v1 - vf;
        fg[2 + _cte_start    + i] = cte1 - ctef;
        fg[2 + _etheta_start + i] = etheta1 - ethetaf;
      }
    }
};

MPC::MPC() {

  //**************************************************************
  //* SET INITIAL VALUES OF VARIABLES
  //**************************************************************
  this->vars.resize(n_vars);

  // the aformentioned states will be initialized when solve() is called

  for (int i = 0; i < n_vars; ++i) {
    this->vars[i] = 0.0;
  }

  //**************************************************************
  //* SET UPPER AND LOWER LIMITS OF VARIABLES
  //**************************************************************

  this->vars_lowerbound.resize(n_vars);
  this->vars_upperbound.resize(n_vars);

    // Set all non-actuators upper and lowerlimits
    // to the max negative and positive values.
    for (int i = 0; i < _angvel_start; i++) {
      vars_lowerbound[i] = -_bound_value;
      vars_upperbound[i] =  _bound_value;
    }

    // The upper and lower limits of angvel are set to -25 and 25
    // degrees (values in radians).
    for (int i = _angvel_start; i < _a_start; i++) {
      vars_lowerbound[i] = -_max_angvel;
      vars_upperbound[i] =  _max_angvel;
    }

    // Acceleration/decceleration upper and lower limits
    for (int i = _a_start; i < n_vars; i++) {
      vars_lowerbound[i] = -_max_accel;
      vars_upperbound[i] =  _max_accel;
    }

  //**************************************************************
  //* SET UPPER AND LOWER LIMITS OF CONSTRAINTS
  //**************************************************************
  this->constraints_lowerbound.resize(n_constraints);
  this->constraints_upperbound.resize(n_constraints);

  // the first constraint for each state veriable
  // refer to the initial state conditions
  // this will be initialized when solve() is called
  // the succeeding constraints refer to the relationship
  // between succeeding states based on our kinematic model of the system

  for (int i = 0; i < n_constraints; ++i) {
    this->constraints_lowerbound[i] = 0.0;
    this->constraints_upperbound[i] = 0.0;
  }
}

MPC::~MPC() {}

vector<double> MPC::solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {

  typedef CPPAD_TESTVECTOR(double) Dvector;
  const double x      = state[0];
  const double y      = state[1];
  const double theta  = state[2];
  const double v      = state[3];
  const double cte    = state[4];
  const double etheta = state[5];

  vars[_x_start]      = x;
  vars[_y_start]      = y;
  vars[_theta_start]  = theta;
  vars[_v_start]      = v;
  vars[_cte_start]    = cte;
  vars[_etheta_start] = etheta;

  constraints_lowerbound[_x_start]      = x;
  constraints_lowerbound[_y_start]      = y;
  constraints_lowerbound[_theta_start]  = theta;
  constraints_lowerbound[_v_start]      = v;
  constraints_lowerbound[_cte_start]    = cte;
  constraints_lowerbound[_etheta_start] = etheta;

  constraints_upperbound[_x_start]      = x;
  constraints_upperbound[_y_start]      = y;
  constraints_upperbound[_theta_start]  = theta;
  constraints_upperbound[_v_start]      = v;
  constraints_upperbound[_cte_start]    = cte;
  constraints_upperbound[_etheta_start] = etheta;

  //**************************************************************
  //* SOLVE
  //**************************************************************

  // object that computes objective and constraints
  FG_eval fconstraints_eval(coeffs);

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
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options,
      vars,
      vars_lowerbound,
      vars_upperbound,
      constraints_lowerbound,
      constraints_upperbound,
      fconstraints_eval,
      solution);

  // comment out the lines below to debug!
  /*
  bool ok = true;
  auto cost = solution.obj_value;

  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  if (ok) {
    std::cout << "OK! Cost:" << cost << std::endl;
  } else {
    std::cout << "SOMETHIn_constraints IS WROn_constraints!" << cost << std::endl;
  }
  */

  //**************************************************************
  //* STORE RELEVANT INFORMATION FROM SOLUTION
  //**************************************************************

  this->future_xs = {};
  this->future_ys = {};

  for (int i = 0; i < _mpc_steps; ++i) {
    this->future_xs.push_back(solution.x[_x_start + i]);
    this->future_ys.push_back(solution.x[_y_start + i]);
  }

  this->acceleration     = solution.x[_a_start];
  this->angular_velocity = solution.x[_angvel_start];
  
  vector<double> speeds;
  double desired_velocity = v + acceleration;
  
  double vl = (2 * desired_velocity - angular_velocity * L) / (2 * R);
  double vr = (2 * desired_velocity + angular_velocity * L) / (2 * R);
  
  speeds[0] = vl;
  speeds[1] = vr;

  return speeds;
}
