#ifndef MPC_H
#define MPC_H

#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <Eigen/Core>
#include <vector>

typedef CPPAD_TESTVECTOR(double) Dvector;

const int _mpc_steps = 20;     // how many states we "lookahead" in the future
const double _dt     = 0.02;   // how much time we expect environment changes
const double L       = 2.67;   // this is the length from front of vehicle to COM
const double R       = 0.0762; // radius of each wheel in meters

// desired values for errors and velocity
const double _ref_cte    = 0;
const double _ref_etheta = 0;
const double _ref_vel    = 0.5; // m/s

const int NUMBER_OF_STATES     = 6; // x, y, theta, v, cte, etheta (error in heading)
const int NUMBER_OF_ACTUATIONS = 2; // steering angle, acceleration

const int n_vars        = _mpc_steps * NUMBER_OF_STATES + (_mpc_steps - 1) * NUMBER_OF_ACTUATIONS; // number of state + actuation variables
const int n_constraints = _mpc_steps * NUMBER_OF_STATES;                                           // number of constraints

// where the first element of each state variable is stored in the vector to be feeded the optimization algorithm
const int _x_start      = 0;
const int _y_start      = _x_start      + _mpc_steps;
const int _theta_start  = _y_start      + _mpc_steps;
const int _v_start      = _theta_start  + _mpc_steps;
const int _cte_start    = _v_start      + _mpc_steps;
const int _etheta_start = _cte_start    + _mpc_steps;
const int _angvel_start = _etheta_start + _mpc_steps;
const int _a_start      = _angvel_start + _mpc_steps - 1;

// weights for cost computations
const double _w_cte      = 100;
const double _w_etheta   = 100;
const double _w_vel      = 1;
const double _w_angvel   = 100;
const double _w_accel    = 50;
const double _w_angvel_d = 10; // weight cost for high difference between consecutive angular velocity actuations
const double _w_accel_d  = 10; // weight cost for high difference between consecutive acceleration actuations

const double _max_angvel  = 3.0;   // Maximal angvel radian (~30 deg)
const double _max_accel   = 1.0;   // Maximal throttle accel
const double _bound_value = 1.0e3; // Bound value for other variables

class MPC {

 public:

  double angular_velocity;
  double acceleration;
  
  Dvector vars;                   // where all the state and actuation variables will be stored
  Dvector vars_lowerbound;        // lower limit for each corresponding variable in vars
  Dvector vars_upperbound;        // upper limit for each corresponding variable in vars
  Dvector constraints_lowerbound; // value constraint for each corresponding constraint expression
  Dvector constraints_upperbound; // value constraint for each corresponding constraint expression

  std::vector<double> future_xs;
  std::vector<double> future_ys;

  MPC();
  virtual ~MPC();

  // this function solves the model given the current state and path curve coefficients.
  vector<double> solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */