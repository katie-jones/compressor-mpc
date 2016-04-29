#ifndef TANK_H
#define TANK_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include "global_constants.h"

/**
 * Class containing tank used in compressor systems.
 */
class Tank {
 public:
  constexpr static int n_states = 1; ///< Number of tank states.
  constexpr static int n_inputs = 1; ///< Number of tank inputs.

  /// Type describing input to tank.
  typedef Eigen::Array<double, n_inputs, 1> TankInput;

  /// Type describing state of tank.
  typedef Eigen::Array<double, n_states, 1> TankState;

  /// Parameters determining dynamics of tank.
  struct Params {
    double pout;
    double volume;
    Vec<8> D;
    double m_out_c;
    Params();
    Params(const Params &x);
  }; 

  /// Specify parameters to use.
  Tank(Params params) : params(params) {}

  /// Use default parameters.
  Tank() : params(Params()) {}

  /**
   * Get derivative of tank.
   * Calculate derivative at state x and input u. Total mass flow into tank is given by mass_flow_compressors.
   */
  TankState GetDerivative(const TankState x, const TankInput u,
                          const double mass_flow_compressors) const;

  /// Return default input to tank.
  static inline TankInput GetDefaultInput() {
    return ((TankInput() << 0.7).finished());
  }
  
  /// Return default state of tank.
  static inline TankState GetDefaultState() {
    return ((TankState() << 1.12).finished());
  }

 private:
  const Params params;
};

#endif
