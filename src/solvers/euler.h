#ifndef EULER_H
#define EULER_H
//
// Filename     : euler.h
// Description  : Simple Euler numerical solvers
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : March 2009
// Revision     : $Id:$
//

#include "baseSolver.h"

///
/// @brief A simple Euler step solver class
///
class Euler : public BaseSolver {

 private:
  
  double h_;
  
 public:
  ///
  /// @brief Main constructor
  ///
  Euler(Organism *O,std::ifstream &IN);
  
  ///
  /// @brief Reads the parameters used by the Euler solver
  ///
  /// This function is responsible for reading parameters used by the
  /// Euler algorithm. The parameter file sent to the simulator
  /// binary looks like:
  ///
  /// <pre> 
  /// Euler
  /// T_start T_end 
  /// printFlag printNum 
  /// h 
  /// </pre> 
  ///
  /// where Euler is the identity string used by BaseSolver::getSolver
  /// to identify that the Euler algorithm should be used. T_start
  /// (T_end) is the start (end) time for the simulation, printFlag is an
  /// integer which sets the output format (read by BaseSolver::print()),
  /// printNum is the number of equally spread time points to be printed. h
  /// is the step size for each Euler step.
  ///
  /// Comments can be included in the parameter file by starting the line with
  /// an #. Caveat: No check on the validity of the read data is applied.
  ///
  /// @see BaseSolver::getSolver()
  /// @see BaseSolver::print()
  ///
  void readParameterFile(std::ifstream &IN);
	
  ///
  /// @brief Runs a simulation of an organism model
  ///
  void simulate(void);
	
  ///
  /// @brief An Euler step
  ///
  void euler();
};

///
/// @brief A Euler step solver class simulating until an equilibrium
/// condition has been met.
///
class EulerEquilibrium : public BaseSolver {

 private:
  
  double h_;
  double threshold_;
  double printInterval_;
  
 public:
  ///
  /// @brief Main constructor
  ///
  EulerEquilibrium(Organism *O,std::ifstream &IN);
  
  ///
  /// @brief Reads the parameters used by the Euler solver
  ///
  /// This function is responsible for reading parameters used by the
  /// Euler algorithm. The parameter file sent to the simulator
  /// binary looks like:
  ///
  /// <pre> 
  /// EulerEquilibrium
  /// T_start Tmax
  /// printFlag printInterval 
  /// h tolDeriv
  /// </pre> 
  ///
  /// where Euler is the identity string used by BaseSolver::getSolver
  /// to identify that the Euler algorithm should be used. T_start
  /// is the start time for the simulation and Tmax is the maximum
  /// allowed simulation time, printFlag is an integer which sets the
  /// output format (read by BaseSolver::print()),
  /// printInterval is the time between printing the outputs. If
  /// printInterval is set to zero, only the final state is printed.
  ///
  /// h is the step size for each Euler step and tolDeriv is the
  /// allowed value of the maximum derivative.
  /// The simulation is stopped when the maximum derivative is below
  /// tolDeriv, or if the simulation time Tmax is reached.
  ///
  ///
  /// Comments can be included in the parameter file by starting the line with
  /// an #. Caveat: No check on the validity of the read data is applied.
  ///
  /// @see BaseSolver::getSolver()
  /// @see BaseSolver::print()
  ///
  void readParameterFile(std::ifstream &IN);
	
  ///
  /// @brief Runs a simulation of an organism model
  ///
  void simulate(void);
	
  ///
  /// @brief An Euler step
  ///
  void eulerEquilibrium();
};


///
/// @brief A simple Euler step solver class, where negative values are truncated to zero and NaN Inf 
/// truncated to maxVal_ (parameter given)
///
/// WARNING! This solver will give numerical errors, and the results
/// cannot be trusted. It is implemented for use as a quick test of a
/// simulation.
///
/// The algorithm takes Euler steps and truncates values at zero and NaN,Inf
/// i.e. negative values will be set to zero and NaN,Inf are set to user provided maxVal.
///
class EulerTruncated : public BaseSolver {
  
 private:
  
  double h_;
  double maxVal_;
  
 public:
  ///
  /// @brief Main constructor
  ///
  EulerTruncated(Organism *O,std::ifstream &IN);
  
  ///
  /// @brief Reads the parameters used by the Euler solver
  ///
  /// This function is responsible for reading parameters used by the
  /// Euler algorithm. The parameter file sent to the simulator
  /// binary looks like:
  ///
  /// <pre> 
  /// EulerTruncated
  /// T_start T_end 
  /// printFlag printNum 
  /// h maxVal
  /// </pre> 
  ///
  /// where EulerTruncated is the identity string used by BaseSolver::getSolver
  /// to identify that the Euler algorithm should be used. T_start
  /// (T_end) is the start (end) time for the simulation, printFlag is an
  /// integer which sets the output format (read by BaseSolver::print()),
  /// printNum is the number of equally spread time points to be printed. h
  /// is the step size for each Euler step. maxVal is the (LARGE) value replacing
  /// NaN and Inf generated in the simulation.
  ///
  /// Comments can be included in the parameter file by starting the line with
  /// an #. Caveat: No check on the validity of the read data is applied.
  ///
  /// @see BaseSolver::getSolver()
  /// @see BaseSolver::print()
  ///
  void readParameterFile(std::ifstream &IN);
	
  ///
  /// @brief Runs a simulation of an organism model
  ///
  void simulate(void);
	
  ///
  /// @brief An Euler step with truncation of negative numbers to zero
  ///
  void eulerTrunc();
};

#endif /* EULER_H */

