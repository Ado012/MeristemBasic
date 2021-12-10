#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "../solvers/baseSolver.h"

///
/// @brief An eight-order Runge-Kutta solver
///
class Dopr853 : public BaseSolver {

 private:
  
  double eps_;
  double h1_;

  static const double c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c14,c15,c16,
    b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,
    er1,er6,er7,er8,er9,er10,er11,er12,
    a21,a31,a32,a41,a43,a51,a53,a54,a61,a64,a65,a71,a74,a75,a76,
    a81,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,a101,a104,a105,
    a106,a107,a108,a109,a111,a114,a115,a116,a117,a118,a119,a1110,
    a121,a124,a125,a126,a127,a128,a129,a1210,a1211,a141,a147,a148,
    a149,a1410,a1411,a1412,a1413,a151,a156,a157,a158,a1511,a1512,
    a1513,a1514,a161,a166,a167,a168,a169,a1613,a1614,a1615,
    d41,d46,d47,d48,d49,d410,d411,d412,d413,d414,d415,d416,d51,d56,
    d57,d58,d59,d510,d511,d512,d513,d514,d515,d516,d61,d66,d67,d68,
    d69,d610,d611,d612,d613,d614,d615,d616,d71,d76,d77,d78,d79,
    d710,d711,d712,d713,d714,d715,d716;
  
 public:
  ///
  /// @brief Main constructor
  ///
  Dopr853(Organism *O,std::ifstream &IN);

	///
	/// @brief Reads the parameters used by the RK5Adaptive algorithm
	///
	/// This function is responsible for reading parameters used by the eight
	/// order Runge-Kutta algorithm. The parameter file sent to the simulator
	/// binary looks like:
	///
	/// <pre> 
	/// Dopr853
	/// T_start T_end 
	/// printFlag printNum 
	/// h1 eps 
	/// </pre> 
	///
	/// where Dopr853 is the identity string used by BaseSolver::getSolver
	/// to identify that the Dopr853 algorithm should be used. T_start
	/// (T_end) is the start (end) time for the simulation, printFlag is an
	/// integer which sets the output format (read by BaseSolver::print()),
	/// printNum is the number of equally spread time points to be printed. h1
	/// is the maximal (and initial) step size for each Dopr853 step, and eps sets
	/// the error threshold (should be <<1.0).
	///
	/// Comments can be included in the parameter file by starting the line with
	/// an #. Caveat: No check on the validity of the read data is applied.
	///
	/// @see BaseSolver::getSolver()
	/// @see BaseSolver::print()
	///
	void readParameterFile(std::ifstream &IN);
	
	void simulate(void);

	///
	/// @brief dopr853 adaptive stepper.
	///
	void rkqs(double hTry, double &hDid, double &hNext,
	std::vector< std::vector<double> > &yScal,
	std::vector< std::vector<double> > &yTemp,
	std::vector< std::vector<double> > &yErr,
	std::vector< std::vector<double> > &yErr2,
	std::vector< std::vector<double> > &k2,
	std::vector< std::vector<double> > &k3,
	std::vector< std::vector<double> > &k4,
	std::vector< std::vector<double> > &k5,
	std::vector< std::vector<double> > &k6,
	std::vector< std::vector<double> > &k7,
	std::vector< std::vector<double> > &k8,
	std::vector< std::vector<double> > &k9,
	std::vector< std::vector<double> > &k10,
	std::vector< std::vector<double> > &yTempRkck);

	///
	/// @brief One dopr853 step
	///
	void rkck(double h,
	std::vector< std::vector<double> > &yOut,
	std::vector< std::vector<double> > &yErr,
	std::vector< std::vector<double> > &yErr2,
	std::vector< std::vector<double> > &k2,
	std::vector< std::vector<double> > &k3,
	std::vector< std::vector<double> > &k4,
	std::vector< std::vector<double> > &k5,
	std::vector< std::vector<double> > &k6,
	std::vector< std::vector<double> > &k7,
	std::vector< std::vector<double> > &k8,
	std::vector< std::vector<double> > &k9,
	std::vector< std::vector<double> > &k10,
	std::vector< std::vector<double> > &yTempRkck);
		
};

///
/// @brief A fifth order Runge-Kutta solver for ODEs
///
class RK5Adaptive : public BaseSolver {

private:

	double eps_;
	double h1_;

public:
	///
	/// @brief Main constructor
	///
	RK5Adaptive(Organism *O,std::ifstream &IN);
	
	///
	/// @brief Reads the parameters used by the RK5Adaptive algorithm
	///
	/// This function is responsible for reading parameters used by the fifth
	/// order Runge-Kutta algorithm. The parameter file sent to the simulator
	/// binary looks like:
	///
	/// <pre> 
	/// RK5Adaptive
	/// T_start T_end 
	/// printFlag printNum 
	/// h1 eps 
	/// </pre> 
	///
	/// where RK5Adaptive is the identity string used by BaseSolver::getSolver
	/// to identify that the RK5Adaptive algorithm should be used. T_start
	/// (T_end) is the start (end) time for the simulation, printFlag is an
	/// integer which sets the output format (read by BaseSolver::print()),
	/// printNum is the number of equally spread time points to be printed. h1
	/// is the maximal (and initial) step size for each RK5 step, and eps sets
	/// the error threshold (should be <<1.0).
	///
	/// Comments can be included in the parameter file by starting the line with
	/// an #. Caveat: No check on the validity of the read data is applied.
	///
	/// @see BaseSolver::getSolver()
	/// @see BaseSolver::print()
	///
	void readParameterFile(std::ifstream &IN);
	
	void simulate(void);
	
	///
	/// @brief Runs a simulation of an organism model where mechanics is assumed
	/// fast
	///
	/// Runs a simulation where the mechanics equations are treated specifically
	/// to always be in equilibrium. At each step the mechanical equations are
	/// run until convergence.
	///
	void simulateWithConvergingMechanics(double cEps, double cTdiff, std::vector<size_t> &mechEq, std::vector<size_t> &mechPos);
	
	///
	/// @brief Updates mechanical equations until convergance 
	///
	/// Runs a simulation on mechanical equations only until it converges. It
	/// takes a system state as input and updates only positional variables
	/// using all mechanical equations until this system converges. During the
	/// process also the neighborhood is allowed to be updated, but no
	/// compartmentChanges are allowed for.
	///
	void updateMechanicsUntilConvergence(double cEps,double cTdiff, std::vector< std::vector<double> >  &yScal, std::vector< std::vector<double> >  &yTemp,
	 std::vector< std::vector<double> >  &yErr, std::vector< std::vector<double> >  &ak2, std::vector< std::vector<double> >  &ak3,
	 std::vector< std::vector<double> >  &ak4, std::vector< std::vector<double> >  &ak5, std::vector< std::vector<double> >  &ak6,
	 std::vector< std::vector<double> >  &yTempRkck, std::vector<size_t> &mechEq, std::vector<size_t> &mechPos);
	
	///
	/// @brief Fifth order Runge-Kutta adaptive stepper.
	///
	void rkqs(double hTry, double &hDid, double &hNext,
	std::vector< std::vector<double> > &yScal,
	std::vector< std::vector<double> > &yTemp,
	std::vector< std::vector<double> > &yErr,
	std::vector< std::vector<double> > &ak2,
	std::vector< std::vector<double> > &ak3,
	std::vector< std::vector<double> > &ak4,
	std::vector< std::vector<double> > &ak5,
	std::vector< std::vector<double> > &ak6,
	std::vector< std::vector<double> > &yTempRkck);

        ///
        /// @brief Fifth order Runge-Kutta adaptive stepper for mechanical update.
        ///
        void rkqsMechanical(double &t, double hTry, double &hDid, double &hNext,
        std::vector< std::vector<double> > &yScal,
        std::vector< std::vector<double> > &yTemp,
        std::vector< std::vector<double> > &yErr,
        std::vector< std::vector<double> > &ak2,
        std::vector< std::vector<double> > &ak3,
        std::vector< std::vector<double> > &ak4,
        std::vector< std::vector<double> > &ak5,
        std::vector< std::vector<double> > &ak6,
        std::vector< std::vector<double> > &yTempRkck,
        std::vector<size_t> &mechEq,
        std::vector<size_t> &posVar);
		
	///
	/// @brief One fifth order Runge-Kutta step
	///
	void rkck(double h,
	std::vector< std::vector<double> > &yOut,
	std::vector< std::vector<double> > &yErr,
	std::vector< std::vector<double> > &ak2,
	std::vector< std::vector<double> > &ak3,
	std::vector< std::vector<double> > &ak4,
	std::vector< std::vector<double> > &ak5,
	std::vector< std::vector<double> > &ak6,
	std::vector< std::vector<double> > &yTempRkck);
	
	        ///
        /// @brief One fifth order Runge-Kutta step for positional variables
        ///
        /// A fifth order Runge-Kutta step. Equivalent to the normal rkck, but only
        /// updates according to the equations provided in mechEq and updates only
        /// variables provided in posVar. Caveat: The check that mechEq only update
        /// positional variables are left to the user!
        ///
        void rkckMechanical(double h,
        std::vector< std::vector<double> > &yOut,
        std::vector< std::vector<double> > &yErr,
        std::vector< std::vector<double> > &ak2,
        std::vector< std::vector<double> > &ak3,
        std::vector< std::vector<double> > &ak4,
        std::vector< std::vector<double> > &ak5,
        std::vector< std::vector<double> > &ak6,
        std::vector< std::vector<double> > &yTempRkck,
        std::vector<size_t> &mechEq,
        std::vector<size_t> &posVar );

	double maxDerivative();
};

///
/// @brief A fifth order Runge-Kutta solver for ODEs
///
class RK5AdaptiveEquilibrium : public BaseSolver {

private:

	double eps_;
	double h1_;
	double threshold_;
	double printInterval_;


public:
	///
	/// @brief Main constructor
	///
	RK5AdaptiveEquilibrium(Organism *O,std::ifstream &IN);
	
	///
	/// @brief Reads the parameters used by the RK5AdaptiveEquilibrium algorithm
	///
	/// This function is responsible for reading parameters used by the fifth
	/// order Runge-Kutta algorithm. The parameter file sent to the simulator
	/// binary looks like:
	///
	/// <pre> 
	/// RK5AdaptiveEquilibrium
	/// T_start T_max
	/// printFlag printInterval
	/// h1 eps maxDeriv
	/// </pre> 
	///
	/// where RK5AdaptiveEquilibrium is the identity string used by BaseSolver::getSolver
	/// to identify that the RK5AdaptiveEquilibrium algorithm should be used. T_start
	/// is the start time for the simulation, printFlag is an
	/// integer which sets the output format (read by BaseSolver::print()), h1
	/// is the maximal (and initial) step size for each RK5 step, and eps sets
	/// the error threshold (should be <<1.0).
	///
	/// The state is printed with a time printInterval between the prints. If printInterval
	/// is set to zero, only the final state is printed.
	///
	/// The simulation continues until either the time T_max is reached or
	/// the maximimum derivative value is lower than maxDeriv.
	///
	/// Comments can be included in the parameter file by starting the line with
	/// an #. Caveat: No check on the validity of the read data is applied.
	///
	/// @see BaseSolver::getSolver()
	/// @see BaseSolver::print()
	///
	void readParameterFile(std::ifstream &IN);
	
	void simulate(void);
	
	///
	/// @brief Runs a simulation of an organism model where mechanics is assumed
	/// fast
	///
	/// Runs a simulation where the mechanics equations are treated specifically
	/// to always be in equilibrium. At each step the mechanical equations are
	/// run until convergence.
	///
	void simulateWithConvergingMechanics(double cEps, double cTdiff, std::vector<size_t> &mechEq, std::vector<size_t> &mechPos);
	
	///
	/// @brief Updates mechanical equations until convergance 
	///
	/// Runs a simulation on mechanical equations only until it converges. It
	/// takes a system state as input and updates only positional variables
	/// using all mechanical equations until this system converges. During the
	/// process also the neighborhood is allowed to be updated, but no
	/// compartmentChanges are allowed for.
	///
	void updateMechanicsUntilConvergence(double cEps,double cTdiff, std::vector< std::vector<double> >  &yScal, std::vector< std::vector<double> >  &yTemp,
	 std::vector< std::vector<double> >  &yErr, std::vector< std::vector<double> >  &ak2, std::vector< std::vector<double> >  &ak3,
	 std::vector< std::vector<double> >  &ak4, std::vector< std::vector<double> >  &ak5, std::vector< std::vector<double> >  &ak6,
	 std::vector< std::vector<double> >  &yTempRkck, std::vector<size_t> &mechEq, std::vector<size_t> &mechPos);
	
	///
	/// @brief Fifth order Runge-Kutta adaptive stepper.
	///
	void rkqs(double hTry, double &hDid, double &hNext,
	std::vector< std::vector<double> > &yScal,
	std::vector< std::vector<double> > &yTemp,
	std::vector< std::vector<double> > &yErr,
	std::vector< std::vector<double> > &ak2,
	std::vector< std::vector<double> > &ak3,
	std::vector< std::vector<double> > &ak4,
	std::vector< std::vector<double> > &ak5,
	std::vector< std::vector<double> > &ak6,
	std::vector< std::vector<double> > &yTempRkck);

        ///
        /// @brief Fifth order Runge-Kutta adaptive stepper for mechanical update.
        ///
        void rkqsMechanical(double &t, double hTry, double &hDid, double &hNext,
        std::vector< std::vector<double> > &yScal,
        std::vector< std::vector<double> > &yTemp,
        std::vector< std::vector<double> > &yErr,
        std::vector< std::vector<double> > &ak2,
        std::vector< std::vector<double> > &ak3,
        std::vector< std::vector<double> > &ak4,
        std::vector< std::vector<double> > &ak5,
        std::vector< std::vector<double> > &ak6,
        std::vector< std::vector<double> > &yTempRkck,
        std::vector<size_t> &mechEq,
        std::vector<size_t> &posVar);
		
	///
	/// @brief One fifth order Runge-Kutta step
	///
	void rkck(double h,
	std::vector< std::vector<double> > &yOut,
	std::vector< std::vector<double> > &yErr,
	std::vector< std::vector<double> > &ak2,
	std::vector< std::vector<double> > &ak3,
	std::vector< std::vector<double> > &ak4,
	std::vector< std::vector<double> > &ak5,
	std::vector< std::vector<double> > &ak6,
	std::vector< std::vector<double> > &yTempRkck);
	
	        ///
        /// @brief One fifth order Runge-Kutta step for positional variables
        ///
        /// A fifth order Runge-Kutta step. Equivalent to the normal rkck, but only
        /// updates according to the equations provided in mechEq and updates only
        /// variables provided in posVar. Caveat: The check that mechEq only update
        /// positional variables are left to the user!
        ///
        void rkckMechanical(double h,
        std::vector< std::vector<double> > &yOut,
        std::vector< std::vector<double> > &yErr,
        std::vector< std::vector<double> > &ak2,
        std::vector< std::vector<double> > &ak3,
        std::vector< std::vector<double> > &ak4,
        std::vector< std::vector<double> > &ak5,
        std::vector< std::vector<double> > &ak6,
        std::vector< std::vector<double> > &yTempRkck,
        std::vector<size_t> &mechEq,
        std::vector<size_t> &posVar );

	double maxDerivative();
};


///
/// @brief A ODE solver using a fifth order Runge-Kutta method with
/// possibility to compare the result to a cost template and read some
/// variables from a file.
///
/// A fifth order Runge-Kutta algorithm that allows for cost calculations and
/// reading some of the varables from a file.
///
class RK5AdaptiveTemplate : public BaseSolver {

private:

	double eps_;
	double h1_;

public:

	RK5AdaptiveTemplate(Organism *O,std::ifstream &IN);
	~RK5AdaptiveTemplate();
	
	///
	/// @brief Reads the parameters used by the RK5AdaptiveTemplate algorithm
	///
	/// This function is responsible for reading parameters used by the fifth
	/// order Runge-Kutta algorithm which allows for cost calculations and
	/// reading some of the varables from a file. The parameter file sent to the
	/// simulator binary looks like:
	///
	/// <pre> 
	/// RK5AdaptiveTemplate
	/// costType
	/// T_start T_end 
	/// printFlag printNum 
	/// h1 eps 
	/// variableInputBoolean
	/// numberOptTargets
	/// indexOfOptTarget_1 indexOfOptTarget_2 ... indexOfOptTarget_n
	/// costTemplateFile 
	/// </pre> 
	///
	/// Explanation of the paramers:
	/// <ul>
	///
	/// <li> RK5AdaptiveTemplate is the identity string used by 
	/// BaseSolver::getSolver to identify that the RK5AdaptiveTemplate algorithm
	/// should be used.
	///
	/// <li> costType is a string identifying which cost function is to be used.
	/// The available cost functions can be found in BaseCostFunction and the 
	/// specifier string should be equal to one of its sub-classes, but with 
	/// the first letter in lower case.
	///
	/// <li> T_start, and T_end, are floating point numbers specifying the 
	/// start and end timepoint for the simulation.
	///
	/// <li> printFlag is an integer which specifies which format the output 
	/// should be in. This is read by BaseSolver::print().
	///
	/// <li> printNum is the number of outputs the simulation will give. The 
	/// printouts are equally spaced between T_start and T_end. It may be worth
	/// knowing that this forces the solver to evaluate the simulation at specific
	/// time points, which may affect the step size the solver may take. 
	///
	/// <li> h1 is the maximal (and initial) step size for each RK5 step.
	///
	/// <li> eps (error per step) sets the error limit. This value should be 
	/// <<1.0. If it is set too high the simulations may be inaccurate, and they
	/// may outright fail if for example a too large step reduces a concentration
	/// to less than zero. If it is too small it will take a long time to finish
	/// the simulation.
	///
	/// <li> variableInputBoolean is a boolean (1 or 0) which specifies if
	/// external variable input should be read (usually set to 0).
	///
	/// <li> numberOptTargets is an integer which specifies how many of the 
	/// variables should be used for cost calculation.
	///
	/// <li> indexOfOptTarget are the indices (integers) for the variables which 
	/// should be used as cost targets. The number of indices specified must be 
	/// equal to numberOptTargets.
	///
	/// <li> costTemplateFile is a string specifying which cost template file to 
	/// be used.
	///
	/// </ul>
	///
	/// Example of file:
	/// <pre>
	/// RK5AdaptiveTemplate
	/// meanSquare		# The cost function to use
	/// 0 500			# range in which to calculate
	/// 1 101			# Which output followed by number of printouts
	/// 1.0 1e-12		# initial/maximum step size for the rk5Adaptive solver followed by error allowance
	/// 0			# not reading variable input
	/// 2			# two variables to be optimized against. 
	/// 0 2			# the index of the variables to be optimized against
	/// costFile.cost		# Cost template file name
	/// </pre>
	///
	// where RK5AdaptiveTemplate is the identity string used by
	// BaseSolver::getSolver to identify that the RK5AdaptiveTemplate algorithm
	// should be used. costType is the functional form (class) used for
	// calculating the 'cost' (objective function value, energy...) of a
	// simulation and possible functions for the moment are meanSquare and
	// meanSquareRelative. T_start (T_end) is the start (end) time for the
	// simulation, printFlag is an integer which sets the output format (read
	// by BaseSolver::print()), printNum is the number of equally spread time
	// points to be printed. h1 is the maximal (and initial) step size for each
	// RK5 step, and eps sets the error limit (should be <<1.0).
	///
	/// Comments can be included in the parameter file by prepending with
	/// an #. Caveat: No check on the validity of the read data is applied.
	///
	/// For optimizations it is advisable to set the printNum value to 0,
	/// as this printout otherwise clutters up the output.
	///
	/// @see BaseSolver::getSolver()
	/// @see BaseSolver::print()
	///	@see BaseCostFunction::createCostFunction()
	///
	void readParameterFile(std::ifstream &IN, int verbose=0);
	
	void simulate(void);
	
	///
	/// @brief Fifth order Runge-Kutta adaptiveTemplate stepper.
	///
	void rkqsTemplate(double hTry, double &hDid, double &hNext,
	std::vector< std::vector<double> > &yScal,
	std::vector< std::vector<double> > &yTemp,
	std::vector< std::vector<double> > &yErr,
	std::vector< std::vector<double> > &ak2,
	std::vector< std::vector<double> > &ak3,
	std::vector< std::vector<double> > &ak4,
	std::vector< std::vector<double> > &ak5,
	std::vector< std::vector<double> > &ak6,
	std::vector< std::vector<double> > &yTempRkck);
	
	///
	/// @brief One fifth order Runge-Kutta step
	///
	/// A fifth order Runge-Kutta step. It uses y,dydt as initial conditions at
	/// time t and takes a step of size h, returning the new function values in
	/// yOut. yErr is the estimated error in each variable (calculated as the
	/// difference of a fourth and fifth order step). The function and all of
	/// the parameters are taken from Numerical Recipes.
	///
	void rkckTemplate(double h,
	std::vector< std::vector<double> > &yOut,
	std::vector< std::vector<double> > &yErr,
	std::vector< std::vector<double> > &ak2,
	std::vector< std::vector<double> > &ak3,
	std::vector< std::vector<double> > &ak4,
	std::vector< std::vector<double> > &ak5,
	std::vector< std::vector<double> > &ak6,
	std::vector< std::vector<double> > &yTempRkck);
};



///
/// @brief This class solves the ODE using a fourth order Runge-Kutta solver
///
class RK4 : public BaseSolver {

private:

	double h_;
	unsigned int numStep_;

public:
	
	///
	/// @brief Main constructor.
	///
	RK4(Organism *O,std::ifstream &IN);
	
	///
	/// @brief Reads the parameters used by the RK4 algorithm
	///
	/// This function is responsible for reading parameters used by the fourth
	/// order Runge-Kutta algorithm. The parameter file sent to the simulator
	/// binary looks like:
	///
	/// <pre> 
	/// RK4 
	/// T_start T_end 
	/// printFlag printNum 
	/// h 
	/// </pre> 
	///
	/// where RK4 is the identity string used by BaseSolver::getSolver to
	/// identify that the RK4 algorithm should be used. T_start (T_end) is the
	/// start (end) time for the simulation, printFlag is an integer which sets
	/// the output format (read by BaseSolver::print()), printNum is the number
	/// of equally spread time points to be printed, and h is the step size for
	/// each RK4 step.
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
	/// @brief Runs a simulation of an organism model where the mechanics is
	/// updated until convergance between each step
	///
	void simulateWithConvergingMechanics(double cEps,double cTmax,
	 std::vector<size_t> &mechEq,
	 std::vector<size_t> &mechPos);
	
	///
	/// @brief The mechanics update function for the
	/// simulateWithConvergingMechanics solver
	///
	void updateMechanicsUntilConvergence(double cEps,double cTmax,
	 std::vector< std::vector<double> > 
	 &yt,
	 std::vector< std::vector<double> > 
	 &dyt,
	 std::vector< std::vector<double> > 
	 &dym,
	 std::vector<size_t> &mechEq,
	 std::vector<size_t> &mechPos );
	///
	/// @brief Fourth order Runge-Kutta stepper
	///
	void rk4(std::vector< std::vector<double> > &yt,
	 std::vector< std::vector<double> > &dyt,
	 std::vector< std::vector<double> > &dym );
	
	///
	/// @brief Fourth order Runge-Kutta stepper for the positional variables
	///
	void rk4Mechanical(std::vector< std::vector<double> > &yt,
	 std::vector< std::vector<double> > &dyt,
	 std::vector< std::vector<double> > &dym,
	 std::vector<size_t> &mechEq,
	 std::vector<size_t> &mechPos);
	double maxDerivative();
};

#endif /* RUNGEKUTTA_H */

