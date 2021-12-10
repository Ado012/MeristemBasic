//
// Filename     : baseSolver.h
// Description  : Base class for solvers 
// Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
// Created      : August 2006
// Revision     : $Id: baseSolver.h 474 2010-08-30 15:23:20Z henrik $
//
#ifndef SOLVER_H
#define SOLVER_H

#include "../organism.h"
#include "../cost/costFunction.h"

///
/// @brief A factory class for classes describing different numerical solvers
/// for the ordinary differential equations
///
/// Detailed description to come...
///
class BaseSolver {

protected:
	Organism *O_;
	BaseCostFunction *C_;
	std::vector< std::vector<double> > y_;
	std::vector< std::vector<double> > dydt_; 
	std::vector< std::vector<double> > previousValue_;
	std::vector< std::vector<double> > futureValue_;
	std::vector< std::vector<std::vector<double> > > yCopy_;
public: double t_; // To be accessible to save the values of the simulator inside the estiöator saved_init container
protected :
	double startTime_;
	double endTime_;
	int printFlag_;
	int numPrint_;
	unsigned int numOk_, numBad_;
	bool debugFlag_;
	std::vector<int> simulationFlag_, simulationIndex_;
	std::vector<double> previousTime_;
	double futureTime_;
	std::string inputTemplateFile_;
	size_t numSimulation_;
	
public:
	BaseSolver();
	BaseSolver(Organism *O,std::ifstream &IN);
	virtual ~BaseSolver();
	
	///
	/// @brief This function implements the factory method for initiating a
	/// numerical solver.
	///
	/// This function reads a solver id in the provided file and initiates a
	/// numerical solver of the id type, which is returned. The Organism (model)
	/// pointer is used to introduce the model to the solver. The initiation is
	/// continued by the individual solver classes since they use different
	/// parameters (see links below). The parameter file sent to the simulator
	/// binary looks like:
	///
	/// @verbatim 
	/// solverId 
	/// parameter_i ...
	/// ...
	/// @endverbatim
	///
	/// where the solverId is the name of the numerical method (class) used. The
	/// parameters used can be found in the links below which lists the
	/// currently available methods/classes.
	///
	/// @see RK5Adaptive::readParameterFile()
	/// @see RK4::readParameterFile()
	/// @see RK5AdaptiveTemplate::readParameterFile()
	///
	static BaseSolver *getSolver(Organism *O, const std::string &file);
	static BaseSolver *getSolver(Organism *O, std::ifstream &IN);
	
	size_t debugCount() const;
	void getInit();
	
	///
	/// @brief General printing function
	///
	/// Caveat: Not yet general, but will be...
	///
	void print(std::ostream &os=std::cout);
	void setSimulationFlag(int verbose=0);
	void printInit(std::ostream &os) const;
	void printDebug(std::ostream &os) const;
	double maxDerivative();
	virtual void readParameterFile(std::ifstream &IN);
	virtual void simulate(void);
	
	inline size_t N() const;
	inline size_t M() const;
	inline double startTime() const;
	inline double endTime() const;
	bool debugFlag() const;
	inline std::vector< std::vector<double> > & y();
	inline std::vector< std::vector<double> > & dydt();
	inline std::vector<double> & y(size_t i);
	inline std::vector<double> & dydt(size_t i);
	inline double y(size_t i,size_t j) const;
	inline double dydt(size_t i,size_t j) const;
	inline int simulationFlag( int i ) const;
	inline const std::vector<int> &simulationFlag() const;
	inline std::string inputTemplateFile() const;
	inline double previousTime( int i ) const;
	inline double futureTime() const;
	inline double previousValue( int i, int j ) const;
	inline double futureValue( int i, int j ) const;
	inline double cost() const;
	inline void printCost(std::ostream &os=std::cerr) const;
	inline void readInit(const std::string &initFile);
	inline Organism *getOrganism();
	inline void setOrganism(Organism *O);
	inline void setY(size_t i,size_t j,double value);
	inline void setDydt(size_t i,size_t j,double value);
	
	inline void setCostTemplateTime( const std::vector<double> &inputVector );
	inline void setCostTemplate( const std::vector< 
															 std::vector< std::vector<double> > > 
															 &costTemplate );
	inline void setCostTemplateFile( const std::string &value );
	inline void setCostList(const std::vector<size_t> &inputVector);
	
	///
	/// @brief Checks that numbers in y_ and dydt_ are not NaN or Inf
	///	
	int numInvalidValues();

	//Functions to be used when template data is read from file
	///
	/// @brief Divides an compartment in the organism and adjusts template
	/// values
	///
	void divideTemplate(int compartment);
	
	///
	/// @brief Initial reading of template data
	///
	/// The template data is read from a file. The number of cells should be the
	/// same as when defined from the organism. The data for the correct time
	/// point must be found.
	///
	void initiateValuesTemplate( std::ifstream &IN,
															 std::vector<int> &divIndex,
															 std::vector<double> &divTime);
	
	///
	/// @brief Update border values for variables which are read from file
	///
	/// This function updates the values of variables that are not simulated but
	/// read from file. It uses variable values for one previous time point and
	/// one future time point, such that variable values can be interpolated
	/// from this information. For the update it copies old future values into
	/// the previous values variables, and if the file is still readable, it
	/// reads another data point. Then data for variables in between the two
	/// time points are interpolated.
	///
	void updateValuesTemplate( std::ifstream &IN,double tmax,
														 std::vector<int> &divIndex,
														 std::vector<double> &divTime);
	///
	/// @brief Update all values in the y matrix using previous/futureValues
	///
	void allValuesTemplate();
	
	///
	/// @brief Returns a value interpolated from values read in a template file
	///
	double valueTemplate( int i, int j);	
};

// Number of cells (rows in the variable matrix)
inline size_t BaseSolver::N() const
{
	return y_.size();
}

// Number of species plus four (cols in the variable matrix)
inline size_t BaseSolver::M() const
{ 
	return (N() ? y_[0].size() : 0);
}

// Variable values stored cell by cell
inline std::vector< std::vector<double> > & BaseSolver::y()
{ 
	return y_;
}

// Derivative values stored cell by cell
inline std::vector< std::vector<double> > & BaseSolver::dydt()
{ 
	return dydt_;
}

inline std::vector<double> & BaseSolver::y(size_t i)
{
	return y_[i];
}

inline std::vector<double> & BaseSolver::dydt(size_t i)
{
	return dydt_[i];
}

inline double BaseSolver::y(size_t i,size_t j) const
{
	return y_[i][j];
}

inline double BaseSolver::dydt(size_t i,size_t j) const
{
	return dydt_[i][j];
}

// Flag marking if a variable is updated by simulation(1) or from file(0)
inline int BaseSolver::simulationFlag( int i ) const 
{
	return simulationFlag_[i];
}

// Returns the complete simulationFlag vector
inline const std::vector<int> &BaseSolver::simulationFlag() const 
{
	return simulationFlag_;
}

// File from which variable values are read (for update)
inline std::string BaseSolver::inputTemplateFile() const 
{
	return inputTemplateFile_;
}

inline double BaseSolver::startTime() const
{
	return startTime_;
}

inline double BaseSolver::endTime() const 
{
	return endTime_;
}

inline bool BaseSolver::debugFlag() const
{
	return debugFlag_;
}

// Closest previous time point for a specific cell when updating with template
inline double BaseSolver::previousTime( int i ) const 
{
	return previousTime_[i];
}

// Closest future time point when updating variables from template file
inline double BaseSolver::futureTime() const 
{ 
	return futureTime_;
}

// Closest previous variable value when updating variable from template file
inline double BaseSolver::previousValue( int i, int j ) const 
{
	return previousValue_[i][j];
}

// Closest future variable value when updating variable from template file
inline double BaseSolver::futureValue( int i, int j ) const
{ 
	return futureValue_[i][j];
}

inline double BaseSolver::cost() const 
{ 
	if(C_!=0)
		return C_->cost();
	//std::cerr << "warning Solver::cost() returns -1\n";
	return -1;
}

inline void BaseSolver::printCost(std::ostream &os) const 
{ 
	if(C_!=0)
		C_->printCost(os);
}

inline void BaseSolver::readInit(const std::string &initFile)
{
	O_->readInit(initFile);
}

inline Organism *BaseSolver::getOrganism()
{
	return O_;
}

inline void BaseSolver::setOrganism(Organism *O)
{
	O_= O;
}

inline void BaseSolver::setY(size_t i,size_t j,double value)
{
	y_[i][j]=value;
}

inline void BaseSolver::setDydt(size_t i,size_t j,double value)
{
	dydt_[i][j]=value;
}

inline void BaseSolver::setCostTemplateFile( const std::string &value ) 
{
	C_->setCostTemplateFile(value); 
}

inline void BaseSolver::
setCostTemplate( const std::vector< std::vector< std::vector<double> > >  &costTemplate )
{
	C_->setCostTemplate(costTemplate); 
}

inline void BaseSolver::setCostTemplateTime(const std::vector<double> &inputVector ) 
{
	C_->setCostTemplateTime(inputVector); 
}

inline void BaseSolver::setCostList(const std::vector<size_t> &inputVector) 
{
	C_->setCostList(inputVector);
}

#endif /* SOLVER_H */

