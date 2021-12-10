#ifndef BASECOMPARTMENTNEIGHBORHOOD_H
#define BASECOMPARTMENTNEIGHBORHOOD_H
//
// Filename     : baseCompartmentNeighborhood.h
// Description  : The common base for classes describing comp updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2005
// Revision     : $Id: baseCompartmentNeighborhood.h 643 2015-08-20 08:08:11Z henrik $
//

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../compartment.h"

class Organism;

///
/// @brief A base class for describing a compartment neighborhood
///
/// The BaseCompartmentNeighborhood class is a base class used to define
/// different types of compartmental neighborhoods. Each class uses a vector
/// of parameters and variable indeces to create and update which compartments
/// are neighbors to each other. The variable indices are devided into
/// multiple layers to allow for different types of contributions for
/// different variables. The baseCompartmentNeighborhood class can be seen as
/// a 'factory' creating different types of neighborhood definitions.
///
class BaseCompartmentNeighborhood {
  
 private:
  
  std::string id_;
  std::vector<double> parameter_;           
  std::vector<std::string> parameterId_;           
  std::vector< std::vector<size_t> > variableIndex_;
  double previousTime_;
  unsigned int numNeighbor_;

 public:
  
	///
	/// @brief Main factory creator
	///
	/// All creation of compartmentNeighborhood classes should be mapped onto
	/// this function. It calls the correct class constructor defined by the id
	/// string.
	///
	/// @param paraValue Vector of parameters used by the neighborhood
	/// definition.
	///
	/// @param indValue  Vector of vector of variable indices used for
	/// neighborhood definition.
	///
	/// @param idValue  Unique string for identifying which
	/// compartmentNeighborhood class that should be used.
	///
	static BaseCompartmentNeighborhood* 
	createCompartmentNeighborhood(std::vector<double> &paraValue, 
																std::vector< std::vector<size_t> > 
																&indValue, const std::string &idValue );

	///
	/// @brief Factory creator function which reads correct data for a
	/// compartmentNeighborhood and calls the main creator function.
	///
	/// It takes an istream and defines appropriate variables to be able to call
	/// the main creator.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
	static BaseCompartmentNeighborhood* 
	  createCompartmentNeighborhood( std::istream &IN ); 
	
	//No constructor defined!
	//BaseCompartmentNeighborhood();

	///
	/// @brief Virtual destructor (empty implementation)
	///
	virtual ~BaseCompartmentNeighborhood();
	
	BaseCompartmentNeighborhood & 
	  operator=( const BaseCompartmentNeighborhood & 
		     baseCompartmentNeighborhoodCopy );
	
	///
	/// @brief Returns the id string
	///
	inline std::string id() const;
	
	///
	/// @brief Returns the number of parameters
	///
	inline size_t numParameter() const;  
	
	///
	/// @brief Returns the number of variable-index levels
	///
	inline size_t numVariableIndexLevel() const;
	
	///
	/// @brief Returns the number of variable-index within one level 
	///
	inline size_t numVariableIndex(size_t level) const;
	
	///
	/// @brief Returns a parameter
	///
	inline double parameter(size_t i) const;
	
	///
	/// @brief Returns a reference of a parameter
	///
	inline double& parameterAddress(size_t i);
	
	///
	/// @brief Returns the id for a parameter
	///
	inline std::string parameterId(size_t i) const;
	
	///
	/// @brief Returns the index of variable used in the update equation
	///
	inline size_t variableIndex(size_t i,size_t j) const;
	
	///
	/// @brief Returns the previous time for update
	///
	inline double previousTime() const;
	
	///	
	/// @brief Returns the current number of neighbors
	///
	inline unsigned int numNeighbor() const;
	
	///
	/// @brief Sets the id string
	///
	inline void setId(const std::string &value);
	
	///
	/// @brief Sets the parameter value with index i
	///
	inline void setParameter(size_t i,double value);
	
	///
	/// @brief Sets all parameters within a compartmentNeighborhood
	///
	inline void setParameter(std::vector<double> &value);
	
	///
	/// @brief Sets the parameter id value with index i
	///
	inline void setParameterId(size_t i,const std::string &value);
	
	///
	/// @brief Sets all parameter id's within a compartmentNeighborhood
	///
	inline void setParameterId(std::vector<std::string> &value);
	
	///
	/// @brief Sets the variable index at specific level and place
	///
	inline void setVariableIndex(size_t i, size_t j,size_t value);
	
	///
	/// @brief Sets the variable indeces for a given level
	///
	inline void setVariableIndex(size_t i, std::vector<size_t> &value);
	
	///
	/// @brief Sets the variable indeces used in a compartmentNeighborhood
	///
	inline void setVariableIndex(std::vector< std::vector<size_t> > &value);
	
	///
	/// @brief Sets the time when updated
	///
	inline void setPreviousTime(double value);  
	
	///
	/// @brief Sets the current number of neighbors
	///
	inline void setNumNeighbor(unsigned int value);  
	
	///
	/// @brief Initial creation of neighborhood before simulation is started.
	///
	/// In this function each class defines its rule for initiating a
	/// neighborhood. It is called before the start of a simulation, defines the
	/// neighbors, and saves the neighborhood relations in the organism class
	/// variable O.neighbor() using O.setNeighbor().
	///
	/// @see Organism::setNeighbor()
	///
	/// @param O Model defined in Organism.
	///
	/// @param y Current state of the system as defined in the solver.
	///
	/// @param t Current time used for the ability to create a @f$\Delta t@f$
	/// between neighborhood updates.
	///
	virtual unsigned int create(Organism &O,std::vector< std::vector<double> > &y,
				    double t=0.0 );
	
	///
	/// @brief Updates neighborhood during simulation.
	///
	/// In this function each class defines its rule for updating a neighborhood
	/// during a simulation. It is called between each time step update of the
	/// system during a simulation, updates the neighborhood, and saves the
	/// neighborhood relations in the organism class variable
	/// Organism::neighbor() using Organism::setNeighbor(). Often it is
	/// implemented by calling the baseCompartmentNeighborhood::create()
	/// function.
	///
	/// @see Organism::setNeighbor()
	///
	/// @param O Model defined in Organism.
	///
	/// @param y Current state of the system as defined in the solver.
	///
	/// @param t Current time used for the ability to create a @f$\Delta t@f$
	/// between neighborhood updates.
	///	
	virtual unsigned int update(Organism &O, 
				    std::vector< std::vector<double> > &y,
				    double t=0.0 );  
	///
	/// @brief Prints the Neighborhood reaction in Cambium format
	///
	/// 
	///
	/// @see Organism::printCambium()
	///
	virtual void printCambium( std::ostream &os ) const;	
};

inline std::string BaseCompartmentNeighborhood::id() const 
{
  return id_;
}

inline size_t BaseCompartmentNeighborhood::numParameter() const 
{
  return parameter_.size();
}

inline size_t BaseCompartmentNeighborhood::
numVariableIndexLevel() const 
{
  return variableIndex_.size();
}

inline size_t BaseCompartmentNeighborhood::
numVariableIndex(size_t level) const 
{
  return variableIndex_[level].size();
}

inline double BaseCompartmentNeighborhood::parameter(size_t i) const 
{
  return parameter_[i];
}

inline double& BaseCompartmentNeighborhood::parameterAddress(size_t i) 
{
  return parameter_[i];
}

inline std::string BaseCompartmentNeighborhood::
parameterId(size_t i) const 
{
  return parameterId_[i];
}

inline size_t BaseCompartmentNeighborhood::
variableIndex(size_t i,size_t j) const 
{
  return variableIndex_[i][j];
}

inline double BaseCompartmentNeighborhood::previousTime() const 
{
  return previousTime_;
}

inline unsigned int BaseCompartmentNeighborhood::
numNeighbor() const 
{
  return numNeighbor_;
}

inline void BaseCompartmentNeighborhood::
setId(const std::string &value) 
{
  id_=value;
}

inline void BaseCompartmentNeighborhood::
setParameter(size_t i,double value) 
{
  parameter_[i]=value;
}

inline void BaseCompartmentNeighborhood::
setParameter(std::vector<double> &value) 
{
  parameter_ = value;
}

inline void BaseCompartmentNeighborhood::
setParameterId(size_t i,const std::string &value) 
{
  parameterId_[i]=value;
}

inline void BaseCompartmentNeighborhood::
setParameterId(std::vector<std::string> &value) 
{
  parameterId_=value;
}

inline void BaseCompartmentNeighborhood::
setVariableIndex(size_t i, size_t j,size_t value) 
{
  variableIndex_[i][j] = value;
}

inline void BaseCompartmentNeighborhood::
setVariableIndex(size_t i, std::vector<size_t> &value) 
{
  variableIndex_[i] = value;
}

inline void BaseCompartmentNeighborhood::
setVariableIndex(std::vector< std::vector<size_t> > &value) 
{
  variableIndex_ = value;
}

inline void BaseCompartmentNeighborhood::
setPreviousTime(double value) 
{
  previousTime_=value;
}

inline void BaseCompartmentNeighborhood::
setNumNeighbor(unsigned int value) 
{
  numNeighbor_=value;
}

#endif


