//
// Filename     : baseReaction.h
// Description  : The common base for classes describing reaction updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2003
// Revision     : $Id: baseReaction.h 646 2015-08-23 18:35:06Z henrik $
//
#ifndef BASEREACTION_H
#define BASEREACTION_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "../compartment.h"
#include "../common/typedefs.h"

class Organism;

///
/// @brief A factory class for classes describing differential equation
/// updates for dynamical variables
///
/// The BaseReaction class is a base (factory) class used when defining
/// different types of "reaction" classes. Each reaction class uses a vector
/// of parameters and variable indices to calculate a derivative of model
/// variables. The variable indices are divided into multiple layers to allow
/// for different types of contributions for different variables. The
/// BaseReaction class can be seen as a 'factory' creating reactions of
/// different types.
///
class BaseReaction {
	
 public:
  
  //BaseReaction();//Caveat: Constructor not defined!
  
  ///
  /// @brief Empty destructor.
  ///
  virtual ~BaseReaction();
  
  ///
  /// @brief Unimplemented operator=
  ///
  BaseReaction & operator=( const BaseReaction & baseReactionCopy );
  
  ///
  /// @brief Main factory creator, all creation should be mapped onto this one
  ///
  /// Given the idValue a reaction of the defined type is returned (using new
  /// Class). It chooses from the user defined list of possible reactions, and
  /// returns the correct one if defined. If the idValue given is not a defined
  /// reaction it exits.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///
  /// @param idValue identification of which reaction that should be created
  ///
  /// @return Returns a pointer to an instance of a reaction class defined by
  /// the idValue string.
  ///
  static BaseReaction* createReaction(std::vector<double> &paraValue, 
				      std::vector< std::vector<size_t> > 
				      &indValue,
				      const std::string &idValue );
  ///
  /// @brief This creator reads from an open file and then calls for the main
  /// factory creator
  ///
  /// @see createReaction(std::vector<double>&,std::vector<
  /// std::vector<size_t> >&,string)
  ///
  static BaseReaction* createReaction( std::istream &IN ); 
  
  ///
  /// @brief Returns the Organism pointer
  ///
  Organism* organism() const;
  
  ///
  /// @brief Returns the id string
  ///
  inline std::string id() const;
  
  ///
  /// @brief Returns the number of parameters
  ///
  inline size_t numParameter() const;  
  
  ///
  /// @brief Returns the number of variable index levels
  ///
  inline size_t numVariableIndexLevel() const;
  
  ///
  /// @brief Returns the number of variable-index within one level 
  ///
  inline size_t numVariableIndex(size_t level) const;
  
  ///
  /// @brief Returns the parameter at index i
  ///
  inline double parameter(size_t i) const;
  
  ///
  /// @brief Returns a reference to parameter at index 
  ///
  inline double& parameterAddress(size_t i);
  
  ///
  /// @brief Returns the id for a parameter
  ///
  inline std::string parameterId(size_t i) const;
  
  ///
  /// @brief Returns the index of variable used in the derivs equation
  ///
  inline size_t variableIndex(size_t i,size_t j) const;
  
  ///
  /// @brief Sets the Organism pointer
  ///
  inline void setOrganism(Organism* value);
  
  ///
  /// @brief Sets the id string
  ///
  inline void setId(const std::string &value);
  
  ///
  /// @brief Sets the parameter with index i to value
  ///
  inline void setParameter(size_t i,double value);
  
  ///
  /// @brief Sets all parameters within a reaction
  ///
  inline void setParameter(std::vector<double> &value);
  
  ///
  /// @brief Sets the parameter id with index i to value
  ///
  inline void setParameterId(size_t i,const std::string &value);
  
  ///
  /// @brief Sets all parameter id's within a reaction
  ///
  inline void setParameterId(std::vector<std::string> &value);
  
  ///
  /// @brief Sets the variable index at specific level and place to value
  ///
  inline void setVariableIndex(size_t i, size_t j,size_t value);
  
  ///
  /// @brief Sets all variable indices for level i
  ///
  inline void setVariableIndex(size_t i, std::vector<size_t> &value);
  
  ///
  /// @brief Sets the variable indeces used in a reaction
  ///
  inline void setVariableIndex(std::vector< std::vector<size_t> > &value);
  
  ///
  /// @brief Gives the total number of species
  ///
  /// can not be inline because it needs the total definition of the organism class, which depends then of this
  /// class,
  /// so, to avoid circular inclusions, we can not use the fields of organism beacause organism is not included
  /// here. We just know it as the name of a class
  /// size_t numSpecies();

  ///
  /// @brief For reactions that need to initiate its parameters (state) before
  /// the simulation
  ///
  /// Occasionally, a reaction have 'states' included as parameters, or affect
  /// model parameters. If these
  /// need to be initiated before the simulation this is done by this
  /// function. Update of the state variables can also be done. Also for
  /// example some approximate stochastic updates can be done in this
  /// function. Model parameters are reachable via organism()->parameter(i).
  ///
  /// @param t The current time (start time)
  /// @param y The current state which can be updated
  ///
  /// @note This is typically only included for very specific reactions.
  ///
  virtual void initiate(double t, DataMatrix &y);
  
  ///
  /// @brief Calculates the derivative given the state in y and adds it to dydt
  ///
  /// This function, together with the constructor are the main functions
  /// needed to create a new reaction class. It derives the time derivative
  /// contribution from a reaction given the state in y via the indices given
  /// in the member variableIndex and parameters given in the member parameter
  /// and adds the result to dydt. This is done for a single Compartment and
  /// this function is called from Organism::derivs(), Species::derivs(), and
  /// Topology::derivs().
  ///
  /// @see Organism::derivs()
  ///
  /// @see Species::derivs()
  ///
  /// @see Topology::derivs()
  ///
  /// @param compartment The compartment (e.g. cell) which should be updated.
  ///
  /// @param species The index of the variable (e.g. molecule) that should be
  /// updated.
  ///
  /// @param y The current state used for the update.
  ///
  /// @param dydt The derivs to which the output of the reaction is added.
  ///
  virtual void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);

  ///
  /// @brief For reactions that need to update its parameters (state) during
  /// the simulation
  ///
  /// Occasionally, a reaction have 'states' included as parameters. If these
  /// need to be updated during the simulation this is done by this
  /// function. Update of the state variables can also be done. Also for
  /// example some approximate stochastic updates can be done in this
  /// function.
  ///
  /// @param h The time step taken by the solver.
  /// @param t The current time
  /// @param y The current state which can be updated
  ///
  /// @note This is typically only included for very specific reactions.
  ///
  virtual void update(double h, double t, DataMatrix &y);

  ///
  /// @brief Calculates the derivative given the state in y and adds it to dydt and its abs value to sdydt
  ///
  /// This derivatives function is needed for those reactions that should be
  /// able to run updates with a noise added. It derives the time derivative 
  /// contribution from a reaction given the state in y via the indices given
  /// in the member variableIndex and parameters given in the member parameter
  /// and adds the result to dydt. At the same time is adds the absolute value of the derivative 
  /// contribution into sdydt for noise additions. 
  ///
  /// @see Organism::derivsWithAbs()
  ///
  /// @see Species::derivsWithAbs()
  ///
  /// @see Topology::derivsWithAbs()
  ///
  /// @param compartment The compartment (e.g. cell) which should be updated.
  ///
  /// @param species The index of the variable (e.g. molecule) that should be
  /// updated.
  ///
  /// @param y The current state used for the update.
  ///
  /// @param dydt The derivs to which the output of the reaction is added.
  ///
  /// @param sdydt The absolute value of the derivs to which the output of the reaction is added.
  ///
  virtual void derivsWithAbs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt);

  ///
  /// @brief Calculates the Jacobian given the state in y and adds it to A
  ///
  /// This function, together with the constructor and derivative
  /// function are the main functions
  /// needed to create a new reaction class if it to be updated using
  /// an implicit numerical solver. It derives the contribution to the Jacobian
  /// from a reaction given the state in y via the indices given
  /// in the member variableIndex and parameters given in the member parameter
  /// and adds the result to A. This is done for a single Compartment and
  /// this function is called from Organism::derivs(), Species::derivs(), and
  /// Topology::derivs().
  ///
  /// @see Organism::Jacobian()
  ///
  /// @see Species::Jacobian()
  ///
  /// @see Topology::Jacobian()
  ///
  /// @param compartment The compartment (e.g. cell) which should be updated.
  ///
  /// @param species The index of the variable (e.g. molecule) that should be
  /// updated.
  ///
  /// @param y The current state used for the update.
  ///
  /// @param A The Jacobian matrix to which the output of the reaction
  /// is added.
  ///
  ///
  /*
  virtual size_t Jacobian(Compartment &compartment,size_t species,DataMatrix &y,JacobianMatrix &A);
  
  ///
  /// @brief Calculates the propensity for a specific reaction.
  ///
  /// The state (and reaction information) is used to calculate the
  /// propensity (the probability for a reaction to happen whithin t
  /// and t+dt). This is used by stochastic solvers. The value is for
  /// the compartment and species index provided. For reactions not
  /// connected to a specific species (e.g. mass action reactions) the
  /// species index should be -1.
  ///
  /// @note If reactions are used where this function is not
  /// implemented, the simulation will terminate.
  ///
  /// @return The propensity value
  ///
  */
  ///
  virtual double propensity(Compartment &compartment,size_t species,DataMatrix &y);

  ///
  /// @brief Makes a discrete update of a single reaction.
  ///
  /// The state is updated discretely according to the reaction rule.
  /// The update is for the compartment and species index
  /// provided. For reactions not connected to a specific species
  /// (e.g. mass action reactions) the species index should be -1, and
  /// species are updated according to the rule(s) defined in the
  /// reaction.
  ///
  /// @note If reactions are used where this function is not
  /// implemented, the simulation will terminate.
  ///
  virtual void discreteUpdate(Compartment &compartment,size_t species,DataMatrix &y);

  ///
  /// @brief Prints the reaction
  ///
  /// The print function can be defined if the reaction needs to be
  /// printed. Typically, this is only done for specific reactions which have
  /// parameters as states.
  ///
  /// @see baseReaction::update()
  ///
  virtual void print( std::ofstream &os );

  ///
  /// @brief Prints the reaction in Cambium format
  ///
  /// 
  ///
  /// @see Organism::printCambium()
  ///
  virtual void printCambium( std::ostream &os, size_t varIndex=-1 ) const;

 protected://MODIFICATION AD111018 Changed from private to protected because these were used in derived classes
  
  // A pointer to the organism class to be able to get all possible data fields 
  Organism *organism_;
  
  // The id_ string is an identification for which reaction that are to be
  // created
  std::string id_;
  
  // The parameter_ vector holds all parameters used by the reaction
  std::vector<double> parameter_;           
  
  // String vector for the parameter names
  std::vector<std::string> parameterId_;           
  
  // Vector of vectors holding indices for variables influencing the reaction
  // derivative
  std::vector< std::vector<size_t> > variableIndex_;  
};

inline std::string BaseReaction::id() const 
{
  return id_;
}

inline size_t BaseReaction::numParameter() const 
{
  return parameter_.size();
}

inline size_t BaseReaction::numVariableIndexLevel() const 
{
  return variableIndex_.size();
}

inline size_t BaseReaction::numVariableIndex(size_t level) const 
{
  return variableIndex_[level].size();
}

inline double BaseReaction::parameter(size_t i) const 
{
  return parameter_[i];
}

inline double& BaseReaction::parameterAddress(size_t i) 
{
  return parameter_[i];
}

inline std::string BaseReaction::parameterId(size_t i) const 
{
  return parameterId_[i];
}

inline size_t BaseReaction::variableIndex(size_t i,size_t j) const 
{
  return variableIndex_[i][j];
}

inline void BaseReaction::setOrganism(Organism* value) 
{
  organism_=value;
}

inline void BaseReaction::setId(const std::string &value) 
{
  id_=value;
}

inline void BaseReaction::setParameter(size_t i,double value) 
{
  parameter_[i]=value;
}

inline void BaseReaction::setParameter(std::vector<double> &value) 
{
  parameter_ = value;
}

inline void BaseReaction::setParameterId(size_t i,
																				 const std::string &value) 
{
  parameterId_[i]=value;
}

inline void BaseReaction::setParameterId(std::vector<std::string> &value) 
{
  parameterId_=value;
}

inline void BaseReaction::setVariableIndex(size_t i, size_t j,size_t value) 
{
  variableIndex_[i][j] = value;
}

inline void BaseReaction::
setVariableIndex(size_t i, std::vector<size_t> &value) 
{
  variableIndex_[i] = value;
}

inline void BaseReaction::
setVariableIndex(std::vector< std::vector<size_t> > &value) 
{
  variableIndex_ = value;
}

#endif


