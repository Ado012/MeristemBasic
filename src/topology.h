/**
 * Filename     : topology.h
 * Description  : A class describing the compartment topology 
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : October 2003
 * Revision     : $Id: topology.h 600 2015-01-02 16:45:56Z henrik $
 */
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>

//#include"compartment/baseCompartmentChange.h"
#include"reactions/baseReaction.h"
#include"compartment.h"
#include"common/typedefs.h"

class BaseReaction;

//! A class for describing the updates of an topology (compartment volumes etc)
/*! The Topology class is describing topology properties of a
  compartment e.g the radius if the compartment is spherical. An
  organism have typically one object of the type (if homolog
  compartments are assumed). The update rules come from the included
  baseReactions.
*/ 
class Topology {
  
 private:
  
  std::string id_;
  std::vector<BaseReaction*> reaction_;           
 // std::vector<BaseCompartmentChange*> compartmentChange_;
  unsigned int numVariable_;
  unsigned int numDimension_;
  
  inline void setNumVariable( int value );
  inline void setNumDimension( int value );
  
 public:
  
  Topology();
  Topology( Topology &topologyCopy );
  Topology( char *inFile );
  Topology( const std::string &inFile );
  Topology( std::ifstream &IN );
  
  ~Topology();
  
  // Get values
  inline std::string id() const;
  inline size_t numReaction() const;  
  inline size_t numCompartmentChange() const;  
  inline unsigned int numVariable() const;  
  inline unsigned int numDimension() const;  
  
  inline BaseReaction* reaction(int i) const;
  inline const std::vector<BaseReaction*> reaction() const;
  //inline BaseCompartmentChange* compartmentChange(int i) const;
  //inline const std::vector<BaseCompartmentChange*> compartmentChange() const;
  
  // Set values
  inline void setId(const std::string &value);
  inline void setReaction(int i,BaseReaction* value);
  inline void setReaction(const std::vector<BaseReaction*> value);
 // inline void setCompartmentChange(int i,BaseCompartmentChange* value);
//  inline void setCompartmentChange(const
                 //  std::vector<BaseCompartmentChange*> value);
  
  //Other functions
  void readTopology( std::ifstream &IN);
  int addReaction( std::istream &IN );
  int addCompartmentChange( std::istream &IN );
  ///
  /// @brief Loops over topology reactions and calls the reaction
  /// derivative function
  ///
  void derivs(Compartment &compartment,DataMatrix &y,DataMatrix &dydt);
  ///
  /// @brief Loops over topology reactions and calls the reaction
  /// derivativeWithAbs function
  ///
  void derivsWithAbs(Compartment &compartment,DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt);
  ///
  /// @brief Loops over topology reactions and calls the reaction
  /// Jacobian function
  ///
  /// Adds all Jacobian elements from all reactions that belong to a
  /// topology. The function returns a nonzero number if future updates
  /// are required.
  ///
  /*
  size_t Jacobian(Compartment &compartment,DataMatrix &y,JacobianMatrix &A);
  */
};

//!Returns the id string
inline std::string Topology::id() const {
  return id_;
}

//!Returns the number of reactions
inline size_t Topology::numReaction() const {
  return reaction_.size();
}

/*
//!Returns the number of compartmentChange updates 
inline size_t Topology::numCompartmentChange() const {
  return compartmentChange_.size();
}
*/

//!Returns the number of topology variables
inline unsigned int Topology::numVariable() const {
  return numVariable_;
}

//!Returns the number of dimensions
inline unsigned int Topology::numDimension() const {
  std::cerr << "Topology::numDimension() not safly defined...\n";
  exit(0);
  return numDimension_;
}

//!Returns a reaction (pointer)
inline BaseReaction* Topology::reaction(int i) const {
  return reaction_[i];
}

/*
//!Returns a compartmentChange (pointer)
inline BaseCompartmentChange* Topology::compartmentChange(int i) const {
  return compartmentChange_[i];
}
*/

//!Returns the reaction pointers vector
inline const std::vector<BaseReaction*> Topology::reaction() const {
  return reaction_;
}

/*
//!Returns the compartmentChange pointers vector
inline const std::vector<BaseCompartmentChange*> 
Topology::compartmentChange() const {
  return compartmentChange_;
}
*/

//!Sets the id string
inline void Topology::setId(const std::string &value) {
  id_=value;
}

//!Sets the reaction(pointer) with index i
inline void Topology::setReaction(int i,BaseReaction* value) {
  reaction_[i]=value;
}

//!Sets all reaction(pointers) for a topology
inline void Topology::setReaction(const std::vector<BaseReaction*> value) {
  reaction_ = value;
}

/*
//!Sets the compartmentChange(pointer) with index i
inline void Topology::
setCompartmentChange(int i, BaseCompartmentChange* value) {
  compartmentChange_[i]=value;
}

//!Sets all reaction(pointers) for a topology
inline void Topology::
setCompartmentChange(const std::vector<BaseCompartmentChange*> value) {
  compartmentChange_ = value;
}
*/
//!Sets the number of topology variables
inline void Topology::setNumVariable(int value ) {
  numVariable_ = value;
}

//!Sets the number of dimensions
inline void Topology::setNumDimension(int value ) {
  numDimension_ = value;
}

#endif
