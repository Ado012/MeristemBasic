/**
 * Filename     : species.h
 * Description  : A class describing a species (typically a molecule)
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : October 2003
 * Revision     : $Id: species.h 600 2015-01-02 16:45:56Z henrik $
 */
#ifndef SPECIES_H
#define SPECIES_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>

#include"reactions/baseReaction.h"
#include"compartment.h"
#include"common/typedefs.h"

///
/// @brief A class for describing a species (typically a molecule)
///
/// The Species class is describing properties of a
///  molecule(species). An organism have one object for each molecule
///  type. The update rules comes from the baseReactions included.
///
class Species {
  
 private:
  
  std::string id_;
  size_t index_;
  std::vector<BaseReaction*> reaction_;
  
 public:
  
  Species();
  Species( const Species & speciesCopy );
  Species( char *inFile );
  Species( const std::string &inFile );
  Species( std::ifstream &IN );
  
  ~Species();
  
  //Species & operator=( const Species & speciesCopy );
  
  // Get values
  //!Returns the id string
  inline std::string id() const;
  //!Returns the index of the species
  inline size_t index() const;
  //!Returns the number of reactions
  inline size_t numReaction() const;  
  
  //!Returns a reaction (pointer)
  inline BaseReaction* reaction(size_t i) const;
  //!Returns the reaction pointers vector
  inline const std::vector<BaseReaction*> reaction() const;
  
  // Set values
  ///
  /// @brief Sets the id string
  ///
  inline void setId(const std::string &value);
  ///
  /// @brief Sets the index of the species
  ///
  inline void setIndex(size_t value);
  ///
  /// @brief Sets the reaction(pointer) with index i
  ///
  inline void setReaction(size_t i,BaseReaction* value);
  ///
  /// @brief Sets all reaction(pointers) for a species
  ///
  inline void setReaction(const std::vector<BaseReaction*> value);
  
  //Other functions
  ///
  /// @brief Reads a species from an open file stream
  ///
  void readSpecies( std::ifstream &IN );
  ///
  /// @brief Reads and adds a reaction to the species from an open file stream
  /// 
  int addReaction( std::ifstream &IN );
  ///
  /// @brief Loops over species reactions and calls the reaction derivative function
  ///
  void derivs(Compartment &compartment,DataMatrix &y,DataMatrix &dydt);
  ///
  /// @brief Loops over species reactions and calls the reaction derivativeWithAbs function
  ///
  void derivsWithAbs(Compartment &compartment,DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt);
  ///
  /// @brief Loops over species reactions and calls the reaction Jacobian function
  ///
  /// Adds all Jacobian elements from all reactions that belong to a
  ///species. The function returns a nonzero number if future updates are
  ///required.
  ///
  /*
  size_t Jacobian(Compartment &compartment,DataMatrix &y,JacobianMatrix &A);
  */
};

inline std::string Species::id() const {
  return id_;
}

inline size_t Species::index() const {
  return index_;
}

inline size_t Species::numReaction() const {
  return reaction_.size();
}

inline BaseReaction* Species::reaction(size_t i) const {
  return reaction_[i];
}

inline const std::vector<BaseReaction*> Species::reaction() const {
  return reaction_;
}

inline void Species::setId(const std::string &value) {
  id_=value;
}

inline void Species::setIndex( size_t value ) {
  index_ = value;
}

inline void Species::setReaction(size_t i,BaseReaction* value) {
  reaction_[i]=value;
}

inline void Species::setReaction( const std::vector<BaseReaction*> value) {
  reaction_ = value;
}

#endif
