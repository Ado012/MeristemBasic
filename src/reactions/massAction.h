//
// Filename     : massAction.h
// Description  : Classes describing mass action reactions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : November 2003
// Revision     : $Id: massAction.h 664 2016-02-26 17:15:54Z andre $
//
#ifndef MASSACTION_H
#define MASSACTION_H

#include "baseReaction.h"

///
/// @brief A collection of one way mass action kinetics reactions, including special enzyme 
/// and Michaelis Menten type of reactions.
///
/// These reactions follows mass action kinetics, i.e. the reaction rate is proportional to the reactants
/// concentrations. The general reaction handles any number of reactants and products, while specialized
/// versions a restricted to specific numbers of reactants/products. These specialized versions have been 
/// implemented to handle combinatoric factors in stochastic simulations where more than one reactant are of
/// the same species. The keyword Enzymatic indicates the possibility to have molecular species increase the 
/// reaction rate proportionaly to their concentrations, but that are not updated (exists before and after the
/// reaction). MM implements a Michaelis-Menten type of dynamics for both the loss of the reactant and increase of
/// the product.
///
namespace MassAction 
{
  ///
  /// @brief A one way mass action reaction
  ///
  /// The reactant indeces are in level zero and products in level one in
  /// variableIndex. One parameter (rate constant) is needed.
  ///
  /// Reactants [R] and products [P] stored in variableIndexLevel 0/1 are
  /// updated. The parameter is k_f in 
  ///
  /// @f[ \frac{d[P]}{dt} = k_f * \prod [R] @f]
  /// @f[ \frac{d[R]}{dt} = - k_f \prod [R] @f]
  ///
  /// In a model file the reaction is defined as:
  ///
  /// @verbatim
  /// massAction 1 2 N_R N_P
  /// k_f
  /// R_1 ... R_{N_R}
  /// P_1 ... P_{N_P}
  /// @endverbatim
  ///
  /// @note In the model file also massAction::general 1 2 ... can be used.
  /// @note In the propensity (for stochastic simulations) combanitoric factors due to same reactant species are not handled.
  /// @note (for developers) )varIndex is not used.
  ///
  class General : public BaseReaction {
    
  public:
    
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    General(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > &indValue );
    
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt );

    ///
    /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
    ///
    /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
    ///
    void derivsWithAbs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt);

    ///
    /// @brief Jacobian function for this reaction class
    ///
    /// @see BaseReaction::Jacobian(Compartment &compartment,size_t species,...)
    ///
    //DISABLED 042521
    /*
    size_t Jacobian(Compartment &compartment,size_t species,DataMatrix &y,JacobianMatrix &A);
    
    ///
    /// @brief Propensity calculation (probability to happen in [t,t+dt]) for this reaction class
    ///
    /// @see BaseReaction::propensity(Compartment &compartment,size_t species,DataMatrix &y)
    ///
    /// @note This function is under construction and currently only works for the case
    /// when the reactants are of different kinds.
    */
    ///
    ///
    double propensity(Compartment &compartment,size_t species,DataMatrix &y);
    
    ///
    /// @brief Discrete update for this reaction used in stochastic simulations
    ///
    /// @see BaseReaction::discreteUpdate(Compartment &compartment,size_t species,DataMatrix &y)
    ///
    void discreteUpdate(Compartment &compartment,size_t species,DataMatrix &y);
  };
  

  ///
  /// @brief A one way mass action reaction with enzyme contribution
  ///
  /// The reactant indeces are in level zero, products in level one, and enzymes
  /// in level two in variableIndex. One parameter (rate constant) is needed.
  ///
  /// Reactants [R], products [P] and enzymes [E] are stored in
  /// variableIndexLevel 0/1/2. [R] and [P] are updated. The parameter is k_f in
  /// 
  /// @f[ \frac{d[P]}{dt} = k_f \prod [E] \prod [R] @f]
  /// @f[ \frac{d[R]}{dt} = - k_f \prod [E] \prod [R] @f]
  /// @f[ \frac{d[E]}{dt} = 0 @f]
  ///
  /// In a model file the reaction is defined as:
  ///
  /// @verbatim
  /// massAction::Enzymatic 1 3 N_R N_P N_E
  /// k_f
  /// R_1 ... R_{N_R}
  /// P_1 ... P_{N_P}
  /// E_1 ... E_{N_E}
  /// @endverbatim
  ///
  /// @note In the model file also massActionEnzymatic 1 3... can be used.
  /// @note Can be used in stochastic simulations.
  /// @note (for developers) varIndex is not used.
  /// 
  class Enzymatic : public BaseReaction {
    
  public:
    
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    Enzymatic(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > &indValue );
    
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt );

    ///
    /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
    ///
    /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
    ///
    void derivsWithAbs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt);
    
    ///
    /// @brief Jacobian function for this reaction class
    ///
    /// @see BaseReaction::Jacobian(Compartment &compartment,size_t species,...)
    ///
    ///
    //DISABLED 042521
    //size_t Jacobian(Compartment &compartment,size_t species,DataMatrix &y,JacobianMatrix &A);
    ///
    /// @brief Prints the reaction in Cambium format
    ///
    /// @see BaseReaction::printCambium()
    ///
    void printCambium( std::ostream &os, size_t varIndex ) const;
  };
  


  
} //end namespace MassAction
#endif

