//
// Filename     : extendedMeristemReactions.h
// Description  : reactions to extend the meristem reaction network
// Author(s)    : Al Do (ado012@ucr.edu)
// Created      : June 2017
// Revision     :
//
#ifndef EXTENDEDMERISTEMREACTIONS_H
#define EXTENDEDMERISTEMREACTIONS_H

#include <cmath>
#include <algorithm>

#include "baseReaction.h"

///
/// @brief The class Grn use a neural network inspired mathematics for gene
/// regulation
///
/// This class uses a neural network inspired update with a sigmoidal function
/// for gene regulation as defined in Mjolsness et al (1991). The update is
/// done according to:
///
/// @f[ \frac{dy_{ij}}{dt} = \frac{1}{p_1}g\left( p_0 + \sum_k p_k y_{ik} \right) @f]
///
/// where g() is a sigmoidal function, the @f$p_k@f$ are the parameters given
/// from position 2 and forward, and the k indices for @f$y_{ik}@f$ are given
/// at the first level of variableIndex. Hence the number of parameters must
/// be two more than the number of indices given at the first level.
///
/// @see Grn::sigmoid(double x) for implementation of sigmoid function.
///
//ADDITION 051017
class crmHill : public BaseReaction {
  
 public:
  float indContrib;
  float chanceBind=0.5;
  float chanceUnbind=0.5;
  float k_bind=0.5;
  float geneCRMSiteChance_Bind; 
  float geneCRMSiteChance_Unbind[5]={0.5,0.5,0.5,0.5,0.5};
  float geneCRMSiteBindBaseChance[5]={0.5,0.5,0.5,0.5,0.5};
  int crmBindingSites=5;
  
  float randvalue, randvalue2, randvalue3;
  
  std::ofstream crmfile;//ADDITION 051517
int crmcheck; //ADDITION 051517
  
  
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
  crmHill(std::vector<double> &paraValue, 
	std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);

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
    /*
  size_t Jacobian(Compartment &compartment,size_t species,DataMatrix &y,JacobianMatrix &A);
*/
  ///
  /// @brief Propensity calculation (probability to happen in [t,t+dt]) for this reaction class
  ///
  /// @see BaseReaction::propensity(Compartment &compartment,size_t species,DataMatrix &y)
  ///
  double propensity(Compartment &compartment,size_t species,DataMatrix &y);

  ///
  /// @brief Discrete update for this reaction used in stochastic simulations
  ///
  /// @see BaseReaction::discreteUpdate(Compartment &compartment,size_t species,DataMatrix &y)
  ///
  void discreteUpdate(Compartment &compartment,size_t species,DataMatrix &y);

  ///
  /// @brief Prints the reaction in Cambium format
  ///
  /// @see BaseReaction::printCambium()
  ///
};

///
/// @brief This class describes a Hill-type production with activators,
///  repressors and one variable restricting @f$V_{max}@f$
///
/// This reaction is a Michaelis-Menten formalism for production with a
/// restricting variable given by
///
/// @f[ \frac{dy_{ij}}{dt} = p_0\cdot y_{Restricting} \frac{y_{ik}^{p_2}}{p_1^{p_2}+y_{ik}^{p_2}}...\frac{p_{1'}^{p_{2'}}}{p_{1'}^{p_{2'}}+y_{ik'}^{p_{2'}}}@f]
///
/// where @f$p_0 \cdot y_{Restricting}@f$ is the maximal rate (@f$V_{max}@f$), and @f$p_1@f$ is the Hill
/// constant (@f$K_{half}@f$), and @f$p_2@f$ is the Hill coefficient (n). The
/// k index is given in the first level of varIndex and corresponds to the
/// activators, and the k' in the second level which corresponds to the
/// repressors.
///
///	There are three variable layers (hillRepressed x 3 x x). @f$y_{Restricting}@f$ is given by the first layer
///
///	The second layer is for activators and the third for repressors. In both of these layers each index corresponds to a pair of parameters
/// (K,n) that are preceded by the @f$p_0@f$ parameter.
///
/// In the model file the reaction will be defined (somewhat different for different number of
/// activators and repressors) e.g. as
///
/// One activator:
/// @verbatim
/// hill 3 3 1 1 0
/// V K n
/// restrictingVariable_index
/// activator_index
/// @endverbatim
///
/// One repressor:
/// @verbatim
/// hillRestricted 3 3 1 0 1
/// V K n
/// restrictingVariable_index
/// repressor_index
/// @endverbatim
///
/// Two activators and one repressor:
/// @verbatim
/// hillrestricted 5 3 1 2 1
/// V K_A1 n_A1 K_A2 n_A2 K_R n_R
/// restrictingVariable_index
/// A1_index A2_index R_index
/// @endverbatim
///
///	The original use for this reaction was to be able to use a dummy variable in order to smoothly regulate the
///	production of a certain compound during optimizations.
///

///
/// @brief Implements a degradation proportional to its own concentration and
/// two additional variables.
///
/// The update to the variable is
///
/// @f[ \frac{dc}{dt} -= k_d c C_1 C_2 @f]
///
/// where @f$ k_d @f$ is the degradation rate, and @f$ C_1,C_2 @f$ are the
/// user supplied variable.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// spatial 1 1 2
/// k_d
/// index1 index2
/// @endverbatim
///
class Spatial : public BaseReaction {
  
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
  Spatial(std::vector<double> &paraValue, 
									 std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class Hill_Weitao : public BaseReaction {
  
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
  Hill_Weitao(std::vector<double> &paraValue, 
									 std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};



class Hill_Weitao2 : public BaseReaction {
  
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
  Hill_Weitao2(std::vector<double> &paraValue, 
									 std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};

class Hill_Weitao3 : public BaseReaction {
  
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
  Hill_Weitao3(std::vector<double> &paraValue, 
									 std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};




class Hill_Gated : public BaseReaction {
  
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
  Hill_Gated(std::vector<double> &paraValue, 
									 std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class MutatableNuclearExport : public BaseReaction {
  
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
  MutatableNuclearExport(std::vector<double> &paraValue, 
									 std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class Multiplicative : public BaseReaction {
  
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
  Multiplicative(std::vector<double> &paraValue, 
									 std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};

class CreationTest : public BaseReaction {
  
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
  /// @see BaseReaction::createReaction(std::vector<double>&,...)
  ///
  CreationTest(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
  
  ///
  /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
  ///
  /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
  ///
 
};

class DegradationTest : public BaseReaction {
  
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
  /// @see BaseReaction::createReaction(std::vector<double>&,...)
  ///
  DegradationTest(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
  
  ///
  /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
  ///
  /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
  ///
 
};

class CreationLimited : public BaseReaction {
  
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
  /// @see BaseReaction::createReaction(std::vector<double>&,...)
  ///
  CreationLimited(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
  
  ///
  /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
  ///
  /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
  ///
 
};



class DegradationOne_alt : public BaseReaction {
  
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
  DegradationOne_alt(std::vector<double> &paraValue, 
								 std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);


};

class DiffusionSimple_alt : public BaseReaction {
  
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
  DiffusionSimple_alt(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);

  ///
  /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
  ///
  /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
  ///

};

class CreationLinear_alt : public BaseReaction {
  
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
  /// @see BaseReaction::createReaction(std::vector<double>&,...)
  ///
  CreationLinear_alt(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt );

  ///
  /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
  ///
  /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
  ///
};

class CLV3_Dynamics : public BaseReaction {
  
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
  CLV3_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class CrmupdateM : public BaseReaction {
  
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
  CrmupdateM(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class CrmupdateD : public BaseReaction {
  
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
  CrmupdateD(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class WUSRNA_Dynamics : public BaseReaction {
  
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
  WUSRNA_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};

class WUSNuc_Dynamics : public BaseReaction {
  
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
  WUSNuc_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};

class WUSCyto_Dynamics : public BaseReaction {
  
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
  WUSCyto_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
	///
	/// @brief Derivative function for this reaction class
	///
	/// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
	///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class CLV3PEPTIDE_Dynamics : public BaseReaction {

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
  CLV3PEPTIDE_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class CLV3_Tracker : public BaseReaction {

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
  CLV3_Tracker(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class CLV3Marker_Dynamics : public BaseReaction {

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
  CLV3Marker_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};

class WUSP_ExportRate : public BaseReaction {

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
  WUSP_ExportRate(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class CkLigand_Dynamics : public BaseReaction {

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
   CkLigand_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};

class CkReceptor_Dynamics : public BaseReaction {

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
  CkReceptor_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


class CkComplex_Dynamics : public BaseReaction {

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
  CkComplex_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);
};


///
/// @brief The class DiffusionSimple represents the simplest diffusion update
/// NOT including topological properties
///
/// Generates a diffusion update where NO topological properties such as
/// volumes are taken into account according to:
///
/// @f[ \frac{dy_{ij}}{dt} = D (y_{nj}-y_{ij}) @f]
///
/// and
///
/// @f[ \frac{dy_{nj}}{dt} = D (y_{nj}-y_{ij}) @f]
///
/// where \e j is the molecular index. The update only applies if the neighbor
/// index \e n is larger than \e i.
///
class DiffusionSimple : public BaseReaction {

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
  DiffusionSimple(std::vector<double> &paraValue,
          std::vector< std::vector<size_t> > &indValue );

  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt);

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

  ///
  double propensity(Compartment &compartment,size_t species,DataMatrix &y);

  ///
  /// @brief Discrete update for this reaction used in stochastic simulations
  ///
  /// @see BaseReaction::discreteUpdate(Compartment &compartment,size_t species,DataMatrix &y)
  ///
  void discreteUpdate(Compartment &compartment,size_t species,DataMatrix &y);

  */
  ///
  /// @brief Prints the reaction in Cambium format
  ///
  /// @see BaseReaction::printCambium()
  ///
  void printCambium( std::ostream &os, size_t varIndex ) const;
};











#endif
