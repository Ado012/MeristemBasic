//
// Filename     : massAction.cc
// Description  : Classes describing mass action reactions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : November 2003
// Revision     : $Id: massAction.cc 664 2016-02-26 17:15:54Z andre $
//
#include"massAction.h"
#include"baseReaction.h"
#include"../organism.h"

namespace MassAction
{
  General::
  General(std::vector<double> &paraValue, 
	  std::vector< std::vector<size_t> > &indValue ) 
  {  
    //
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "MassAction::General::General() "
		<< "Uses only one parameter k_f\n";
      exit(0);
    }
    if( indValue.size() !=2 ) {
      std::cerr << "MassAction::General::General() "
		<< "Two levels of variable indeces are used.\n"
		<< "One for reactants and one for products\n";
      exit(0);
    }
    if( indValue[0].size()<1 ) {
      std::cerr << "MassAction::General::General() "
		<< "If no reactants given, there will be no reaction...\n";
      exit(0);
    }
    //
    // Set the variable values
    //
    setId("massAction::general");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_f";
    setParameterId( tmp );
  }
  
  void General::
  derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
  {  
    double rate = parameter(0);
    if( numVariableIndex(0) )
      for( size_t i=0 ; i< numVariableIndex(0) ; i++ )
	rate *= y[compartment.index()][variableIndex(0,i)];
    else //No reaction defined...
      return;
    
    if (rate<=0.0) //No update
      return;
    
    for( size_t i=0 ; i< numVariableIndex(0) ; ++i )
      dydt[compartment.index()][variableIndex(0,i)] -= rate;
    
    for( size_t i=0 ; i< numVariableIndex(1) ; ++i )
      dydt[compartment.index()][variableIndex(1,i)] += rate;
  }

  void General::
  derivsWithAbs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt, DataMatrix &sdydt) 
  {  
    double rate = parameter(0);
    if( numVariableIndex(0) )
      for( size_t i=0 ; i< numVariableIndex(0) ; i++ )
	rate *= y[compartment.index()][variableIndex(0,i)];
    else //No reaction defined...
      return;
    
    if (rate<=0.0) //No update
      return;
    
    for( size_t i=0 ; i< numVariableIndex(0) ; ++i ) {
      dydt[compartment.index()][variableIndex(0,i)] -= rate;
      sdydt[compartment.index()][variableIndex(0,i)] += std::fabs(rate);
    }
    
    for( size_t i=0 ; i< numVariableIndex(1) ; ++i ) {
      dydt[compartment.index()][variableIndex(1,i)] += rate;
      sdydt[compartment.index()][variableIndex(1,i)] += std::fabs(rate);
    }
  }
  
  //DISABLED 042521
  /*
  size_t General::
  Jacobian(Compartment &compartment,size_t varIndex,DataMatrix &y,JacobianMatrix &A) 
  {  
    double rate = parameter(0);
    if( numVariableIndex(0) )
      for (size_t k=0; k< numVariableIndex(0); ++k)
	rate *= y[compartment.index()][variableIndex(0,k)];
    else //No reaction defined...
      return 0;
    if (rate<=0.0) //No update
      return 1;
    //
    // Update all A
    //
    for (size_t k=0; k<numVariableIndex(0); ++k ) {
      double element = rate / y[compartment.index()][variableIndex(0,k)];
      size_t j = compartment.index()*y[0].size() + variableIndex(0,k); 
      // Reactants
      for (size_t kk=0; kk<numVariableIndex(0); ++kk ) {
	size_t i = compartment.index()*y[0].size() + variableIndex(0,kk); 
	A(i,j) -= element;
      }
      // products
      for (size_t kk=0; kk<numVariableIndex(1); ++kk ) {
	size_t i = compartment.index()*y[0].size() + variableIndex(1,kk); 
	A(i,j) += element;
      }
    }
    return 1;
  }
  */


  double General::
  propensity(Compartment &compartment,size_t varIndex,DataMatrix &y)
  {
    double rate = parameter(0);
    if( numVariableIndex(0) )
      for (size_t k=0; k<numVariableIndex(0); ++k)
	rate *= y[compartment.index()][variableIndex(0,k)];
    else //No reaction defined...
      return 0.0;
    if (rate<=0.0) //No update
      return 0.0;
    return rate;
  }
  
  void General::
  discreteUpdate(Compartment &compartment,size_t varIndex,DataMatrix &y)
  {
    if( numVariableIndex(0) )
      for (size_t k=0; k<numVariableIndex(0); ++k)
	y[compartment.index()][variableIndex(0,k)] -= 1.0;
    if( numVariableIndex(1) )
      for (size_t k=0; k<numVariableIndex(1); ++k)
	y[compartment.index()][variableIndex(1,k)] += 1.0;
  }
  

  
  Enzymatic::
  Enzymatic(std::vector<double> &paraValue, 
	    std::vector< std::vector<size_t> > 
	    &indValue )
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "MassAction::Enzymatic::Enzymatic() "
		<< "Uses only one parameter k_f\n";
      exit(0);
    }
    if( indValue.size() !=3 ) {
      std::cerr << "MassAction::Enzymatic::Enzymatic() "
		<< "Three levels of variable indices are used.\n"
		<< "One for reactants, one for products and one for enzymes\n";
      exit(0);
    }
    if( indValue[0].size()<1 ) {
      std::cerr << "MassAction::Enzymatic::Enzymatic() "
		<< "If no reactants given, there will be no reaction...\n";
      exit(0);
    }
    if( indValue[2].size()<1 ) {
      std::cerr << "MassAction::Enzymatic::Enzymatic() "
		<< "If no enzymes given, there will be no reaction...\n";
      exit(0);
    }
    
    // Set the variable values
    //
    setId("massAction::enzymatic");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_f";
    setParameterId( tmp );
  }
  
  void Enzymatic::
  derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
  {  
    double rate = parameter(0);
    if( numVariableIndex(0) )
      for( size_t i=0 ; i< numVariableIndex(0) ; i++ )
	rate *= y[compartment.index()][variableIndex(0,i)];
    else //No reaction defined...
      return;
    
    if( numVariableIndex(2) )
      for( size_t i=0 ; i< numVariableIndex(2) ; i++ )
	rate *= y[compartment.index()][variableIndex(2,i)];
    else //No reaction defined...
      return;
    
    if (rate<=0.0) //No update
      return;
    
    for( size_t i=0 ; i< numVariableIndex(0) ; i++ )
      dydt[compartment.index()][variableIndex(0,i)] -= rate;
    
    for( size_t i=0 ; i< numVariableIndex(1) ; i++ )
      dydt[compartment.index()][variableIndex(1,i)] += rate;
  }

  void Enzymatic::
  derivsWithAbs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt, DataMatrix &sdydt) 
  {  
    double rate = parameter(0);
    if( numVariableIndex(0) )
      for( size_t i=0 ; i< numVariableIndex(0) ; i++ )
	rate *= y[compartment.index()][variableIndex(0,i)];
    else //No reaction defined...
      return;
    
    if( numVariableIndex(2) )
      for( size_t i=0 ; i< numVariableIndex(2) ; i++ )
	rate *= y[compartment.index()][variableIndex(2,i)];
    else //No reaction defined...
      return;
    
    if (rate<=0.0) //No update
      return;
    
    for( size_t i=0 ; i< numVariableIndex(0) ; i++ ) {
      dydt[compartment.index()][variableIndex(0,i)] -= rate;
      sdydt[compartment.index()][variableIndex(0,i)] += rate;
    }
    
    for( size_t i=0 ; i< numVariableIndex(1) ; i++ ) {
      dydt[compartment.index()][variableIndex(1,i)] += rate;
      sdydt[compartment.index()][variableIndex(1,i)] += rate;
    }
  }
  
  //DISABLED 042521
  /*
  size_t Enzymatic::
  Jacobian(Compartment &compartment,size_t varIndex,DataMatrix &y,JacobianMatrix &A) 
  {  
    double rate = parameter(0);
    if( numVariableIndex(0) )
      for (size_t k=0; k<numVariableIndex(0); ++k)
	rate *= y[compartment.index()][variableIndex(0,k)];
    else //No reaction defined...
      return 0;
    if( numVariableIndex(2) )
      for (size_t k=0 ; k<numVariableIndex(2); ++k)
	rate *= y[compartment.index()][variableIndex(2,k)];
    else //No reaction defined...
      return 0;
    if (rate<=0.0) //No update
      return 1;
    //
    // Update all A
    //
    // Dependence on reactants
    for (size_t k=0; k<numVariableIndex(0); ++k ) {
      double element = rate / y[compartment.index()][variableIndex(0,k)];
      size_t j = compartment.index()*y[0].size() + variableIndex(0,k); 
      // Reactants
      for (size_t kk=0; kk<numVariableIndex(0); ++kk ) {
	size_t i = compartment.index()*y[0].size() + variableIndex(0,kk); 
	A(i,j) -= element;
      }
      // products
      for (size_t kk=0; kk<numVariableIndex(1); ++kk ) {
	size_t i = compartment.index()*y[0].size() + variableIndex(1,kk); 
	A(i,j) += element;
      }
    }
    // Dependence on enzymes
    for (size_t k=0; k<numVariableIndex(2); ++k ) {
      double element = rate / y[compartment.index()][variableIndex(2,k)];
      size_t j = compartment.index()*y[0].size() + variableIndex(2,k); 
      // Reactants
      for (size_t kk=0; kk<numVariableIndex(0); ++kk ) {
	size_t i = compartment.index()*y[0].size() + variableIndex(0,kk); 
	A(i,j) -= element;
      }
      // products
      for (size_t kk=0; kk<numVariableIndex(1); ++kk ) {
	size_t i = compartment.index()*y[0].size() + variableIndex(1,kk); 
	A(i,j) += element;
      }
    }
    return 1;
  }
*/

  void Enzymatic::printCambium( std::ostream &os, size_t varIndex ) const
  {
    std::string varName=organism()->variableId(varIndex);
    std::string varNameCre=organism()->variableId(variableIndex(0,0));
    os << "Arrow[Organism[" << id() << "],{";
    //reactants
    for (size_t i=0; i<numVariableIndex(0); i++) {
      if (i!=0) {
	os << ",";
      }
      std::string varNameCre=organism()->variableId(variableIndex(0,i));
      os << varNameCre << "[i]";
    }
    os << "},{";
    //Enzymes
    for (size_t i=0; i<numVariableIndex(2); i++) {
      if (i!=0) {
	os << ",";
      }
      std::string varNameCre=organism()->variableId(variableIndex(2,i));
      os << varNameCre << "[i]";
    }
    os << "},{";
    // Products
    for (size_t i=0; i<numVariableIndex(1); i++) {
      if (i!=0) {
	os << ",";
      }
      std::string varNameCre=organism()->variableId(variableIndex(2,i));
      os << varNameCre << "[i]";
    }
    os << "},Parameters[";
    os << parameter(0);
    os << "], ParameterNames[" << parameterId(0) << "], VarIndices[{";
    //reactants
    for (size_t i=0; i<numVariableIndex(0); i++) {
      if (i!=0) {
	os << " ";
      }
      os << variableIndex(0,i);
    }
    os << "}{";
    //products
    for (size_t i=0; i<numVariableIndex(1); i++) {
      if (i!=0) {
	os << " ";
      }
      os << variableIndex(1,i);
    }
    os << "}{";
    //enzymes
    for (size_t i=0; i<numVariableIndex(2); i++) {
      if (i!=0) {
	os << " ";
      }
      os << variableIndex(2,i);
    }
    os << "}],"; 
    os << "Solve{";
    //reactants
    for (size_t i=0; i<numVariableIndex(0); i++) {
      std::string varNameCre=organism()->variableId(variableIndex(0,i));
      os << varNameCre << "[i]\' = ";
    }
    //products
    for (size_t i=0; i<numVariableIndex(1); i++) {
      std::string varNameCre=organism()->variableId(variableIndex(0,i));
      os << "-" << varNameCre << "[i]\' = ";
    }
    os << "p_0 ";
    //reactants
    for (size_t i=0; i<numVariableIndex(0); i++) {
      std::string varNameCre=organism()->variableId(variableIndex(0,i));
      os << varNameCre << "[i] ";
    }
    //Enzymes
    for (size_t i=0; i<numVariableIndex(2); i++) {
      std::string varNameCre=organism()->variableId(variableIndex(2,i));
      os << varNameCre << "[i]";
    }
    os << ", "
       << "Using[" << organism()->topology().id() << "::Cell, i]}";
    os << "]" << std::endl;
  }




} //end namespace MassAction
