//
// Filename     : baseReaction.cc
// Description  : A base class describing variable updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Modder       : Albert Do 
// Created      : October 2003
// Modified     : September 2017
// Revision     : $Id: baseReaction.cc 664 2016-02-26 17:15:54Z andre $
//

#include <vector>

#include "../organism.h"
#include "baseReaction.h"
#include <boost/algorithm/string.hpp>

#include "massAction.h"
#include "extendedMeristemReactions.h"
//#include "additionalReactions.h" ADDITION 050917: Trouble adding reactions to seperate file so just merged to
//old code for time being. 

BaseReaction::~BaseReaction()
{
}

BaseReaction *
BaseReaction::createReaction(std::vector<double> &paraValue,
			     std::vector< std::vector<size_t> > &indValue, 
			     const std::string &idValue ) 
{

	std::cerr << idValue << "\n";
	//Grn type of (creation) reactions

    //Extended meristem reactions
    //extendedMeristemReactions.h,extendedMeristemReactions.cc
    if(boost::iequals(idValue,"hill_Weitao"))
		return new Hill_Weitao(paraValue,indValue);    
	else if(boost::iequals(idValue,"hill_Weitao2"))  
		return new Hill_Weitao2(paraValue,indValue);    
	else if(boost::iequals(idValue,"hill_Weitao3"))  
		return new Hill_Weitao3(paraValue,indValue);    
	else if(boost::iequals(idValue,"mutatableNuclearExport"))  
		return new MutatableNuclearExport(paraValue,indValue);    
	else if(boost::iequals(idValue,"spatial"))  
		return new Spatial(paraValue,indValue);    
	else if(boost::iequals(idValue,"hill_Gated"))  
		return new Hill_Gated(paraValue,indValue);    
	else if(boost::iequals(idValue,"multiplicative"))  
		return new Multiplicative(paraValue,indValue);    
	else if(boost::iequals(idValue,"creationTest"))  
		return new CreationTest(paraValue,indValue);  
	else if(boost::iequals(idValue,"degradationOne_alt"))  
		return new DegradationOne_alt(paraValue,indValue);
	else if(boost::iequals(idValue,"diffusionSimple_alt"))  
		return new DiffusionSimple_alt(paraValue,indValue);    
	else if(boost::iequals(idValue,"creationLinear_alt"))  
		return new CreationLinear_alt(paraValue,indValue);  
	else if(boost::iequals(idValue,"creationLimited"))  
		return new CreationLimited(paraValue,indValue);  
	else if(boost::iequals(idValue,"degradationTest"))  
		return new DegradationTest(paraValue,indValue);  
	else if(boost::iequals(idValue,"cLV3_Dynamics"))  
		return new CLV3_Dynamics(paraValue,indValue); 
	else if(boost::iequals(idValue,"crmupdatem"))  
		return new CrmupdateM(paraValue,indValue);   	
	else if(boost::iequals(idValue,"crmupdated"))  
		return new CrmupdateD(paraValue,indValue);   
	else if(boost::iequals(idValue,"wUSRNA_Dynamics"))  
		return new WUSRNA_Dynamics(paraValue,indValue);
	else if(boost::iequals(idValue,"wUSNuc_Dynamics"))  
		return new WUSNuc_Dynamics(paraValue,indValue);
	else if(boost::iequals(idValue,"wUSCyto_Dynamics"))  
		return new WUSCyto_Dynamics(paraValue,indValue);
    else if(boost::iequals(idValue,"CLV3PEPTIDE_Dynamics"))
        return new CLV3PEPTIDE_Dynamics(paraValue,indValue);
    else if(boost::iequals(idValue,"cLV3_Tracker"))
        return new CLV3_Tracker(paraValue,indValue);
    else if(boost::iequals(idValue,"wUSP_ExportRate"))
        return new WUSP_ExportRate(paraValue,indValue);
    else if(boost::iequals(idValue,"ckReceptor_Dynamics"))
        return new CkReceptor_Dynamics(paraValue,indValue);//Added 051020
    else if(boost::iequals(idValue,"ckLigand_Dynamics"))
        return new CkLigand_Dynamics(paraValue,indValue);
    else if(boost::iequals(idValue,"ckComplex_Dynamics"))
        return new CkComplex_Dynamics(paraValue,indValue);

    //Mass action types of reactions
    //massAction.h,massAction.cc
    else if( boost::iequals(idValue,"massAction::Enzymatic") || boost::iequals(idValue,"massActionEnzymatic"))
      return new MassAction::Enzymatic(paraValue,indValue);


	//
	// Transport reactions
	//
	//transport/diffusion.h
	else if(boost::iequals(idValue,"diffusionSimple"))
		return new DiffusionSimple(paraValue,indValue);


	// Default, if nothing found
	//
	else {
	  std::cerr << "\nBaseReaction::createReaction() WARNING: Reactiontype " 
		    << idValue << " not known, no reaction created.\n\7";
	  exit(-1);
	}
}

BaseReaction* 
BaseReaction::createReaction(std::istream &IN ) 
{
	std::string idVal;
  size_t pNum,levelNum;
  IN >> idVal;
  IN >> pNum;
  IN >> levelNum;
  std::vector<size_t> varIndexNum( levelNum );
  for( size_t i=0 ; i<levelNum ; i++ )
    IN >> varIndexNum[i];
  
  std::vector<double> pVal( pNum );
  for( size_t i=0 ; i<pNum ; i++ )
    IN >> pVal[i];
  
  std::vector< std::vector<size_t> > varIndexVal( levelNum );
  for( size_t i=0 ; i<levelNum ; i++ )
    varIndexVal[i].resize( varIndexNum[i] );
  
  for( size_t i=0 ; i<levelNum ; i++ )
    for( size_t j=0 ; j<varIndexNum[i] ; j++ )
      IN >> varIndexVal[i][j];
  
  return createReaction(pVal,varIndexVal,idVal);
}

void BaseReaction::initiate(double t,DataMatrix &y) 
{
}

void BaseReaction::
derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt) 
{
  std::cerr << "BaseReaction::derivs() should not be used. "
	    << "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

void BaseReaction::update(double h, double t,DataMatrix &y) 
{
}

void BaseReaction::
derivsWithAbs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt) 
{
  std::cerr << "BaseReaction::derivsWithAbs() should not be used. "
	    << "Should always be mapped onto one of the real reaction types." << std::endl
	    << "This indicates that the derivsWithAbs is not defined for one of the reactions" << std::endl
	    << " in the model file and that a solver with noise (e.g. heunito) might not be used." << std::endl;
  exit(0);
}  

/*
size_t BaseReaction::
Jacobian(Compartment &compartment,size_t species,DataMatrix &y,JacobianMatrix &A)
{
  std::cerr << "BaseReaction::Jacobian() should not be used. "
	    << "Should always be mapped onto one of the real types." << std::endl
	    << "This indicates that the Jacobian is not defined for one of the reactions" << std::endl
	    << " and that an implicit solver might not be used." << std::endl;
  exit(0);
}
*/

double BaseReaction::
propensity(Compartment &compartment,size_t species,DataMatrix &y)
{
  std::cerr << "BaseReaction::propensity() should not be used. "
	    << "Should always be mapped onto one of the real types." 
	    << std::endl
	    << "This indicates that the propensity is not defined for one of the reactions" 
	    << std::endl
	    << " and that a stochastic solver can not be used." << std::endl;
  exit(0);
}

void BaseReaction::
discreteUpdate(Compartment &compartment,size_t species,DataMatrix &y)
{
  std::cerr << "BaseReaction::discreteUpdate() should not be used. "
	    << "Should always be mapped onto one of the real types." 
	    << std::endl
	    << "This indicates that the propensity is not defined for one of the reactions" 
	    << std::endl
	    << " and that a stochastic solver can not be used." << std::endl;
  exit(0);
}

void BaseReaction::print( std::ofstream &os ) 
{
  std::cerr << "BaseReaction::print(ofstream) should not be used. "
	    << "Should always be mapped onto one of the real types.\n";
  exit(0);
}

void BaseReaction::printCambium( std::ostream &os, size_t varIndex ) const
{
  std::cerr << "BaseReaction::printCambium(ofstream) should not be used. " << std::endl
	    << "Should always be mapped onto one of the real types." << std::endl
	    << "Creating unresolve reaction in output." << std::endl;
  os << "Arrow[Organism[" << id() << "],{},{},{},Parameters[";
  for (size_t i=0; i<numParameter(); ++i) {
    if (i!=0)
      os << ", ";
    os << parameter(i);
  }
  os << "]" << std::endl;
}

Organism* BaseReaction::organism() const 
{
  return organism_;
}

/*size_t BaseReaction::numSpecies(){
	return organism()->numSpecies();
}*/
