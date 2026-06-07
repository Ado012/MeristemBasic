/**
 * Filename     : topology.cc
 * Description  : A class describing a topology (size,volume updates)
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : October 2003
 * Revision     : $Id: topology.cc 600 2015-01-02 16:45:56Z henrik $
 */

#include<iostream>
#include<fstream>
#include "topology.h"

Topology::Topology() {
  //std::cerr << "Topology::Topology().\n";
  setNumVariable(0);
  setNumDimension(0);
}

Topology::Topology( Topology & topologyCopy ) { 

  //std::cerr << "Topology::Topology(Topology).\n";
  setId( topologyCopy.id() );
  setNumVariable( topologyCopy.numVariable() );
  setReaction( topologyCopy.reaction() );
//  setCompartmentChange( topologyCopy.compartmentChange() );
}

Topology::Topology( char *inFile ) {
  
  std::ifstream IN( inFile );
  if( !IN ) {
    std::cerr << "Topology::Topology() - "
	      << "Cannot open file " << inFile << "\n\n\7";exit(-1);}
  readTopology(IN);
}

Topology::Topology( const std::string &inFile ) { 
  
  const char *tmp = inFile.c_str();
  std::ifstream IN( tmp );
  if( !IN ) {
    std::cerr << "Topology::Topology() - "
	      << "Cannot open file " << inFile << "\n\n\7";exit(-1);}
  readTopology(IN);
}

Topology::Topology( std::ifstream &IN ) {
  
  readTopology(IN);
}

// Topology::Topology(std::vector<BaseReaction*> &reactionValue, 
// 		   std::string idValue ) {
//   std::cerr << "Topology::Topology(std::vector<BaseReaction*> "
// 	    << "&reactionValue,... not defined.\n";
// }

Topology::~Topology() {
  //std::cerr << "Topology::~Topology() not defined.\n";
}

//!Read a toppology from an open filestream
void Topology::readTopology( std::ifstream &IN ) {
  
  //Remove all reactions and compartmentChanges before adding
  reaction_.resize(0);
 // compartmentChange_.resize(0);
  
  int numVariableVal,numReactionVal,numCompartmentChangeVal;
  std::string idVal;
  IN >> idVal;
  setId(idVal);
  IN >> numVariableVal;
  setNumVariable( numVariableVal );
  IN >> numReactionVal;
  IN >> numCompartmentChangeVal;
  
	for( int i=0 ; i<numReactionVal ; i++ ) 
		if( addReaction(IN))
      std::cerr << "Topology::readTopology() Warning Adding reaction failed "
								<< "for topology " << id() << "\n\7";
	
    //DISABLED 042521
    /*
	for( int i=0 ; i<numCompartmentChangeVal ; i++ ) 
		if (addCompartmentChange(IN) )
      std::cerr << "Topology::readTopology() Warning Adding compartmentChange"
								<< " failed "
								<< "for topology " << id() << "\n\7";  
                                */
}

//!Adds a reaction to the list from an open file
int Topology::addReaction( std::istream &IN ) {
  if( !IN )
    return -1;
  reaction_.push_back( BaseReaction::createReaction(IN) );
  return 0;
}

/*
int Topology::addCompartmentChange( std::istream &IN ) {
	
	if( !IN ) {
		std::cerr << "Topology::addCompartmentChange\n";
		return -1;
	}
  compartmentChange_.push_back( BaseCompartmentChange::
				createCompartmentChange(IN) );
  return 0;
}
*/


void Topology::derivs(Compartment &compartment,DataMatrix &y,DataMatrix &dydt ) 
{  
  for( size_t i=0 ; i<numReaction() ; i++ )
    reaction(i)->derivs(compartment,compartment.topologyStart(),y,dydt);
}

void Topology::derivsWithAbs(Compartment &compartment,DataMatrix &y,DataMatrix &dydt, DataMatrix &sdydt) 
{  
  for( size_t i=0 ; i<numReaction() ; i++ )
    reaction(i)->derivsWithAbs(compartment,compartment.topologyStart(),y,dydt,sdydt);
}

/*
size_t Topology::Jacobian(Compartment &compartment,DataMatrix &y,JacobianMatrix &A) 
{  
  size_t updateFlag=0;
  for( size_t i=0 ; i<numReaction() ; i++ )
    updateFlag += reaction(i)->Jacobian(compartment,compartment.topologyStart(),y,A);
  return updateFlag;
}
*/
