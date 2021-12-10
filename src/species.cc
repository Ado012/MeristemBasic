/**
 * Filename     : species.cc
 * Description  : A class describing a species (typically a molecule)
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : October 2003
 * Revision     : $Id: species.cc 600 2015-01-02 16:45:56Z henrik $
 */

#include<iostream>
#include<fstream>
#include "species.h"

Species::Species() 
{
  //std::cerr << "Species::Species() not defined.\n";
  setIndex(static_cast<size_t>(-1));
}

Species::Species( const Species & speciesCopy ) 
{   
  setId( speciesCopy.id() );
  setIndex( speciesCopy.index() );
  setReaction( speciesCopy.reaction() );
}

Species::Species( char *inFile ) 
{  
  std::ifstream IN( inFile );
  if( !IN ) {
    std::cerr << "Species::Species(char*) - "
	      << "Cannot open file " << inFile << "\n\n\7";exit(-1);}
  readSpecies(IN);
}

Species::Species( const std::string &inFile ) 
{   
  const char *tmp = inFile.c_str();
  std::ifstream IN( tmp );
  if( !IN ) {
    std::cerr << "Species::Species() - "
	      << "Cannot open file " << inFile << "\n\n\7";exit(-1);}
  readSpecies(IN);
}

Species::Species( std::ifstream &IN ) 
{  
  readSpecies(IN);
}

Species::~Species()
{
}

void Species::readSpecies( std::ifstream &IN ) 
{  
  //Remove all reactions before adding
  if( numReaction() ) 
    reaction_.resize(0);
  	
  size_t indexVal,numReactionVal;
  std::string idVal;
  IN >> idVal;
  setId(idVal);
  IN >> indexVal;
  setIndex(indexVal);
  IN >> numReactionVal;
  
  
  for( size_t i=0 ; i<numReactionVal ; i++ ) {
    if( addReaction(IN) )
      std::cerr << "Species::readSpecies() Warning Adding reaction failed for "
		<< "species " << index() << " " << id() << "\n";
  }
}

int Species::addReaction( std::ifstream &IN ) 
{
  if( !IN )
    return -1;
  reaction_.push_back( BaseReaction::createReaction(IN) );
  return 0;
}

void Species::derivs(Compartment &compartment,DataMatrix &y,DataMatrix &dydt)
{
  for( size_t i=0 ; i<numReaction() ; i++ )
    reaction(i)->derivs(compartment,index(),y,dydt);

}

void Species::derivsWithAbs(Compartment &compartment,DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt)
{
  for( size_t i=0 ; i<numReaction() ; i++ )
    reaction(i)->derivsWithAbs(compartment,index(),y,dydt,sdydt);
}

/*
size_t Species::Jacobian(Compartment &compartment,DataMatrix &y,JacobianMatrix &A)
{
  size_t updateFlag=0;
  for( size_t i=0 ; i<numReaction() ; i++ )
    updateFlag += reaction(i)->Jacobian(compartment,index(),y,A);
  return updateFlag;
}
*/
