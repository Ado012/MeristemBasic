//
// Filename     : organism.cc
// Description  : A class describing an organism 
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2003
// Revision     : $Id: organism.cc 646 2015-08-23 18:35:06Z henrik $
//

#include<cmath>
#include <cstring>

#include "common/myFiles.h"
#include"common/myRandom.h"
#include"organism.h"
//#include"val.h"
#include <iostream>

Organism::Organism() 
{
  setNumTopology(0);
  setNumNeighborhood(0);
}

//Organism::Organism( const Organism & organismCopy ) 
//{  
//std::cerr << "Organism::Organism(const Organism &)"
//    << " Not implemented yet...\n";
//exit(0);
//}

Organism::Organism( const char *modelFile, const char *initFile, 
		    int verbose ) 
{
  readOrganism(modelFile,initFile,verbose);
}

Organism::Organism( const std::ifstream& IN, const char *initFile, int verbose){
	readOrganism((std::ifstream &) IN, initFile,verbose);
}

Organism::Organism( const std::string &modelFile, 
		    const std::string &initFile, 
		    int verbose ) 
{  
    //Addition: 051121 Limitation: based on parameter file name
     int namelength = modelFile.length(); //determine input length
     int namegrabposition = namelength - 18; //determine which portion to grab by account for .model and filename

    std::string modelName = modelFile.substr(namegrabposition,12);
  filename_ = modelName;

  readOrganism(modelFile,initFile,verbose);
}

Organism::~Organism() 
{
  // Delete reactions.
  std::vector<BaseReaction *>::iterator i;
  for (i = reaction_.begin(); i != reaction_.end(); ++i)
    delete *i;
  
  // Delete neighborhood.
  if( numNeighborhood() )
    delete neighborhood_;
}

void Organism::readOrganism(std::ifstream & IN, const char *initFile, int verbose )
{
 //Read topology, species, and reaction models from file
  if( verbose )
    std::cerr << "Organism::readOrganism(char*,char*) reading the model part of the file";
  readModel(IN);
  if( verbose )
    std::cerr << "Model Part Read!\n";
  
  //Read initial values from file and create compartments
  if( std::strcmp(initFile,"") ) {
    if( verbose )
      std::cerr << "Organism::readOrganism(char*,char*) reading init file " 
		<< initFile << ".\n";
    readInit(initFile);
    if( verbose )
      std::cerr << "Init file read!\n";
  }
  else {
    std::cerr << "Organism::readOrganism( char*, char* ) "
	      << "Warning! No initial configuration read.\n";
  }
	
}
  

void Organism::readOrganism( const char *modelFile, const char *initFile, 
			     int verbose ) 
{  
  //Read topology, species, and reaction models from file
  if( verbose )
    std::cerr << "Organism::readOrganism(char*,char*) reading model file " 
	      << modelFile << "\n";
  readModel(modelFile);
  if( verbose )
    std::cerr << "Model File Read!\n";
  
  //Read initial values from file and create compartments
  if( std::strcmp(initFile,"") ) {
    if( verbose )
      std::cerr << "Organism::readOrganism(char*,char*) reading init file " 
		<< initFile << ".\n";
    readInit(initFile);
    if( verbose )
      std::cerr << "Init file read!\n";
  }
  else {
    std::cerr << "Organism::readOrganism( char*, char* ) "
	      << "Warning! No initial configuration read.\n";
  }
}

void Organism::readOrganism( const std::string &modelFile, 
			     const std::string &initFile, 
			     int verbose ) 
{
  if( verbose )
    std::cerr << "Organism::readOrganism( std::string,std::string,int)\n";
  const char* mFile=modelFile.c_str();
  const char* iFile=initFile.c_str();
  readOrganism(mFile,iFile,verbose);
}

void Organism::readModel(std::ifstream &IN) {
  
  unsigned int numTopologyVal,numSpeciesVal,numReactionVal,
    numNeighborhoodVal;
  std::string idVal;
  IN >> idVal;
  IN >> numTopologyVal;
  IN >> numSpeciesVal;
  IN >> numReactionVal;
  IN >> numNeighborhoodVal;
  
  if (numTopologyVal != 0 && numTopologyVal != 1) {
    std::cerr << "Organism::readModel(ifstream) - "
	      << "Zero or one topology is allopwed for right now! " 
	      << numTopologyVal << ".\n\n\7";exit(-1);}
  if( numNeighborhoodVal != 0 && numNeighborhoodVal != 1 ) {
    std::cerr << "Organism::readModel(ifstream) - "
	      << "Zero or one neighborhood is allowed for right now! " 
	      << numNeighborhoodVal << ".\n\n\7";exit(-1);}
  setId( idVal );
  setNumTopology( numTopologyVal );
  setNumNeighborhood( numNeighborhoodVal );
  //
  // Read topology
  //
  if( numTopology() ) {
    std::cerr << "org:readModel topology_.readTopology(IN)\n";
    topology_.readTopology(IN);
    std::cerr << "Topology " << topology().id() << " added with "
	      << topology().numVariable() << " variables, " 
	      << topology().numReaction() << " reactions and "
          << "DISABLED 042521 topology().numCompartmentChange()"
	      << " compartment changes\n";
  }
  //
  // Read Species
  //
  //Remove any present species before adding
  if( numSpecies() ) 
    species_.resize(0);
  
  for( unsigned int i=0 ; i<numSpeciesVal ; i++ ) {
    addSpecies(IN);  
    std::cerr << "Species " << species(i).id() << " added with index "
	      << species(i).index() << " and " << species(i).numReaction()
	      << " reactions.\n"; 
    // Also add organism pointer to the reactions.
    // Should preferable be done further down in the hiarchy (as when adding reactions below)
    // but organsim lost at species level...
    for (size_t k=0; k<species(i).numReaction(); ++k) {
      species(i).reaction(k)->setOrganism(this);
    }
  }
  //
  // Read Reactions
  //
  //Remove any present reactions before adding
  if( numReaction() ) 
    reaction_.resize(0);
  
  
  for (size_t i=0; i<numReactionVal; ++i) {
    if (addReaction(IN)) {
      std::cerr << "Organism::ReadModel(ifstream) "
		<< "Warning Adding reaction failed for "
		<< "organism " << id() << " (index " << i << ")" << std::endl;
    }
    else {
      std::cerr << "Reaction " << reaction(numReaction()-1)->id() << " added." << std::endl;
    }
  }
  //
  // Read Neighborhood
  //
  if( numNeighborhood() )
    if( addNeighborhood(IN) )
      std::cerr << "Organism::ReadModel(ifstream) "
		<< "Warning Adding neighborhood failed for "
		<< "organism " << id() << ".\n";
}

void Organism::readModel(const char *fileName) 
{
  std::string tmp(fileName);
  readModel(tmp);
}

void Organism::readModel(const std::string &fileName) 
{
  std::istream *IN = myFiles::openFile(fileName);
  if( !IN ) {
    std::cerr << "Organism::readModel(char*) - "
	      << "Cannot open file " << fileName << "\n\n\7";exit(-1);}
  readModel((std::ifstream &) *IN);
  delete IN;
}

void Organism::readInit(std::ifstream &IN) 
{  
  //Remove all compartments before adding
  if( numCompartment() )
    compartment_.resize(0);
  
  size_t numCompartmentVal,numVariableVal;
  size_t numModelVariable=numTopologyVariable()+numSpecies();
  IN >> numCompartmentVal;
  IN >> numVariableVal;
  if (numVariableVal!=numModelVariable) {
    std::cerr << "Organism::readInit(std::ifstream ) "
	      << "Warning! Different number of variables in model and init file.\n";
    std::cerr << numVariableVal << " in init file and " 
	      << numTopologyVariable() <<  "+" << numSpecies() 
	      << " in organism.\n";
    if (numVariableVal>numModelVariable) {
      std::cerr << "Disregarding final " << numVariableVal-numModelVariable 
		<< " columns in init file." << std::endl;
		}
    if (numVariableVal<numModelVariable) {
      std::cerr << "Adding " << numModelVariable-numVariableVal 
		<< " zeros as initial values." << std::endl;
    }
  }
  std::vector<double> varTmp( numModelVariable );
  for( size_t i=0 ; i<numCompartmentVal ; i++ ) {
    if (numModelVariable==numVariableVal) {
      for( size_t j=0 ; j<numVariableVal ; j++ ) {
	IN >> varTmp[j];
      }
    }
    else if (numModelVariable>numVariableVal) {
      for( size_t j=0 ; j<numVariableVal ; j++ ) {
	IN >> varTmp[j];
      }
      // Add additional zeros
      for( size_t j=numVariableVal ; j<numModelVariable ; j++ ) {
	varTmp[j] = 0.0;
      }
    }
    else if (numModelVariable<numVariableVal) {
      for( size_t j=0 ; j<numModelVariable ; j++ ) {
	IN >> varTmp[j];
      }
      // Read (without saving) the rest of the values in the init file
      double dTmp;
      for( size_t j=numModelVariable ; j<numVariableVal ; j++ ) {
	IN >> dTmp;
      }
    }			
    Compartment tmpCompartment(this,topologyPointer(), varTmp, i, "");
    addCompartment( tmpCompartment );
  }
}

void Organism::readInit(const char *fileName) 
{
  std::string tmp(fileName);
  readInit(tmp);
}

void Organism::readInit(const std::string &fileName) 
{
  std::istream *IN = myFiles::openFile(fileName);
  if( !IN ) {
    std::cerr << "Organism::readInit(std::string) - "
	      << "Cannot open file " << fileName << "\n\n\7";exit(-1);}
  readInit((std::ifstream &) *IN);
  delete IN;
}

void Organism::setInit(const DataMatrix &input) 
{  
  //Remove all compartments before adding
  if( numCompartment() )
    compartment_.resize(0);
  
  size_t numCompartmentVal,numVariableVal;
  numCompartmentVal = input.size();
  numVariableVal=input[0].size();
  if( numVariableVal != numTopologyVariable()+numSpecies() ) {
    std::cerr << "Organism::setInit(std::ifstream ) "
	      << "Wrong number of variables in init file.\n";
    std::cerr << numVariableVal << " in init file and " 
	      << numTopologyVariable() <<  "+" << numSpecies() 
	      << " in organism.\n";
    exit(-1);
  }
  
  for( size_t i=0 ; i<numCompartmentVal ; i++ ) {
    Compartment tmpCompartment(this,topologyPointer(), input[i], i, "");
    addCompartment( tmpCompartment );
  }
}

int Organism::addReaction( std::istream &IN ) 
{
  if( !IN ) return -1;
  reaction_.push_back( BaseReaction::createReaction(IN) );
  reaction_[reaction_.size()-1]->setOrganism(this);
  return 0;
}

std::string Organism::variableId(size_t varIndex) const
{
  if (varIndex<numTopologyVariable()) {
    // Quick fix before topology Id's set name convention
    if (varIndex==numTopologyVariable()-1) {
      return "cellSize";
    }
    else {
      ostringstream out;
      out << "cellPos_";
      out << varIndex;
      return out.str();
    }
  }
  else if (varIndex<numVariable()) {
    return species(varIndex-numTopologyVariable()).id();
  }
  else {
    return "varIdNull";
  }
}
    

void Organism::initiateParameter() 
{  
  parameter_.resize(0);
  
  if( numTopology() ) {
    //collect topology reaction parameters
    for( size_t i=0 ; i<topology_.numReaction() ; i++  )
      for( size_t j=0 ; j<topology_.reaction(i)->numParameter() ; j++  )
	parameter_.push_back( &( topology_.reaction(i)->parameterAddress(j) ) );
    //collect topology compartment update parameters DISABLED 4/25/21
    /*
    for( size_t i=0 ; i<topology_.numCompartmentChange() ; i++  )
      for( size_t j=0 ; j<topology_.compartmentChange(i)->numParameter() ; j++  )
	parameter_.push_back( &( topology_.compartmentChange(i)->parameterAddress(j) ) );
    */
  }  
  
  //collect species-reaction parameters
  for( size_t k=0 ; k<numSpecies() ; k++  )
    for( size_t i=0 ; i<species(k).numReaction() ; i++  )
      for( size_t j=0 ; j<species(k).reaction(i)->numParameter() ; j++  )
	parameter_.push_back( &( species(k).reaction(i)->parameterAddress(j) ));
  
  //collect reaction parameters
  for( size_t i=0 ; i<numReaction() ; i++  )
    for( size_t j=0 ; j<reaction(i)->numParameter() ; j++  )
      parameter_.push_back( &( reaction(i)->parameterAddress(j) ) );  
  
  //collect neighborhood parameters
  if( numNeighborhood() ) {
    for( size_t j=0 ; j<neighborhood().numParameter() ; j++  )
      parameter_.push_back( &( neighborhood().parameterAddress(j) ) );  
  }
}

void Organism::neighborhoodCreate(DataMatrix &y,double t) 
{ 
  if( numNeighborhood() )
    neighborhoodPointer()->create((*this),y,t); 
}

void Organism::neighborhoodUpdate(DataMatrix &y,double t ) 
{
  if( numNeighborhood() )
    neighborhoodPointer()->update((*this),y,t);
}

void Organism::removeCompartment(size_t i) 
{  
  int verbose=0;
  size_t N=numCompartment();
  size_t NN=N-1;
  if( verbose )
    std::cerr << "Organism::removeCompartment " << i << "\n";
  //Check index
  assert( i<=NN );
  
  //Remove neighbors pointing to removed compartment 
  for( size_t k=0 ; k<compartment(i).numNeighbor() ; k++ ) {
    int removed = compartment( compartment(i).neighbor(k) ).
      removeNeighbor(i);
    if( verbose ) {
      if( removed != 1 ) {
	std::cerr << "Organism::removeCompartment(int) "
		  << "Wrong when removing neighbor index, ( "
		  << removed << " index " << i << " removed from " 
		  << compartment(i).neighbor(k) << "," << numCompartment() 
		  << " )\n";
	for( size_t kk=0 ; kk<compartment(i).numNeighbor() ; kk++ )
	  std::cerr << compartment(i).neighbor(kk) << " ";
	std::cerr << "\n";
	for( size_t kk=0 ; kk<compartment(compartment(i).neighbor(k)).numNeighbor()
	       ; kk++ )
	  std::cerr << compartment(compartment(i).neighbor(k)).neighbor(kk) 
		    << " ";
	std::cerr << "\n";
	//exit(-1);
      }
    }
  }
  
  //(only if not last compartment)
  if( i != NN ) {
    //Change neighbor index for compartments pointing to last compartment
    for( size_t k=0 ; k<compartment(NN).numNeighbor() ; k++ ) {
      unsigned int moved = compartment( compartment(NN).neighbor(k) ).
	changeNeighborIndex(NN,i); 
      if( verbose ) {
	if( moved != 1 ) {
	  std::cerr << "Organism::removeCompartment(int) "
		    << "Wrong when moving neighbor index, ( "
		    << moved << " index " << i << " moved from " << NN << " in " 
		    << compartment(NN).neighbor(k) << "," << numCompartment()
		    << " )\n";
	  for( size_t kk=0 ; kk<compartment(NN).numNeighbor() ; kk++ )
	    std::cerr << compartment(NN).neighbor(kk) << " ";
	  std::cerr << "\n";
	  for( size_t kk=0 ; kk<compartment(compartment(NN).neighbor(k)).
		 numNeighbor()
		 ; kk++ )
	    std::cerr << compartment(compartment(NN).neighbor(k)).neighbor(kk) 
		      << " ";
	  std::cerr << "\n";
	  //exit(-1);
	}
      }
    }
    //Copy last compartment element into i and change index
    setCompartment(i,compartment(NN));
    compartment(i).setIndex(i);
  }
  //Truncate compartment vector
  compartment_.pop_back();  
  assert( numCompartment() == NN );
}

void Organism::derivs(DataMatrix &y,DataMatrix &dydt) 
{  
  //Check sizes for y an dydt
  assert( numCompartment() == y.size() );
  assert( numCompartment() == dydt.size() );
  
  //Check compartment indeces since these are used in the derivatives
  for( size_t i=0 ; i<numCompartment() ; i++ )
    assert( i==compartment(i).index() );
  
  //Check species indeces since these are used in the derivatives
  if( numCompartment() )
    for( size_t i=0 ; i<numSpecies() ; i++ )
      assert( i == species(i).index()-compartment(0).speciesStart() );
  
  // Zero all derivative variables before update
  //
  for( size_t i=0 ; i<dydt.size() ; ++i )
    for( size_t j=0 ; j<dydt[i].size() ; ++j )
      dydt[i][j] = 0.;
  //
  // Start update Topology variables
  //
  if( numTopologyVariable() ) {
    if( topology().numReaction() )
      for( size_t i=0 ; i<numCompartment() ; i++ )
	topology().derivs(compartment(i),y,dydt);
  }
  //
  // Update all species
  //
  if( numSpecies() ) {
    //update
    for( size_t s=0 ; s<numSpecies() ; s++ ) {
      if( species(s).numReaction() )
	for( size_t i=0 ; i<numCompartment() ; i++ )
	  species(s).derivs(compartment(i),y,dydt);  
    }
  }
  //
  // Update all additional reactions
  //
  if( numReaction() ) {
    for( size_t r=0 ; r<numReaction() ; r++ )
      for( size_t i=0 ; i<numCompartment() ; i++ )
	reaction(r)->derivs(compartment(i),static_cast<size_t>(-1),y,dydt);
  }
}



void Organism::derivsWithAbs(DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt) 
{  
  //Check sizes for y an dydt
  assert( numCompartment() == y.size() );
  assert( numCompartment() == dydt.size() );
  assert( numCompartment() == sdydt.size() );
  
  //Check compartment indeces since these are used in the derivatives
  for( size_t i=0 ; i<numCompartment() ; i++ )
    assert( i==compartment(i).index() );
  
  //Check species indeces since these are used in the derivatives
  if( numCompartment() )
    for( size_t i=0 ; i<numSpecies() ; i++ )
      assert( i == species(i).index()-compartment(0).speciesStart() );
  
  // Zero all derivative variables (and absolute values) before update
  //
  for( size_t i=0 ; i<dydt.size() ; ++i ) {
    for( size_t j=0 ; j<dydt[i].size() ; ++j ) {
      dydt[i][j] = 0.;
      sdydt[i][j] = 0.;
    }
  }
  //
  // Start update Topology variables
  //
  if( numTopologyVariable() ) {
    if( topology().numReaction() )
      for( size_t i=0 ; i<numCompartment() ; i++ )
	topology().derivsWithAbs(compartment(i),y,dydt,sdydt);
  }
  //
  // Update all species
  //
  if( numSpecies() ) {
    //update
    for( size_t s=0 ; s<numSpecies() ; s++ ) {
      if( species(s).numReaction() )
	for( size_t i=0 ; i<numCompartment() ; i++ )
	  species(s).derivsWithAbs(compartment(i),y,dydt,sdydt);  
    }
  }
  //
  // Update all additional reactions
  //
  if( numReaction() ) {
    for( size_t r=0 ; r<numReaction() ; r++ )
      for( size_t i=0 ; i<numCompartment() ; i++ )
	reaction(r)->derivsWithAbs(compartment(i),static_cast<size_t>(-1),y,dydt,sdydt);
  }
}

//DISABLED 042521
/*
size_t Organism::Jacobian(DataMatrix &y,JacobianMatrix &A)
{
  //size_t verbose=0;
  //Check sizes for y and A
  assert( A.size1() == y.size()*y[0].size() );
  size_t updateFlag=0;
  
  //if (updateFlag==0)
  // return updateFlag;
  //updateFlag = 0;
  //
  // Zero all Jacobian (A) elements before adding
  //
  A.clear();
  assert( A.size1() == y.size()*y[0].size() );
  //
  // Update topology variables
  //
  if( numTopologyVariable() ) {
    //update
    if( topology().numReaction() )
      for( size_t i=0 ; i<numCompartment() ; i++ )
	updateFlag += topology().Jacobian(compartment(i),y,A);
  }
  //
  // Update for all reactions for all species
  //
  if( numSpecies() ) {
    //update
    for( size_t s=0 ; s<numSpecies() ; s++ ) {
      if( species(s).numReaction() )
	for( size_t i=0 ; i<numCompartment() ; i++ )
	  updateFlag += species(s).Jacobian(compartment(i),y,A);  
    }
  }
  //
  // Update for all additional reactions
  //
  if( numReaction() ) {
    for( size_t r=0 ; r<numReaction() ; r++ )
      for( size_t i=0 ; i<numCompartment() ; i++ )
	reaction(r)->Jacobian(compartment(i),static_cast<size_t>(-1),y,A);
  }  
  //  std::cerr << "Organism::updateJacobian() Matrix after update:" << std::endl;
  //for (size_t i=0; i<A.size1(); ++i) {
  //for (size_t j=0; j<A.size2(); ++j)
  //  std::cerr << A(i,j) << " ";
  //std::cerr << std::endl;
  //}
  //std::cerr << std::endl;
  return updateFlag;
}
*/

void Organism::derivsMechanical(DataMatrix &y,DataMatrix &dydt,
				std::vector<size_t> &mechEq,
				std::vector<size_t> &posVar ) 
{  
  // Positional variables
  //////////////////////////////////////////////////////////////////////
  //set all derivative variables to zero
  for( size_t i=0 ; i<numCompartment() ; i++ )
    for( size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      dydt[i][j] = 0.;
    }
  //update
  for( size_t i=0 ; i<numCompartment() ; i++ )
    for( size_t k=0 ; k<mechEq.size() ; k++ ) {
      size_t j=mechEq[k];
      topology().reaction(j)->
	derivs(compartment(i),compartment(i).topologyStart(),y,dydt);
    }
  //topology().derivs(compartment(i),y,dydt);
}

void Organism::derivsTemplate(DataMatrix &y,DataMatrix &dydt,
			      std::vector<int> &simFlag ) 
{  
  //Check sizes for y an dydt
  if( numCompartment() != y.size() || numCompartment() != dydt.size() 
      || numVariable() != simFlag.size() ) {
    std::cerr << "Organism::derivsTemplate() Wrong size of y/dydt/simFlag.\n";
    std::cerr << y.size() << " " << dydt.size() << " " << numCompartment() 
	      << "\n" << simFlag.size() << " " << numVariable() << "\n";
    exit(-1);
  }
  
  //Check compartment indeces since these are used in the derivatives
  for( size_t i=0 ; i<numCompartment() ; i++ )
    if( i != compartment(i).index() ) {
      std::cerr << "Organism::derivs() Wrong index in compartment " << i 
		<< " (" << compartment(i).index() << ")\n";
      exit(-1);
    }
  //Check species indeces since these are used in the derivatives
  if( numCompartment() ) {
    for( size_t i=0 ; i<numSpecies() ; i++ )
      if( i != species(i).index()-compartment(0).speciesStart() ) {
	std::cerr << "Organism::derivs() Wrong index in species " << i << " ("
		  << species(i).index() << ")\n";
	exit(-1);
      }
  }
  //Check size of simFlag vector
  if( numCompartment() ) {
    if( compartment(0).numVariable() != simFlag.size() ) {
      std::cerr << "Organism::derivsTemplate() Wrong size of simFlag "
		<< "vector.\n";
      exit(-1);
    }
  }
  //Check that embedding and topology variables are flagged in simFlag
  size_t speciesStart = numTopologyVariable();
  //   if( numCompartment() ) {
  //     for( int i=0 ; i<speciesStart ; i++ )
  //       if( simFlag[i] ) { 
  // 	std::cerr << "Organism::derivsTemplate() Embedding or Topoloogy " 
  // 		  << "variable marked for simulation.\n";
  // 	exit(-1);
  //       }
  //   }
  //Update Species marked for simulation
  //////////////////////////////////////////////////////////////////////
  //all derivative variables to zero
  for( size_t j=speciesStart ; 
       j<compartment(0).numVariable() ; j++ )
    if( simFlag[j] )
      for( size_t i=0 ; i<numCompartment() ; i++ )
	dydt[i][j] = 0.;
  //update
  for( size_t s=0 ; s<numSpecies() ; s++ )
    if( simFlag[species(s).index()] )
      for( size_t i=0 ; i<numCompartment() ; i++ )
	species(s).derivs(compartment(i),y,dydt);  
  
  //All additional reactions
  //////////////////////////////////////////////////////////////////////
  for( size_t r=0 ; r<numReaction() ; r++ )
    for( size_t i=0 ; i<numCompartment() ; i++ )
      reaction(r)->derivs(compartment(i),static_cast<size_t>(-1),y,dydt);
}

void Organism::
reactionInitiate( double t,
		  std::vector< std::vector<double> > &y)
{
  // Initiate Topology reactions
  if( numTopologyVariable() && topology().numReaction() )
    for (size_t rCount=0; rCount<topology().numReaction(); rCount++ )
      topology().reaction(rCount)->initiate(t,y);
  // Initiate all species reactions
  if( numSpecies() ) {
    for( size_t s=0 ; s<numSpecies() ; s++ ) {
      if( species(s).numReaction() )
	for (size_t rCount=0; rCount<species(s).numReaction(); rCount++ )
	  species(s).reaction(rCount)->initiate(t,y);
    }
  }
  // Initiate additional (multispecies) reactions
  if( numReaction() ) {
    for( size_t r=0 ; r<numReaction() ; r++ )
      reaction(r)->initiate(t,y);
  }
}

void Organism::
reactionUpdate( double h, double t,
		std::vector< std::vector<double> > &y)
{
  // Update Topology variables
  if( numTopologyVariable() && topology().numReaction() )
    for (size_t rCount=0; rCount<topology().numReaction(); rCount++ )
      topology().reaction(rCount)->update(h,t,y);
  // Update all species reactions
  if( numSpecies() ) {
    for( size_t s=0 ; s<numSpecies() ; s++ ) {
      if( species(s).numReaction() )
	for (size_t rCount=0; rCount<species(s).numReaction(); rCount++ )
	  species(s).reaction(rCount)->update(h,t,y);
    }
  }
  // Update additional reactions
  if( numReaction() ) {
    for( size_t r=0 ; r<numReaction() ; r++ )
      reaction(r)->update(h,t,y);
  }
}

//DISABLED 042521
/*
void Organism::
compartmentChange(DataMatrix &y,DataMatrix &dydt,
		  double t) 
{  
  for( size_t k=0 ; k<topology().numCompartmentChange() ; k++ ) {
    for( size_t i=0 ; i<numCompartment() ; i++ ) {
      if( topology().compartmentChange(k)->flag(compartment(i),y,dydt) ) {
	topology().compartmentChange(k)->update((*this),compartment(i),
						y,dydt,t);
	if( topology().compartmentChange(k)->numChange()==-1 )
	  i--;
      }
    }
  }
}
*/

void Organism::
createReactionList(std::vector<BaseReaction*> &reactionList,
		   std::vector<size_t> &speciesList)
{
  // Get all species reactions
  if( numSpecies() ) {
    for( size_t s=0 ; s<numSpecies() ; s++ ) {
      if( species(s).numReaction() )
	for (size_t rCount=0; rCount<species(s).numReaction(); rCount++ ) {
	  reactionList.push_back( species(s).reaction(rCount) );
	  speciesList.push_back( species(s).index() );
	}
    }
  }
  // Update additional reactions
  size_t tmp = static_cast<size_t>(-1);
  if( numReaction() ) {
    for( size_t r=0 ; r<numReaction() ; r++ ) {
      reactionList.push_back( reaction(r) );
      speciesList.push_back( tmp );
    }	  
  }
}


void Organism::scaleSpace( double factor) 
{
  if( factor<=0.0 ) {
    std::cerr << "Organism::scaleSpace - Wrong scaling factor given!\n";
    exit(-1);
  }
  
  for( size_t i=0 ; i<numCompartment() ; i++ ) {
    unsigned int dim=numDimension();
    double factor2 = std::pow(factor,(int) (dim-1));
    double factor3 = std::pow(factor,(int) dim);
    
    //Scale positions
    for( unsigned int j=0 ; j<dim ; j++ )
      compartment(i).setVariable(j,factor*compartment(i).variable(j));
    //Scale the volume
    compartment(i).setVariable(dim,factor3*compartment(i).variable(dim));
    
    //Scale possible neighbor areas
    if( compartment(i).numNeighborArea() ) {
      for( unsigned int j=0 ; j<compartment(i).numNeighborArea() ; j++ )
	compartment(i).
	  setNeighborArea(j,factor2*compartment(i).neighborArea(j));
    }
  }
}

unsigned int Organism::
findPeaksGradientAscent(const DataMatrix &y,
			size_t col, std::vector<size_t> &cellMax,
			std::vector<size_t> &flag ) 
{  
  size_t N = numCompartment();
  //Check size of y
  assert( y.size() == N );
  assert( y[0].size()>col );
  if( cellMax.size() )
    cellMax.resize(0);
  if( flag.size() != N ) flag.resize(N);
  for( size_t i=0 ; i<N ; i++ )
    flag[i]=0;
  
  std::vector<size_t> cellTmp;//Values before threshold check
  std::vector<unsigned int> numTmp;//times cellTmp been visited
  std::vector<size_t> walkTmp;//positions for a walk (start point)
  
  size_t count=1;
  //Find the maxima from each pixel
  for( size_t iStart=0 ; iStart<numCompartment() ; iStart++ ) {
    size_t i=iStart;
    double value,newValue;
    walkTmp.resize(1);
    walkTmp[0]=i;
    //find the max by walking uphill (greedy)
    if( !flag[i] ) {
      do {
	newValue=value=y[i][col];
	size_t newI=i;
	//Check all neighboring cells
	for(size_t k=0 ; k<compartment(i).numNeighbor() ; k++ ) {
	  size_t j = compartment(i).neighbor(k);
	  if( y[j][col]>newValue ) {
	    newValue=y[j][col];
	    newI=j;
	  }
	}
	i=newI;
	walkTmp.push_back( i );
      } while( newValue>value && !flag[i] );
    }
    //Collect the path data and add one visit for the maximum
    if( !flag[i] ) { //new maximum
      cellTmp.push_back( i );
      numTmp.push_back(1);
      unsigned int n=count++;//cellTmp.size();
      for( size_t a=0 ; a<walkTmp.size() ; a++ )
	flag[ walkTmp[a] ] = n;
    }
    else { //old maximum or background
      size_t n = flag[i];
      for( size_t a=0 ; a<walkTmp.size() ; a++ )
	flag[ walkTmp[a] ] = n;
      if( flag[i]>0 )//old maxima
	numTmp[n-1]++;
    }
  }
  //No threshold checking...
  unsigned int threshold = 1;
  double valThreshold = 0.0;
  //Get the maxima visited more than threshold times and with an intensity
  //value higher than threshold
  std::vector<int> clusterNum;
  for( size_t n=0 ; n<cellTmp.size() ; n++ )
    if( numTmp[n]>=threshold && y[ cellTmp[n] ][col]>valThreshold ) {
      cellMax.push_back( cellTmp[n] ); 
      clusterNum.push_back( n+1 );
    }
  
  //Save the basins of attraction
  //   boa.resize( cellMax.size() );
  //   for( int i=0 ; i<numCompartment() ; i++ )
  //     for( int n=0 ; n<cellMax.size() ; n++ )
  //       if( flag[i] == clusterNum[n] )
  // 	boa[n].push_back(i);
  
  return static_cast<unsigned int>(cellMax.size());
}


void Organism::printModel(std::ostream &os) const 
{



  os << id() << " " << numTopology() << " " << numSpecies() << " " 
     << numReaction() << " " << numNeighborhood() << "\n\n";
  //
  // Topology
  //
  if( numTopology() ) {
    os << topology_.id() << " " << numTopologyVariable() << " " 
       << topology_.numReaction() << " " 
       << "DISABLED 042521 topology_.numCompartmentChange()" << "\n";
    // reactions
    for( size_t i=0 ; i<topology_.numReaction() ; i++ ) {
      os << topology_.reaction(i)->id() << " " 
				 << topology_.reaction(i)->numParameter() << " " 
				 << topology_.reaction(i)->numVariableIndexLevel() << " ";
      for( size_t k=0 ; k<topology_.reaction(i)->numVariableIndexLevel() ; k++ )
				os << topology_.reaction(i)->numVariableIndex(k) << " ";
      os << "\n";
      for( size_t j=0 ; j<topology_.reaction(i)->numParameter() ; j++ )
				os << topology_.reaction(i)->parameter(j) << " ";
      os << "\n";
      for( size_t k=0 ; k<topology_.reaction(i)->numVariableIndexLevel() ; k++ ) {
				for( size_t j=0 ; j<topology_.reaction(i)->numVariableIndex(k) ; j++ )
					os << topology_.reaction(i)->variableIndex(k,j) << " ";
				os << "\n";
      }
    }
    //DISABLED
    /*
    // compartmentChanges
    for( size_t i=0 ; i<topology_.numCompartmentChange() ; i++ ) {
      os << topology_.compartmentChange(i)->id() << " " 
				 << topology_.compartmentChange(i)->numParameter() << " " 
				 << topology_.compartmentChange(i)->numVariableIndexLevel() << " ";
      for( size_t k=0 ; k<topology_.compartmentChange(i)->numVariableIndexLevel() ; k++ )
				os << topology_.compartmentChange(i)->numVariableIndex(k) << " ";
      os << "\n";
      for( size_t j=0 ; j<topology_.compartmentChange(i)->numParameter() ; j++ )
				os << topology_.compartmentChange(i)->parameter(j) << " ";
      os << "\n";
      for( size_t k=0 ; k<topology_.compartmentChange(i)->numVariableIndexLevel() ; k++ ) {
				for( size_t j=0 ; j<topology_.compartmentChange(i)->numVariableIndex(k) ; j++ )
					os << topology_.compartmentChange(i)->variableIndex(k,j) << " ";
				os << "\n";
      }
    }
    */
    os << "\n";
  }
  

  
  //
  // All species
  //
  for( size_t s=0 ; s<numSpecies() ; s++ ) {
    os << species(s).id() << " " << species(s).index() << " "
       << species(s).numReaction() << "\n";
    // reactions
    for( size_t i=0 ; i<species(s).numReaction() ; i++ ) {
      os << species(s).reaction(i)->id() << " " 
	 << species(s).reaction(i)->numParameter() << " " 
	 << species(s).reaction(i)->numVariableIndexLevel() << " ";
      for( size_t k=0 ; k<species(s).reaction(i)->numVariableIndexLevel() ; k++ )
	os << species(s).reaction(i)->numVariableIndex(k) << " ";
      os << "\n";      
      for( size_t j=0 ; j<species(s).reaction(i)->numParameter() ; j++ )
	os << species(s).reaction(i)->parameter(j) << " ";
      os << "\n";
      for( size_t k=0 ; k<species(s).reaction(i)->numVariableIndexLevel() ; k++ ) {
				for( size_t j=0 ; j<species(s).reaction(i)->numVariableIndex(k) ; j++ )
					os << species(s).reaction(i)->variableIndex(k,j) << " ";
				os << "\n";
      }
    }
    os << "\n";
  }
  
 
  
  //
  // Additional reactions
  //
  for( size_t i=0 ; i<numReaction() ; i++ ) {
    os << reaction(i)->id() << " " << reaction(i)->numParameter() << " " 
       << reaction(i)->numVariableIndexLevel() << " ";
    for( size_t k=0 ; k<reaction(i)->numVariableIndexLevel() ; k++ )
      os << reaction(i)->numVariableIndex(k) << " ";
    os << "\n";      
    for( size_t j=0 ; j<reaction(i)->numParameter() ; j++ )
      os << reaction(i)->parameter(j) << " ";
    os << "\n";
    for( size_t k=0 ; k<reaction(i)->numVariableIndexLevel() ; k++ ) {
      for( size_t j=0 ; j<reaction(i)->numVariableIndex(k) ; j++ )
	os << reaction(i)->variableIndex(k,j) << " ";
      os << "\n";
    }
  }
  

  //
  // Neighborhood
  //
  if( numNeighborhood() ) {
    os << neighborhood().id() << " " 
       << neighborhood().numParameter()
       << "\n";
    for( size_t i=0 ; i<neighborhood().numParameter() ; i++ )
      os << neighborhood().parameter(i) << " ";
    os << "\n";
    for( size_t k=0 ; k<neighborhood().numVariableIndexLevel() ; k++ ) {
      for( size_t j=0 ; j<neighborhood().numVariableIndex(k) ; j++ )
				os << neighborhood().variableIndex(k,j) << " ";
      os << "\n";      
    }
    //os << neighborhood().inFile() << "\n";
  }
} 

void Organism::printModelCambium( std::ostream &os) const 
{
  os << "(*Cambium model description file generated via organism solver*)" << std::endl
   << "(*Note: this is a test version.*)" << std::endl
   << "(*Contact: emj@uci.edu (Cambium), henrik.jonsson@slcu.cam.ac.uk*)" << std::endl
   << "(* http://emj.ics.uci.edu/software/ http://dev.thep.lu.se/organism *)" <<std::endl
   << "(* *)" << std::endl
   << "(* This is a translation of model: " << id() << "*)" << std::endl;
  
  //
  // Topology
  //
  if( numTopology() ) {
    os << "(* Topology reactions [" << topology_.id() << "] *)" << std::endl;
    // reactions
    for( size_t i=0 ; i<topology_.numReaction() ; i++ ) {
      topology_.reaction(i)->printCambium(os);
    }
    os << std::endl;
    //DISABLED 042521
    /*
    // compartmentChanges
    os << "(* Topology compartment changes *)" <<std::endl;
    for( size_t i=0 ; i<topology_.numCompartmentChange() ; i++ ) {
      topology_.compartmentChange(i)->printCambium(os);
    }
    os << std::endl;
    */
  }
  //
  // All species
  //
  os << "(* Molecular species reactions *)" <<std::endl;
  for( size_t s=0 ; s<numSpecies() ; s++ ) {
    os << "(* " << species(s).id() << " *)" << std::endl;
    // reactions
    for( size_t i=0 ; i<species(s).numReaction() ; i++ ) {
      species(s).reaction(i)->printCambium(os,s+numTopologyVariable());
    }
  }
  os << std::endl;
  //
  // Additional reactions
  //
  if (numReaction() ) {
    os << "(* Multi-species reactions *)" << std::endl;
    for( size_t i=0 ; i<numReaction() ; i++ ) {
      reaction(i)->printCambium(os);
    }
    os << std::endl;
  }
  //
  // Neighborhood
  //
  if( numNeighborhood() ) {
    os << "(* Neighborhood reactions *)" << std::endl;
    neighborhood().printCambium(os);  
  }
}

void Organism::printParameter(std::ostream &os) const 
{
  for( size_t i=0 ; i<numParameter() ; i++ )
    os << parameter(i) << " ";
  os << std::endl;
}

void Organism::printVariable(std::ostream &os) const 
{
  os << numCompartment() << " " << numVariable() << "\n";
  
  for( size_t i=0 ; i<numCompartment() ; i++ ) {
    for( size_t j=0 ; j<numVariable() ; j++ )
      os << compartment(i).variable(j) << " ";
    os << "\n";
  }
}

void Organism::printNeighbor(std::ostream &os) const 
{  
  for( size_t i=0 ; i<numCompartment() ; i++ ) {
    os << compartment(i).index() << " ";
    for( size_t k=0 ; k<compartment(i).numNeighbor() ; k++ ) {
      os << compartment(i).neighbor(k) << " ";
      if( compartment(i).numNeighborArea()>k )
	os << "(" << compartment(i).neighborArea(k) << ") ";
    }
    os << "\n";
  }
}
