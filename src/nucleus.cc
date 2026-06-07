/**
 * Filename     : nucleus.cc
 * Description  : A class describing a compartment within a Compartment
 * Author(s)    : Albert Do (ado012@ucr.edu)
 * Created      : March 2017
 * Revision     : v0.1
 */

#include"nucleus.h"
#include"topology.h"
#include"organism.h"
#include <boost/algorithm/string.hpp>
#include <math.h>

Nucleus::Nucleus() {

  setIndex(static_cast<size_t>(-1));
  setVolume( 1.0 );
}


Nucleus::Nucleus( const Nucleus & compartmentCopy ) {
  
  setIndex( compartmentCopy.index() );
  setId( compartmentCopy.id() );
  setTopology( compartmentCopy.topology() );
  setOrganism( compartmentCopy.organism() );
  setVariable( compartmentCopy.variable() );
  setNeighbor( compartmentCopy.neighbor() );
  setNeighborArea( compartmentCopy.neighborArea() );
  setVolume( compartmentCopy.volume() );
  setParameter( compartmentCopy.parameter() );
}

Nucleus::Nucleus( const Organism* organismVal,
													const Topology* topologyVal, 
													size_t numSpeciesVal,size_t indexVal, 
													const std::string &idVal) 
{
  setIndex( indexVal );
  setId( idVal );
  setTopology( topologyVal );
  setOrganism( organismVal );
  
  setVolume( 1.0 );
  std::vector<double> tmpVarVal( numTopologyVariable()+numSpeciesVal );
  setVariable( tmpVarVal );
}

Nucleus::Nucleus( const Organism* organismVal,
													const Topology* topologyVal,
													const std::vector<double> &variableVal,
													size_t indexVal, 
													const std::string &idVal ) 
{
  setIndex( indexVal );
  setId( idVal );
  setTopology( topologyVal );
  setOrganism( organismVal );
  setVariable( variableVal );
  
  std::vector<double> tmpD;
  std::vector<size_t> tmpI;
  setParameter(tmpD);
  setNeighbor(tmpI);
  setNeighborArea(tmpD);
  
}

Nucleus::~Nucleus() 
{
  //std::cerr << "Nucleus::~Nucleus()\n";
}

size_t Nucleus::numDimension() const 
{
  return topology_->numDimension();
}

size_t Nucleus::numTopologyVariable() const 
{
  return topology_->numVariable();
}

size_t Nucleus::neighborIndex(size_t compartmentIndex) const
{
	for (size_t k=0; k<numNeighbor(); ++k) {
		if (neighbor(k)==compartmentIndex)
			return k;
	}
	return static_cast<size_t>(-1);
}


//Nucleus has only one neighbor: To be reworked
/*
const Nucleus & Nucleus::neighborRef(size_t i)
{
	return organism_->compartment(neighbor(i));
}

unsigned int Nucleus::removeNeighbor(size_t index) 
{  
  if( numNeighborLevel()<1 ) return 0;
  unsigned int removed=0;
  std::vector<size_t> newNeigh;
  for( size_t n=0 ; n<numNeighbor() ; n++ ) {
    if( neighbor(n)==index )
      removed++;
    else
      newNeigh.push_back( neighbor(n) );
  }
  setNeighbor(newNeigh);
  return removed;
}

unsigned int Nucleus::changeNeighborIndex(size_t oldIndex,
																							size_t newIndex) 
{  
  if( numNeighborLevel()<1 ) return 0;
  unsigned int moved=0;
  for( size_t n=0 ; n<numNeighbor() ; n++ ) {
    if( neighbor(n)==oldIndex ) {
      setNeighbor(n,newIndex);
      moved++;
    }
  }
  return moved;
}

unsigned int Nucleus::removeNeighborAtLevel(size_t l,size_t index) 
{  
  if( numNeighborLevel()<l+1 ) return 0;
  unsigned int removed=0;
  std::vector<size_t> newNeigh;
  for( size_t n=0 ; n<numNeighborAtLevel(l) ; n++ ) {
    if( neighborAtLevel(l,n)==index )
      removed++;
    else
      newNeigh.push_back( neighborAtLevel(l,n) );
  }
  setNeighborAtLevel(l,newNeigh);
  return removed;
}

unsigned int Nucleus::
changeNeighborIndexAtLevel(size_t l,size_t oldIndex,size_t newIndex) 
{
  if( numNeighborLevel()<l+1 ) return 0;
  unsigned int moved=0;
  for( size_t n=0 ; n<numNeighborAtLevel(l) ; n++ ) {
    if( neighborAtLevel(l,n)==oldIndex ) {
      setNeighborAtLevel(l,n,newIndex);
      moved++;
    }
  }
  return moved;
}

*/

void Nucleus::setTopology(const Topology* value) 
{
  topology_ = value;
}

void Nucleus::setOrganism(const Organism* value) 
{
  organism_ = value;
}

double Nucleus::getVolume()
{
  size_t volumeIndex = organism_->numTopologyVariable()-1; //also dimension in sphere
  double r = variable_[volumeIndex];
  if (boost::iequals(topology_->id(),"sphere")) {
    if (volumeIndex==1) {
      return volume_=r;
    }
    else if (volumeIndex==2) {
      return volume_ = M_PI*r*r;
    }
    else if (volumeIndex==2) {
      return volume_ = 1.3333333*M_PI*r*r*r;
    }
    else {
      std::cerr << "Nucleus::getVolume() wrong dimension!" << std::endl;
      exit(EXIT_FAILURE);      
    }
  }
  else {	
    return volume_ = r;
  }
}

double Nucleus::getVolume(std::vector<double> &y)
{
  size_t volumeIndex = organism_->numTopologyVariable()-1; //also dimension in sphere
  double r = y[volumeIndex];
  if (boost::iequals(topology_->id(),"sphere")) {
    if (volumeIndex==1) {
      return volume_= r;
    }
    else if (volumeIndex==2) {
      return volume_ = M_PI*r*r;
    }
    else if (volumeIndex==2) {
      return volume_ = 1.3333333*M_PI*r*r*r;
    }
    else {
      std::cerr << "Nucleus::getVolume() wrong dimension!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else {	
    return volume_ = r;
  }
}
