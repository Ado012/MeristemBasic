/**
 * Filename     : compartment.cc
 * Description  : A class describing a compartment with neighborhood
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : October 2003
 * Revision     : $Id: compartment.cc 637 2015-07-07 12:55:36Z henrik $
 */

#include"compartment.h"
#include"topology.h"
#include"organism.h"
#include <boost/algorithm/string.hpp>
#include <math.h>

Compartment::Compartment() {

  setIndex(static_cast<size_t>(-1));
  setVolume( 1.0 );
}


Compartment::Compartment( const Compartment & compartmentCopy ) {
  
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

Compartment::Compartment( const Organism* organismVal,
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

Compartment::Compartment( const Organism* organismVal,
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

Compartment::~Compartment() 
{
  //std::cerr << "Compartment::~Compartment()\n";
}

size_t Compartment::numDimension() const 
{
  return topology_->numDimension();
}

size_t Compartment::numTopologyVariable() const 
{
  return topology_->numVariable();
}

size_t Compartment::neighborIndex(size_t compartmentIndex) const
{
	for (size_t k=0; k<numNeighbor(); ++k) {
		if (neighbor(k)==compartmentIndex)
			return k;
	}
	return static_cast<size_t>(-1);
}

const Compartment & Compartment::neighborRef(size_t i)
{
	return organism_->compartment(neighbor(i));
}

unsigned int Compartment::removeNeighbor(size_t index) 
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

unsigned int Compartment::changeNeighborIndex(size_t oldIndex,
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

unsigned int Compartment::removeNeighborAtLevel(size_t l,size_t index) 
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

unsigned int Compartment::
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

void Compartment::setTopology(const Topology* value) 
{
  topology_ = value;
}

void Compartment::setOrganism(const Organism* value) 
{
  organism_ = value;
}

double Compartment::getVolume()
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
      std::cerr << "Compartment::getVolume() wrong dimension!" << std::endl;
      exit(EXIT_FAILURE);      
    }
  }
  else {	
    return volume_ = r;
  }
}

double Compartment::getVolume(std::vector<double> &y)
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
      std::cerr << "Compartment::getVolume() wrong dimension!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else {	
    return volume_ = r;
  }
}
