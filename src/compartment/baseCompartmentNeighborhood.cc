//
// Filename     : baseCompartmentNeighborhood.cc
// Description  : A base class describing variable updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2003
// Revision     : $Id: baseCompartmentNeighborhood.cc 646 2015-08-23 18:35:06Z henrik $
//

#include<vector>
#include<boost/algorithm/string.hpp>


#include"baseCompartmentNeighborhood.h"

#include"compartmentNeighborhood.h"

BaseCompartmentNeighborhood::~BaseCompartmentNeighborhood(){}

BaseCompartmentNeighborhood *
BaseCompartmentNeighborhood::
createCompartmentNeighborhood(std::vector<double> &paraValue,
															std::vector< std::vector<size_t> > &indValue, 
															const std::string &idValue ) {
  
  // Distance rules
  //(compartmentNeighborhood.h,compartmentNeighborhood.cc)
  if( boost::iequals(idValue,"neighborhoodLatticeCigar" )) 
    return new NeighborhoodLatticeCigar(paraValue,indValue);
  else if( boost::iequals(idValue,"neighborhoodLatticeCigarMultipleWall" )) 
    return new NeighborhoodLatticeCigarMultipleWall(paraValue,indValue);

  else if( boost::iequals(idValue,"neighborhoodDistanceSphere") ||
	   boost::iequals(idValue,"Sphere::neighborhoodDistance") ) 
    return new Sphere::NeighborhoodDistance(paraValue,indValue);
  else if( boost::iequals(idValue,"neighborhoodDistanceSpherePeriodic" )) 
    return new NeighborhoodDistanceSpherePeriodic(paraValue,indValue);
  else if( boost::iequals(idValue,"neighborhoodDistanceSphereWithWalls" )) 
    return new NeighborhoodDistanceSphereWithWalls(paraValue,indValue);
  else if( boost::iequals(idValue,"neighborhoodDistanceSphereBud" ))
    return new NeighborhoodDistanceSphereBud(paraValue,indValue);
  else if( boost::iequals(idValue,"neighborhoodDistanceEllipse" ))
    return new NeighborhoodDistanceEllipse(paraValue,indValue);
  else if( boost::iequals(idValue,"neighborhoodDistanceCigar" ))
    return new NeighborhoodDistanceCigar(paraValue,indValue);
  else if ( boost::iequals(idValue,"neighborhoodDistanceSpherePeriodicBox"))
    return new NeighborhoodDistanceSpherePeriodicBox(paraValue, indValue);

  // From indices
  // (compartmentNeighborhood.h, compartmentNeighborhood.cc)
  else if (boost::iequals(idValue,"neighborhoodIndex"))
    return new NeighborhoodIndex(paraValue, indValue);
  // Null
  else if (boost::iequals(idValue,"nullNeighborhood"))
    return new NullNeighborhood(paraValue, indValue);

	
  // From file
  //(compartmentNeighborhood.h,compartmentNeighborhood.cc)
  else if( boost::iequals(idValue,"neighborhoodFromFileInitial" )) {
    std::cerr << "\nBaseCompartmentNeighborhood::"
	      << "createCompartmentNeighborhood(pVal,indVal,idVal) "
	      << "Should never find idValue=" << idValue << "\n";
    exit(EXIT_FAILURE);
  }
  //Default, if nothing found
  else {
    std::cerr << "\nBaseCompartmentNeighborhood::"
	      << "createCompartmentNeighborhood() "
	      << "WARNING: CompartmentNeighborhoodtype " << idValue 
	      << " not known, no compartmentNeighborhood created.\n\7";
    exit(EXIT_FAILURE);
  }
}

BaseCompartmentNeighborhood * 
BaseCompartmentNeighborhood::
createCompartmentNeighborhood(std::istream &IN ) {
  
 
  std::string idVal;
  size_t pNum,levelNum;
  IN >> idVal;
  
  if( idVal == "neighborhoodFromFileInitial" ) {
    IN >> pNum;
    std::vector<double> pVal( pNum );
    for( size_t i=0 ; i<pNum ; i++ )
      IN >> pVal[i];
    std::string fileName;
    IN >> fileName;
    return new NeighborhoodFromFileInitial(pVal,fileName);
  }
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
	
  return createCompartmentNeighborhood(pVal,varIndexVal,idVal);
}

unsigned int BaseCompartmentNeighborhood::
create(Organism &O,std::vector< std::vector<double> > &y,double t ) {
  std::cerr << "BaseCompartmentNeighborhood::create() should not be used. "
	    << "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

unsigned int BaseCompartmentNeighborhood::
update(Organism &O,std::vector< std::vector<double> > &y,double t ) {
  std::cerr << "BaseCompartmentNeighborhood::update() should not be used. "
	    << "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

void BaseCompartmentNeighborhood::
printCambium( std::ostream &os ) const
{
  std::cerr << "BaseCompartmentNeighborhood::printCambium(ofstream) should not be used. " << std::endl
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
