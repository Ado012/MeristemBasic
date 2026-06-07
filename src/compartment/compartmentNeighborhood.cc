//
// Filename     : compartmentNeighborhood.cc
// Description  : Classes describing compartment neighborhoods
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2005
// Revision     : $Id: compartmentNeighborhood.cc 638 2015-07-08 07:30:02Z henrik $
//

#include<cassert>
#include<cmath>
#include <set>

#include"compartmentNeighborhood.h"
#include"baseCompartmentNeighborhood.h"
#include"../common/myRandom.h"

NeighborhoodLatticeCigarMultipleWall::
NeighborhoodLatticeCigarMultipleWall(std::vector<double> &paraValue, 
																		 std::vector< std::vector<size_t> > 
																		 &indValue ) 
{
  // Do some checks on the parameters and variable indeces
  if( paraValue.size() < 10 ) {
    std::cerr << "NeighborhoodLatticeCigarMultipleWall::"
							<< "NeighborhoodLatticeCigarMultipleWall() "
							<< "Must use at least one wall!\n"
							<< "five parameters numCell, cellSize, updateFlag, "
							<< "dt_min, b and b_frac followed by wall coordinates.\n";
    exit(0);
  }
  if( indValue.size() ) {
    std::cerr << "NeighborhoodLatticeCigarMultipleWall::"
							<< "NeighborhoodLatticeCigarMultipleWall() "
							<< "No variable index is used.\n";
    exit(0);
  }

  // Check for consistency in parameter values
	if( paraValue[2]!=0.0 && paraValue[2]!=1.0 ) {
    std::cerr << "NeighborhoodLatticeCigarMultipleWall::"
							<< "NeighborhoodLatticeCigarMultipleWall() "
							<< "parameter(1) (updateFlag) must be 0 or 1.\n";
    exit(0);
  }

	// Set the variable values
  setId("neighborhoodLatticeCigarMultipleWall");
  setParameter(paraValue);  
  setVariableIndex(indValue);

	//Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "numCell";
  tmp[1] = "cellSize";
  tmp[2] = "updateFlag";
  tmp[3] = "dt_min";
  tmp[4] = "b";
  tmp[5] = "b_frac";
  setParameterId( tmp );
}

unsigned int NeighborhoodLatticeCigarMultipleWall::
create(Organism &O, std::vector< std::vector<double> > &y, double t ) 
{   
  if( O.numTopologyVariable() != 5 && 
      O.numTopologyVariable() != 9 ) {
    std::cerr << "NeighborhoodLatticeCigar::create() "
							<< "Wrong number of topology variables.\n";
    exit(-1);
  }
  size_t N = O.numCompartment();
  size_t numWall = static_cast<size_t>( (numParameter() - 6)/4 );  
  std::vector< std::vector<size_t> > neighbor(N), neighborWall(N);
  unsigned int numNeigh=0;
  double cellSize = parameter(1);
  double cellSizeInverse = 1/cellSize;
  int  numCell = static_cast<int>(parameter(0));
  size_t dimension=  static_cast<size_t>( (O.numTopologyVariable()-1)/2 );
  //only in the case where we havv x''=F 
  if(O.numTopologyVariable() == 9 ) 
		dimension = 2;
  size_t totalNumCells = static_cast<size_t>( std::pow((double) numCell, (int) dimension) );
  
  std::vector<double> yc(dimension);
  std::vector< size_t > ic1( N );//given a compartment index i, gives lattice cell index
  std::vector< std::vector< size_t > > latticeCell(totalNumCells);
  double distanceThreshold = 0.5*(std::sqrt(2)+1.0)*cellSize + parameter(4);
	
  //This is done only in create(), not in update() 
  wallNeighborLattice.resize(numWall);
  for(int ly = 0; ly < numCell; ++ly) {
    for(int lx = 0; lx < numCell; ++lx) {
      for(size_t w = 0; w < numWall; ++w) {
				//Wall properties
				std::vector<double> n(dimension),  x1(dimension),x2(dimension);
				x1[0] = parameter(6+w*4);
				x1[1] = parameter(7+w*4);
				x2[0] = parameter(8+w*4);
				x2[1] = parameter(9+w*4);
				double length = 0;
				for(size_t dim = 0; dim < dimension; ++dim) {
					n[dim] = x2[dim] - x1[dim];
					length += ( x2[dim] -x1[dim])*( x2[dim] -x1[dim]);
				}
				length = std::sqrt(length);
				for(size_t dim = 0; dim < dimension; ++dim) 
					n[dim] = n[dim]/length;
				//lattice cell properties
				size_t latticeIndex = lx + ly*numCell;
				std::vector<double> cellCOM(dimension);
				cellCOM[0] = cellSize*( lx % numCell + 0.5);
				cellCOM[1] = cellSize*( ly % numCell + 0.5);
				//compute closest distance between cell and wall if it's close enough add as neighbor
				double dist = (x1[1] - cellCOM[1])*n[0] - 
					(x1[0] - cellCOM[0])*n[1];
				double s = n[0]<n[1] ? 
					(cellCOM[0]-x1[0]-n[1]*dist)/n[0] : 
					(cellCOM[1]-x1[1]+n[0]*dist)/n[1] ;
				if( s<0.0 || s>length ) {
					double dist1= std::sqrt( (x1[0]-cellCOM[0])*(x1[0]-cellCOM[0])+
																	 (x1[1]-cellCOM[1])*(x1[1]-cellCOM[1]) );
					double dist2= std::sqrt( (x2[0]-cellCOM[0])*(x2[0]-cellCOM[0])+
																	 (x2[1]-cellCOM[1])*(x2[1]-cellCOM[1]) );
					dist = dist1<dist2 ? dist1 : dist2;
				}
				if( std::fabs(dist) < distanceThreshold )
					wallNeighborLattice[w].push_back(latticeIndex);
      }
      
    }
  }
	for (size_t i =0; i<wallNeighborLattice.size(); ++i) {
 		std::cerr << "Wall " << i  << " is in boxes:\n\t";
 		for (size_t j =0; j<wallNeighborLattice[i].size(); ++j)
 			std::cerr << wallNeighborLattice[i][j] << " ";
 		std::cerr << "\n";
 	}
  //set up lattice, and place each cell point in a lattice cell
  for( size_t d = 0; d<dimension; d++){     
    int powNumc =  int( std::pow((double) numCell,(int) d));
    for(size_t i = 0; i<N;i++){
      yc[d] = 0.5*(y[i][ d ]+y[i][ d+dimension ]);
      if( (yc[d] < 0 || yc[d] > cellSize*numCell) ){
				std::cerr << "NeighborhoodLatticeCigar::create(): "
									<< "compartment is outside lattice\n";
				exit(-1);
      }
      ic1[i] +=  powNumc*static_cast<size_t>( yc[d]*cellSizeInverse );
    }
  }
  
  for(size_t i = 0; i<N;i++){
    latticeCell[ ic1[i] ].push_back(i);
  }
  //for each cell compartment, add the cells in neighboring lattice cells as neighbors
  for( size_t i=0 ; i<N ; i++ ){
    std::set<size_t> cellIndices;
    std::set<size_t>::iterator cellIndex;
    //std::cerr << "ic1 " << ic1[i] << ": ";
    bool flag = 1;
    std::vector<int> k(dimension);
    if( dimension == 2) {
      for(k[1] = -numCell; k[1] <= numCell; k[1]+=numCell) {
				for(k[0] = -1; k[0] <= 1; k[0]++){
					cellIndices.insert(ic1[i]+ k[0] + k[1]);
				}
      }
    } 
    else if( dimension == 3){
      for(k[2] = - numCell*numCell; k[2] <= numCell*numCell; k[2] += numCell*numCell) {
				for(k[1] = -numCell; k[1] <= numCell; k[1]+=numCell) {
					for(k[0] = -1; k[0] <= 1; k[0]++){
						cellIndices.insert(ic1[i]+ k[0] + k[1] + k[2]);
					}
				}
      }
    }
    else {
      std::cerr << "NeighborhoodLatticeCigar::create() "
								<< "Wrong number of dimensions.\n";
      exit(-1);
    }
    for(cellIndex = cellIndices.begin(); cellIndex != cellIndices.end(); cellIndex++) {
      if( (*cellIndex >= 0 && *cellIndex <= totalNumCells) 
					&& latticeCell[*cellIndex].size()) { 
				for(size_t n = 0; n < latticeCell[*cellIndex].size(); n++){
					std::cerr << *cellIndex << "\n";
					flag = 1;
					size_t LCell = latticeCell[*cellIndex][n];
					if(neighbor[i].size()){
						for(size_t l = 0; l<neighbor[i].size(); l++) {
							if(LCell == neighbor[i][l])
								flag = 0;
						}
					}
					if(LCell > i && flag) {
						//bool control = check(i,LCell, dimension, y);
						//bool control = check2(y[i],y[LCell]);
						bool control = 1;
						if( control){
							neighbor[i].push_back(LCell );
							neighbor[LCell].push_back(i);
							numNeigh++;
						}
					}
				}
      }
    }
  }
  
  //Adds walls into cigar neighbors
  for(size_t w = 0; w < numWall; ++w) {
    for(size_t m = 0; m < wallNeighborLattice[w].size(); ++m) {
      for(size_t n = 0 ; 
					n < latticeCell[ wallNeighborLattice[w][m] ].size() ; ++n) {
				neighborWall[latticeCell[ wallNeighborLattice[w][m] ][n] ].push_back(w); 
      }
    }
  }  
  //Copy the new neighbours into the cells (Caveat: Not a good solution)
  for(size_t i=0 ; i<N ; i++ ) {
    O.compartment(i).setNeighbor( neighbor[i] );
    O.compartment(i).setNeighborAtLevel( 1, neighborWall[i] );
  }
  
  //Set the previos time variable to current t value
  setPreviousTime(t);
  setNumNeighbor(numNeigh);
  
  return numNeigh;
}

unsigned int NeighborhoodLatticeCigarMultipleWall::
update(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{
  if( parameter(2) == 1.0 && t-previousTime()>parameter(3) ) {
    
    if( O.numTopologyVariable() != 5 && 
	O.numTopologyVariable() != 9 ) {
      std::cerr << "NeighborhoodLatticeCigar::create() "
		<< "Wrong number of topology variables.\n";
      exit(-1);
    }
    size_t N = O.numCompartment();
    size_t numWall = static_cast<size_t>( (numParameter() - 6)/4 );  
    std::vector< std::vector<size_t> > neighbor(N), neighborWall(N);
    unsigned int numNeigh=0;
    double cellSize = parameter(1);
    double cellSizeInverse = 1/cellSize;
    int  numCell = static_cast<int>(parameter(0));
    size_t dimension=  static_cast<size_t>( (O.numTopologyVariable()-1)/2 );
    //only in the case where we havv x''=F 
    if(O.numTopologyVariable() == 9 ) 
      dimension = 2;
    size_t totalNumCells = static_cast<size_t>( std::pow((double) numCell, (int) dimension) );
    
    std::vector<double> yc(dimension);
    std::vector< size_t > ic1( N );//given a compartment index i, gives lattice cell index
    std::vector< std::vector< size_t > > latticeCell(totalNumCells);
    
    //set up lattice, and place each cell point in a lattice cell
    for( size_t d = 0; d<dimension; d++){     
      int powNumc =  int( std::pow((double) numCell, (int) d) );
      for(size_t i = 0; i<N;i++){
	yc[d] = 0.5*(y[i][ d ]+y[i][ d+dimension ]);
	if( (yc[d] < cellSize || yc[d] > cellSize*(numCell-1)) ){
	  std::cerr << "NeighborhoodLatticeCigarMultipleWall::create(): "
		    << "compartment is outside lattice\n";
	  std::cerr << y[i][0] << " " << y[i][1] << "\n"
		    << y[i][2] << " " << y[i][3] << "\n";
	  exit(-1);
	}
	ic1[i] +=  powNumc*static_cast<size_t>( yc[d]*cellSizeInverse );
      }
    }
    
    for(size_t i = 0; i<N;i++){
      latticeCell[ ic1[i] ].push_back(i);
    }
    
    //for each cell compartment, add the cells in neighboring lattice cells as neighbors
    for( size_t i=0 ; i<N ; i++ ){
      std::set<size_t> cellIndices;
      std::set<size_t>::iterator cellIndex;
      //std::cerr << "ic1 " << ic1[i] << ": ";
      bool flag = 1;
      std::vector<int> k(dimension);
      if( dimension == 2) {
	for(k[1] = -numCell; k[1] <= numCell; k[1]+=numCell) {
	  for(k[0] = -1; k[0] <= 1; k[0]++){
	    cellIndices.insert(ic1[i]+ k[0] + k[1]);
	  }
	}
      } 
      else if( dimension == 3){
	for(k[2] = - numCell*numCell; k[2] <= numCell*numCell; k[2] += numCell*numCell) {
	  for(k[1] = -numCell; k[1] <= numCell; k[1]+=numCell) {
	    for(k[0] = -1; k[0] <= 1; k[0]++){
	      cellIndices.insert(ic1[i]+ k[0] + k[1] + k[2]);
	    }
	  }
	}
      }
      else {
	std::cerr << "NeighborhoodLatticeCigar::create() "
		  << "Wrong number of dimensions.\n";
	exit(-1);
      }
      for(cellIndex = cellIndices.begin(); cellIndex != cellIndices.end(); cellIndex++) {
	if( (*cellIndex >= 0 && *cellIndex <= totalNumCells) 
	    && latticeCell[*cellIndex].size()) { 
	  
	  for(size_t n = 0; n < latticeCell[*cellIndex].size(); n++){
	    
	    flag = 1;
	    size_t LCell = latticeCell[*cellIndex][n];
	    //std::cerr << LCell << " " << *cellIndex << "\n";
	    if(neighbor[i].size()){
	      for(size_t l = 0; l<neighbor[i].size(); l++) {
		if(LCell == neighbor[i][l])
		  flag = 0;
	      }
	    }
	    if(LCell > i && flag) {
	      //bool control = check(i,LCell, dimension, y);
	      //bool control = check2(y[i],y[LCell]);
	      bool control = 1;
	      if( control){
		//std::cerr << i << " and " << LCell <<"\n";
		//std::cerr << ic1[i] << " " << ic1[LCell] << "\n";
		neighbor[i].push_back(LCell );
		neighbor[LCell].push_back(i);
		numNeigh++;
	      }
	    }
	  }
	}
      }
    }
    
    for(size_t w = 0; w < numWall; ++w) {
      for(size_t m = 0; m < wallNeighborLattice[w].size(); ++m) {
	for(size_t n = 0; n < latticeCell[ wallNeighborLattice[w][m] ].size(); ++n) {
	  neighborWall[latticeCell[ wallNeighborLattice[w][m] ][n] ].push_back(w); 
	}
      }
    }
    
    // for ( size_t i = 0; i < latticeCell.size(); ++i) {
    // 			for ( size_t j=0; j < latticeCell[i].size(); ++j) 
    // 				std::cerr << latticeCell[i][j] << " ";
    // 			if (latticeCell[i].size())
    // 				std::cerr << "\n";
    // 		}
    
    
    //Copy the new neighbours into the cells (Caveat: Not a good solution)
    for(size_t i=0 ; i<N ; i++ ) {
      O.compartment(i).setNeighbor( neighbor[i] );
      O.compartment(i).setNeighborAtLevel( 1, neighborWall[i] );
    }
    
    //Set the previos time variable to current t value
    setPreviousTime(t);
    setNumNeighbor(numNeigh);
    
    return numNeigh;
  }
  return numNeighbor();
}

NeighborhoodLatticeCigar::
NeighborhoodLatticeCigar(std::vector<double> &paraValue, 
												 std::vector< std::vector<size_t> > 
												 &indValue ) 
{
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=6 ) {
    std::cerr << "NeighborhoodLatticeCigar::"
							<< "NeighborhoodLatticeCigar() "
							<< "Uses six parameters numCell, cellSize, updateFlag, "
							<< "dt_min, b and b_frac.\n";
    exit(0);
  }
  if( indValue.size() ) {
    std::cerr <<  "NeighborhoodLatticeCigar::"
							<<  "NeighborhoodLatticeCigar() "
							<< "No variable index is used.\n";
    exit(0);
  }

  // Check for consistency in parameter values
	if( paraValue[2]!=0.0 && paraValue[2]!=1.0 ) {
    std::cerr << "NeighborhoodLatticeCigar::"
							<< "NeighborhoodLatticeCigar() "
							<< "parameter(1) (updateFlag) must be 0 or 1.\n";
    exit(0);
  }

  // Set the variable values
  setId("neighborhoodLatticeCigar");
  setParameter(paraValue);  
  setVariableIndex(indValue);

	// Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "numCell";
  tmp[1] = "cellSize";
  tmp[2] = "updateFlag";
  tmp[3] = "dt_min";
  tmp[4] = "b";
  tmp[5] = "b_frac";
  setParameterId( tmp );
}

unsigned int NeighborhoodLatticeCigar::
create(Organism &O, std::vector< std::vector<double> > &y, double t ) 
{   
  if( O.numTopologyVariable() != 5 && 
      O.numTopologyVariable() != 7 &&
      O.numTopologyVariable() != 9 ) {
    std::cerr << "NeighborhoodLatticeCigar::create() "
							<< "Wrong number of topology variables.\n";
    exit(-1);
  }
  size_t N = O.numCompartment();
  std::vector< std::vector<size_t> > neighbor(N);
  unsigned int numNeigh=0;
  double cellSize = parameter(1);
  double cellSizeInverse = 1/cellSize;
  int numCell = static_cast<int>(parameter(0));
  size_t dimension=  static_cast<size_t>( (O.numTopologyVariable()-1)/2 );
  if(O.numTopologyVariable() == 9 ) 
		dimension = 2;
  int totalNumCells = static_cast<int>( std::pow( (double) numCell, (int) dimension) );
  //std::cerr << numCell << " " << cellSize << " " 
  //	    << dimension << " " << totalNumCells << "\n";
  std::vector<double> yc(dimension);
  std::vector< size_t > ic1( N );//given a compartment index i, gives lattice cell index
  std::vector< std::vector< size_t > > latticeCell(totalNumCells);
	
  //set up lattice, and place each cell point in a lattice cell
  for( size_t d = 0; d<dimension; d++){     
    int powNumc =  int( std::pow((double) numCell, (int) d) );
    for(size_t i = 0; i<N;i++){
      yc[d] = 0.5*(y[i][ d ]+y[i][ d+dimension ]);
      if( (yc[d] < cellSize || yc[d] > cellSize*(numCell-1)) ){
				std::cerr << "NeighborhoodLatticeCigar::create(): "
									<< "compartment is outside lattice\n";
				exit(-1);
      }
      ic1[i] +=  powNumc*static_cast<size_t>( yc[d]*cellSizeInverse );
    }
  }
  
  for(size_t i = 0; i<N;i++){
    latticeCell[ ic1[i] ].push_back(i);
  }
  
  //for each cell compartment, add the cells in neighboring lattice cells as neighbors
  for( size_t i=0 ; i<N ; i++ ){
    std::set<int> cellIndices;
    std::set<int>::iterator cellIndex;
    //std::cerr << "ic1 " << ic1[i] << ": ";
    bool flag = 1;
    std::vector<int> k(dimension);
    if( dimension == 2) {
      for(k[1] = -numCell; k[1] <= numCell; k[1]+=numCell) {				
				for(k[0] = -1; k[0] <= 1; k[0]++){
					int ind = ic1[i]+ k[0] + k[1];
					//if (ind < totalNumCells && ind >= 0)
						cellIndices.insert(ind);
				}
      }
    } 
    else if( dimension == 3){
      for(k[2] = - numCell*numCell; k[2] <= numCell*numCell; k[2] += numCell*numCell) {
				for(k[1] = -numCell; k[1] <= numCell; k[1]+=numCell) {
					for(k[0] = -1; k[0] <= 1; k[0]++){
						cellIndices.insert(ic1[i]+ k[0] + k[1] + k[2]);
						//if(ic2[i] != ic1[i]) {
						// cellIndices.insert(ic2[i]+ k[0] + k[1]+ k[2]);
						//    }
					}
				}
      }
    }
    else {
      std::cerr << "NeighborhoodLatticeCigar::create() "
								<< "Wrong number of dimensions.\n";
      exit(-1);
    }
    for(cellIndex = cellIndices.begin(); cellIndex != cellIndices.end(); ++cellIndex) {
      if( (*cellIndex >= 0 && *cellIndex <= totalNumCells) 
					&& latticeCell[*cellIndex].size()) { 
				for(size_t n = 0; n < latticeCell[*cellIndex].size(); n++){
					flag = 1;
					//std::cerr << *cellIndex << " " << n << "\n";
					size_t LCell = latticeCell[*cellIndex][n];
					if(neighbor[i].size()){
						for(size_t l = 0; l<neighbor[i].size(); l++) {
							if(LCell == neighbor[i][l])
								flag = 0;
						}
					}
					if(LCell > i && flag) {
						//bool control = check(i,LCell, dimension, y);
						//bool control = check2(y[i],y[LCell]);
						//std::cerr << ic1[i] << " " << *cellIndex << "\t" << i << " " << LCell << "\t" 
						//      << flag << " " << control << "\n";
						bool control = 1;
						if( control){
							neighbor[i].push_back(LCell );
							neighbor[LCell].push_back(i);
							numNeigh++;
						}
					}
				}
      }
    }
  }
  
  //Copy the new neighbours into the cells (Caveat: Not a good solution)
  for(size_t i=0 ; i<N ; i++ )
    O.compartment(i).setNeighbor( neighbor[i] );
	
  //Set the previos time variable to current t value
  setPreviousTime(t);
  setNumNeighbor(numNeigh);
	
  return numNeigh;
}

unsigned int NeighborhoodLatticeCigar::
update(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  if( parameter(2) == 1.0 && t-previousTime()>parameter(3) )
    return create(O,y,t);
  
  return numNeighbor();
}

bool NeighborhoodLatticeCigar::
check(size_t i, size_t j, size_t dimension,
			std::vector< std::vector<double> > &y ) 
{  
  double threshold = 4.0*parameter(4)*parameter(4)
    *parameter(5)*parameter(5);
  //List columns used
  std::vector<size_t> x1Col(dimension),x2Col(dimension);
  for( size_t d=0 ; d<dimension ; d++ ) {
    x1Col[d] = d;
    x2Col[d] = d+dimension;
  }
  std::vector<double> xc(dimension),nc(dimension),
    xcJ(dimension),ncJ(dimension);
  
  //Set cell parameters
  for( size_t d=0 ; d<dimension ; d++ ) {
    xc[d] = 0.5*(y[i][ x1Col[d] ]+y[i][ x2Col[d] ]);
    nc[d] = xc[d]-y[i][ x1Col[d] ];
  }
  double a2 = 0.0;
  for( size_t d=0 ; d<dimension ; d++ )
    a2 += nc[d]*nc[d];
  
  for( size_t d=0 ; d<dimension ; d++ ) {
    xcJ[d] = 0.5*(y[j][ x1Col[d] ]+y[j][ x2Col[d] ]);
    ncJ[d] = xcJ[d]-y[j][ x1Col[d] ];
  }
  double aJ2 = 0.0;
  for( size_t d=0 ; d<dimension ; d++ )
    aJ2 += ncJ[d]*ncJ[d];
	
  //Check distances between com
  //double distanceCOM=(xc[0]-xcJ[0])*(xc[0]-xcJ[0])+
  //(xc[1]-xcJ[1])*(xc[1]-xcJ[1]);
  //if( std::sqrt(distanceCOM)>std::sqrt(a2)+std::sqrt(aJ2)+1.1 ) return false;
	
  //Check distances from end points to neighbor line/endpoint
  std::vector<double> x1Pos(dimension),x2Pos(dimension),
    nSol(dimension);      
	
	//--------------------------------------------------
  //x1Cell
	//--------------------------------------------------
  for( size_t d=0 ; d<dimension ; d++ )
    x1Pos[d] = y[i][x1Col[d]];
  double t = 0.0;
  for( size_t d=0 ; d<dimension ; d++ )
    t += ncJ[d]*(x1Pos[d]-xcJ[d]);
  t /= aJ2;
  if( t<-1.0 ) {//outside neigh in x1 direction
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = y[j][x1Col[d]];
  }
  else if( t>1.0 ) {//outside neigh in x2 direction
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = y[j][x2Col[d]];
  }
  else {//on line
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = xcJ[d]+t*ncJ[d];
  }
  double dSol = 0.0;
  for( size_t d=0 ; d<dimension ; d++ ) {
    nSol[d] = x2Pos[d]-x1Pos[d];
    dSol += nSol[d]*nSol[d];
  }
  
  if( dSol<threshold ) return true;
  //dSol = std::sqrt( dSol );
	//--------------------------------------------------
  //x2Cell
	//--------------------------------------------------
  for( size_t d=0 ; d<dimension ; d++ )
    x1Pos[d] = y[i][x2Col[d]];
  t = 0.0;
  for( size_t d=0 ; d<dimension ; d++ )
    t += ncJ[d]*(x1Pos[d]-xcJ[d]);
  t /= aJ2;
  if( t<-1.0 ) {//outside neigh in x1 direction
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = y[j][x1Col[d]];
  }
  else if( t>1.0 ) {//outside neigh in x2 direction
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = y[j][x2Col[d]];
  }
  else {//on line
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = xcJ[d]+t*ncJ[d];
  }
  dSol = 0.0;
  for( size_t d=0 ; d<dimension ; d++ ) {
    nSol[d] = x2Pos[d]-x1Pos[d];
    dSol += nSol[d]*nSol[d];
  }
  if( dSol<threshold ) return true;
	
	//--------------------------------------------------
  //x1Neigh
	//--------------------------------------------------
  for( size_t d=0 ; d<dimension ; d++ )
    x1Pos[d] = y[j][x1Col[d]];
  t = 0.0;
  for( size_t d=0 ; d<dimension ; d++ )
    t += nc[d]*(x1Pos[d]-xc[d]);
  t /= a2;
  if( t<-1.0 ) {//outside cell in x1 direction
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = y[i][x1Col[d]];
  }
  else if( t>1.0 ) {//outside cell in x2 direction
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = y[i][x2Col[d]];
  }
  else {//on line
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = xc[d]+t*nc[d];
  }
  dSol = 0.0;
  for( size_t d=0 ; d<dimension ; d++ ) {
    nSol[d] = x2Pos[d]-x1Pos[d];
    dSol += nSol[d]*nSol[d];
  }
  if( dSol<threshold ) return true;
	//--------------------------------------------------
  //x2Neigh
	//--------------------------------------------------
  for( size_t d=0 ; d<dimension ; d++ )
    x1Pos[d] = y[j][x2Col[d]];
  t = 0.0;
  for( size_t d=0 ; d<dimension ; d++ )
    t += nc[d]*(x1Pos[d]-xc[d]);
  t /= a2;
  if( t<-1.0 ) {//outside cell in x1 direction
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = y[i][x1Col[d]];
  }
  else if( t>1.0 ) {//outside cell in x2 direction
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = y[i][x2Col[d]];
  }
  else {//on line
    for( size_t d=0 ; d<dimension ; d++ )
      x2Pos[d] = xc[d]+t*nc[d];
  }
  dSol = 0.0;
  for( size_t d=0 ; d<dimension ; d++ ) {
    nSol[d] = x2Pos[d]-x1Pos[d];
    dSol += nSol[d]*nSol[d];
  }
  if( dSol<threshold ) return true;
	
  return false;
}

bool NeighborhoodLatticeCigar::
check2(std::vector<double> &yi, std::vector<double> &yj ) {
  
  size_t dimension=2;    
  double threshold = 2*parameter(4)*parameter(5);
  threshold *= threshold;
  
  std::vector<double> u(dimension),v(dimension),w(dimension);
	
  for( size_t dim=0 ; dim<dimension ; ++dim ) {
    u[dim] = yi[dim+dimension]-yi[dim];
    v[dim] = yj[dim+dimension]-yj[dim];
    w[dim] = yi[dim]-yj[dim];
  }
  double a=0.0,b=0.0,c=0.0,d=0.0,e=0.0;
  for( size_t dim=0 ; dim<dimension ; ++dim ) {
    a += u[dim]*u[dim];
    b += u[dim]*v[dim];
    c += v[dim]*v[dim];
    d += u[dim]*w[dim];
    e += v[dim]*w[dim];
  } 
  double D=a*c-b*b;
  double sN,sD=D,tN,tD=D;
	
  if( D<1e-10 ) {
    //almost parallel, use yi[d] to prevent division by zero
    sN = 0.0;
    sD = 1.0;
    tN = e;
    tD = c;
  }
  else {
    sN = b*e-c*d;
    tN = a*e-b*d;
    if( sN<0.0 ) {
      sN = 0.0;
      tN = e;
      tD = c;
    }
    else if( sN>sD ) {
      sN = sD;
      tN = e+b;
      tD = c;
    }
  }
  if( tN<0.0 ) {
    tN = 0.0;
    if( -d<0 )
      sN = 0.0;
    else if( -d>a )
      sN = sD;
    else {
      sN = -d;
      sD = a;
    }
  }
  else if( tN>tD ) {
    tN = tD;
    if( (-d+b)<0.0 )
      sN = 0.0;
    else if( (-d+b)>a )
      sN = sD;
    else {
      sN = -d+b;
      sD = a;
    }
  }
  double sc = sN/sD; 
  double tc = tN/tD;
  
  std::vector<double> dP(dimension);
  for( size_t dim=0 ; dim<dimension ; ++dim )
    dP[dim] = w[dim]+sc*u[dim]-tc*v[dim];
  //double dist = dP[0]*dP[0]+dP[1]*dP[1];
  
  if( dP[0]*dP[0]+dP[1]*dP[1] > threshold )
    return false;
  return true;  
}

Sphere::NeighborhoodDistance::
NeighborhoodDistance(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > 
		     &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=3 ) {
    std::cerr << "Sphere::NeighborhoodDistance::"
	      << "NeighborhoodDistance() "
	      << "Uses three parameters r_frac, upDateFlag and dt_min." << std::endl;
    exit(EXIT_FAILURE);
  }
  if( indValue.size() ) {
    std::cerr << "Sphere::NeighborhoodDistance::"
	      << "NeighborhoodDistance() "
	      << "No variable index is used." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // Check for consistency in parameter values
  if( paraValue[0]<0.0 ) {
    std::cerr << "Sphere::NeighborhoodDistance::"
	      << "NeighborhoodDistance() "
	      << "parameter(0) (r_frac) must be larger than zero.\n";
    exit(EXIT_FAILURE);
  }
  if( paraValue[1]!=0.0 && paraValue[1]!=1.0 ) {
    std::cerr << "Sphere::NeighborhoodDistance::"
	      << "NeighborhoodDistance() "
	      << "parameter(1) (updateFlag) must be 0 or 1.\n";
    exit(EXIT_FAILURE);
  }

  // Set the variable values
  setId("Sphere::NeighborhoodDistance");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "r_frac";
  tmp[1] = "updateFlag";
  tmp[2] = "dt_min";
  setParameterId( tmp );
}

unsigned int Sphere::NeighborhoodDistance::
create(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  size_t N = O.numCompartment();
  std::vector< std::vector<size_t> > neighbor(N);
  double d,r;
  unsigned int numNeigh=0;
  
  size_t rIndex = O.numTopologyVariable()-1;
  size_t eStart = 0;
  
  for(size_t i=0 ; i<N ; i++ )
    for(size_t j=i+1 ; j<N ; j++ ) {
      r = y[i][rIndex] + y[j][rIndex];
      d=0.;
      for( size_t dim=eStart ; dim<rIndex ; dim++ )
	d += (y[i][dim]-y[j][dim])*(y[i][dim]-y[j][dim]);
      d = std::sqrt( d );
      
      if( d<=r*parameter(0) ) {
        //Add neighbor
        neighbor[i].push_back(j);
        neighbor[j].push_back(i);
        numNeigh++;
      }
    }
  //Copy the new neighbours into the cells (Caveat: Not a good solution)
  for(size_t i=0 ; i<N ; i++ )
    O.compartment(i).setNeighbor( neighbor[i] );
	
  //Set the previosu time variable to current t value
  setPreviousTime(t);
  setNumNeighbor(numNeigh);
	
  return numNeigh;
}

unsigned int Sphere::NeighborhoodDistance::
update(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  if( parameter(1) == 1.0 && t-previousTime()>parameter(2) )
    return create(O,y,t);
	
  return numNeighbor();
}

NeighborhoodDistanceSpherePeriodic::
NeighborhoodDistanceSpherePeriodic(std::vector<double> &paraValue, 
																	 std::vector< std::vector<size_t> > 
																	 &indValue ) 
{
  // Do some checks on the parameters and variable indices
  if( paraValue.size()!=3 ) {
    std::cerr << "NeighborhoodDistanceSpherePeriodic::"
							<< "NeighborhoodDistanceSpherePeriodic() "
							<< "Uses three parameters r_frac, upDateFlag and dt_min.\n";
    exit(0);
  }
  if( indValue.size() >1 ) {
    std::cerr << "NeighborhoodDistanceSpherePeriodic::"
							<< "NeighborhoodDistanceSpherePeriodic() "
							<< "Max one set of indeces is used.\n";
    exit(0);
  }

  // Check for consistency in parameter values
  if( paraValue[0]<0.0 ) {
    std::cerr << "NeighborhoodDistanceSpherePeriodic::"
							<< "NeighborhoodDistanceSpherePeriodic() "
							<< "parameter(0) (r_frac) must be larger than zero.\n";
    exit(0);
  }
  if( paraValue[1]!=0.0 && paraValue[1]!=1.0 ) {
    std::cerr << "NeighborhoodDistanceSpherePeriodic::"
							<< "NeighborhoodDistanceSpherePeriodic() "
							<< "parameter(1) (updateFlag) must be 0 or 1.\n";
    exit(0);
  }

  // Set the variable values
  setId("neighborhoodDistanceSpherePeriodic");
  setParameter(paraValue);  
  setVariableIndex(indValue);

  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "r_frac";
  tmp[1] = "updateFlag";
  tmp[2] = "dt_min";
  setParameterId( tmp );
}

unsigned int NeighborhoodDistanceSpherePeriodic::
create(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  size_t N = O.numCompartment();
  std::vector< std::vector<size_t> > neighbor(N);
  double d,r;
  unsigned int numNeigh=0;
  
  size_t rIndex = O.numTopologyVariable()-1;
  size_t eStart = 0;
  
  for(size_t i=0 ; i<N ; i++ )
    for(size_t j=i+1 ; j<N ; j++ ) {
      r = y[i][rIndex] + y[j][rIndex];
      d=0.;
      for( size_t dim=eStart ; dim<rIndex ; dim++ )
				d += (y[i][dim]-y[j][dim])*(y[i][dim]-y[j][dim]);
      d = std::sqrt( d );
      
      if( d<=r*parameter(0) ) {
        //Add neighbor
        neighbor[i].push_back(j);
        neighbor[j].push_back(i);
        numNeigh++;
      }
    }
	
  //Add periodic boundary conditions
  if( numVariableIndexLevel() && numVariableIndex(0) ) {
    for( size_t k=0 ; k<numVariableIndex(0) ; k++ ) {
      size_t col=variableIndex(0,k);
      if( col>=rIndex ) {
				std::cerr << "NeighborhoodDistanceSpherePeriodic::create Warning "
									<< "Periodic coordinate larger than dimension\n";
      }
      else {
				addPeriodicBoundaries(y,neighbor,col,rIndex);
      }
    }
  }
  //Copy the new neighbours into the cells 
	//(Caveat: Not a good solution)
  for(size_t i=0 ; i<N ; i++ )
    O.compartment(i).setNeighbor( neighbor[i] );
  
  //Set the previous time variable to current t value
  setPreviousTime(t);
  
  return numNeigh;
}

unsigned int NeighborhoodDistanceSpherePeriodic::
update(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  if( parameter(1) == 1.0 && t-previousTime()>parameter(2) )
    return create(O,y,t);
	
  return numNeighbor();
}

unsigned int NeighborhoodDistanceSpherePeriodic::
addPeriodicBoundaries(std::vector< std::vector<double> > &y, 
											std::vector< std::vector<size_t> > &neighbor,
											size_t dIndex,size_t rIndex) 
{  
  if( y.empty() ) return 0;
  assert( y.size()==neighbor.size() );
  assert( y[0].size()>dIndex && y[0].size()>rIndex );
  
  // Find maximal and minimal positions in dimension dIndex
  double max=y[0][dIndex]+y[0][rIndex],
    min = y[0][dIndex]-y[0][rIndex];
  for( size_t i=1 ; i<y.size() ; i++ ) {
    if( y[i][dIndex]+y[i][rIndex] > max )
      max = y[i][dIndex]+y[i][rIndex];
    if( y[i][dIndex]-y[i][rIndex] < min )
      min = y[i][dIndex]-y[i][rIndex];
  }
  
  // Extract potential cells for periodic boundary
  std::vector<size_t> potentialLow,potentialHigh;
  for( size_t i=0 ; i<y.size() ; i++ ) {
    if( y[i][dIndex]+2*y[i][rIndex] >= max )
      potentialHigh.push_back(i);
    else if( y[i][dIndex]-2*y[i][rIndex] <= min )
      potentialLow.push_back(i);
  }
  
  // Move low cells above top ones
  double delta=max-min;
  for( size_t k=0 ; k<potentialLow.size() ; k++ ) {
    size_t i=potentialLow[k];
    y[i][dIndex] += delta;
  }
  
  // Find (between potLow and potHigh) and add new neighbors
  double rFactor=parameter(0);
  unsigned int numAddedNeighbors=0;
  for( size_t k=0 ; k<potentialLow.size() ; k++ ) {
    size_t i=potentialLow[k];
    for( size_t m=0 ; m<potentialHigh.size() ; m++ ) {
      size_t j=potentialHigh[m];
      double distance=0.0;
      for( size_t d=0 ; d<rIndex ; d++ )
				distance += (y[i][d]-y[j][d])*(y[i][d]-y[j][d]);
      distance = sqrt(distance);
      if( distance<=rFactor*(y[i][rIndex]+y[i][rIndex]) ) {
				neighbor[i].push_back(j);
				neighbor[j].push_back(i);		 
				numAddedNeighbors++;
      }
    }
  }
  // Move the cells back...
  for( size_t k=0 ; k<potentialLow.size() ; k++ ) {
    size_t i=potentialLow[k];
    y[i][dIndex] -= delta;
  } 
  
  return numAddedNeighbors;
}

NeighborhoodDistanceSphereWithWalls::
NeighborhoodDistanceSphereWithWalls(std::vector<double> &paraValue, 
																		std::vector< std::vector<size_t> > 
																		&indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=3 && paraValue.size() !=5 ) {
    std::cerr << "NeighborhoodDistanceSphereWithWalls::"
							<< "NeighborhoodDistanceSphereWithWalls() "
							<< "Uses 3 (or 5) parameters r_frac, upDateFlag and "
							<< "dt_min (and cellVolume and wallWidth).\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 2 ) {
    std::cerr << "NeighborhoodDistanceSphereWithWalls::"
							<< "NeighborhoodDistanceSphereWithWalls() "
							<< "Two variable indeces used, cellMarker wallMarker.\n";
    exit(0);
  }

  // Check for consistency in parameter values
  if( paraValue[0]<0.0 ) {
    std::cerr << "NeighborhoodDistanceSphereWithWalls::"
							<< "NeighborhoodDistanceSphereWithWalls() "
							<< "parameter(0) (r_frac) must be larger than zero.\n";
    exit(0);
  }
  if( paraValue[1]!=0.0 && paraValue[1]!=1.0 ) {
    std::cerr << "NeighborhoodDistanceSphereWithWalls::"
							<< "NeighborhoodDistanceSphereWithWalls() "
							<< "parameter(1) (updateFlag) must be 0 or 1.\n";
    exit(0);
  }

	// Set the variable values
  setId("neighborhoodDistanceSphereWithWalls");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "r_frac";
  tmp[1] = "updateFlag";
  tmp[2] = "dt_min";
  if( numParameter() == 5 ) {
    tmp[3] = "cellVolume";
    tmp[4] = "wallWidth";
  }
  setParameterId( tmp );
}

unsigned int NeighborhoodDistanceSphereWithWalls::
create(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{
  size_t N = O.numCompartment();
  assert( N==y.size() );
  std::vector< std::vector<size_t> > cellCellNeighbor(N),
    cellWallNeighbor(N),wallWallNeighbor;
  double d,r;
  unsigned int numNeigh=0;
  size_t cellIndex = variableIndex(0,0);
  size_t wallIndex = variableIndex(0,1);
  
  size_t rIndex = O.numTopologyVariable()-1;
  if( rIndex != 2 ) {
    std::cerr << "NeighborhoodDistanceSphereWithWalls::create() "
							<< "Only allowed for 2D sofar...\n";
    exit(0);
  }
  size_t eStart = 0;
  
  for(size_t i=0 ; i<N ; i++ ) {
    if( y[i][cellIndex] != 1.0 || y[i][wallIndex] != 0.0 ) {
      std::cerr << "NeighborhoodDistanceSphereWithWalls::create() "
								<< "Cell or Wall marker wrong in cell " << i << "\n";
      exit(0);
    }
    for(size_t j=i+1 ; j<N ; j++ ) {
      r = y[i][rIndex] + y[j][rIndex];
      d=0.;
      for( size_t dim=eStart ; dim<rIndex ; dim++ )
				d += (y[i][dim]-y[j][dim])*(y[i][dim]-y[j][dim]);
      d = std::sqrt( d );
      
      if( d<=r*parameter(0) ) {
        //Add cell-cell neighbors
        cellCellNeighbor[i].push_back(j);
        cellCellNeighbor[j].push_back(i);
        numNeigh++;
				//Add wall compartment
				// y
				size_t iWall = y.size();
				y.push_back( y[i] );
				y[iWall][0] = 0.5*(y[i][0]+y[j][0]);
				y[iWall][1] = 0.5*(y[i][1]+y[j][1]);
				y[iWall][cellIndex]=0.0;
				y[iWall][wallIndex]=1.0;
				//O.compartment
				Compartment tmpCompartment( O.compartment(i) );
				O.addCompartment( tmpCompartment );
				O.compartment(iWall).setIndex(iWall);
				O.compartment(iWall).setVariable(0,y[iWall][0]);
				O.compartment(iWall).setVariable(1,y[iWall][1]);
				O.compartment(iWall).setVariable(cellIndex,y[iWall][cellIndex]);
				O.compartment(iWall).setVariable(wallIndex,y[iWall][wallIndex]);
				//Add cell-wall neighbors
				cellWallNeighbor.resize( y.size() );
				cellWallNeighbor[i].push_back(iWall);
				cellWallNeighbor[j].push_back(iWall);
				cellWallNeighbor[iWall].push_back(i);
				cellWallNeighbor[iWall].push_back(j);
      }
    }
  }
  assert(O.numCompartment()==y.size());
  //Add extra boundary walls
  //Not yet...
  
  // Get wall-wall neighbors by angular sorting from cells
	//--------------------------------------------------------------------
  wallWallNeighbor.resize(y.size()-N);
  double PI=3.141592654;
  for(size_t i=0 ; i<N ; i++ ) {
    std::vector<double> phi( cellCellNeighbor[i].size() );
    for( size_t k=0 ; k<cellCellNeighbor[i].size() ; k++ ) {
      size_t j = cellCellNeighbor[i][k];
      double dx = y[j][0]-y[i][0];
      double dy = y[j][1]-y[i][1];
      if( dx==0.0 ) {
				if( dy>0 )
					phi[k] = 0.5*PI;
				else
					phi[k] = 1.5*PI;
      }
      else if( dy==0.0 ) {
				if( dx>0.0 ) 
					phi[k]=0.0;
				else
					phi[k]=PI;
      }
      else {
				phi[k] = std::atan( std::fabs(dy/dx) );
				if( dx>0 && dy<0 )
					phi[k] = 2*PI-phi[k];
				else if( dx<0 && dy>0 )
					phi[k] = PI-phi[k];
				else if( dx<0 && dy<0 )
					phi[k] = PI+phi[k];
      }
    }
    // Sort after angles
    std::vector<size_t> kVector( phi.size() );
    for( size_t k=0 ; k<kVector.size() ; k++ )
      kVector[k] = k;
    
    for( size_t k=0 ; k<phi.size() ; k++ )
      for( size_t kk=k+1 ; kk<phi.size() ; kk++ )
				if( phi[kk]>phi[k] ) {
					double tmp=phi[kk];
					phi[kk]=phi[k];
					phi[k]=tmp;
					size_t tmpI=kVector[kk];
					kVector[kk]=kVector[k];
					kVector[k]=tmpI;
				}
    
    // Add the wall-wall neighbors    
    for( size_t k=0 ; k<kVector.size()-1 ; k++ ) {
      double deltaPhi = std::fabs(phi[k+1]-phi[k]);
      if( deltaPhi>PI ) deltaPhi = 2*PI-deltaPhi;
      if( deltaPhi<(0.4*PI) ) {
				wallWallNeighbor[ cellWallNeighbor[i][kVector[k]]-N ].
					push_back( cellWallNeighbor[i][kVector[k+1]] );
				wallWallNeighbor[ cellWallNeighbor[i][kVector[k+1]]-N ].
					push_back( cellWallNeighbor[i][kVector[k]] );
      }
    }      
    if( kVector.size()>2 ) {
      size_t k=kVector.size()-1;
      double deltaPhi = std::fabs(phi[0]-phi[k]);
      if( deltaPhi>PI ) deltaPhi = 2*PI-deltaPhi;
      if( deltaPhi<(0.4*PI) ) {
				wallWallNeighbor[ cellWallNeighbor[i][kVector[k]]-N ].
					push_back( cellWallNeighbor[i][kVector[0]] );
				wallWallNeighbor[ cellWallNeighbor[i][kVector[0]]-N ].
					push_back( cellWallNeighbor[i][kVector[k]] );
      }
    }
  }
  // Copy the new neighbours into the cells 
	// (Caveat: Not a good solution)
  // Cell compartments
  for(size_t i=0 ; i<N ; i++ ) {
    // Adjacent compartments (cell-wall) in level 0
    O.compartment(i).setNeighbor( cellWallNeighbor[i] );
    // Cell-cell neighbors at level 1
    O.compartment(i).setNeighborAtLevel(1,cellCellNeighbor[i] );
  }
  // Wall compartments
  for(size_t i=N ; i<y.size() ; i++ ) {
    std::vector<size_t> tmpV(cellWallNeighbor[i]);
    tmpV.insert( tmpV.end(), wallWallNeighbor[i-N].begin(),
		 wallWallNeighbor[i-N].end() );
    O.compartment(i).setNeighbor( tmpV );
  }
  
  // If spatial variables are defined add those
  if( numParameter()==5 ) {
    // Volumes: 
    // cells: volume=cellVolume=parameter(3)
    for(size_t i=0 ; i<N ; i++ ) {
      y[i][rIndex] = parameter(3);
      O.compartment(i).setVariable(rIndex,y[i][rIndex]);
    }
    // walls: volume=sqrt(cellVolume)*wallWidth=sqrt(param(3))*param(4)
    for(size_t i=N ; i<y.size() ; i++ ) {
      y[i][rIndex] = std::sqrt(parameter(3))*parameter(4);
      O.compartment(i).setVariable(rIndex,y[i][rIndex]);
    }
    // NeighborAreas:
    // cellWall: area=sqrt(cellVolume)=sqrt(parameter(3))
    // wallWall: area=wallWidth=parameter(4)
    double cellWallArea=std::sqrt(parameter(3));
    double wallWallArea=parameter(4);
    for(size_t i=0 ; i<N ; i++ ) {
      std::vector<double> tmpV(O.compartment(i).numNeighbor(),cellWallArea);
      O.compartment(i).setNeighborArea(tmpV);
    }
    for(size_t i=N ; i<y.size() ; i++ ) {
      std::vector<double> tmpV(cellWallNeighbor[i].size(),cellWallArea);
      tmpV.insert( tmpV.end(), wallWallNeighbor[i-N].size(),wallWallArea);
      O.compartment(i).setNeighborArea(tmpV);
    }
  }
  //Set the previous time variable to current t value
  setPreviousTime(t);
  
  return numNeigh;
}

unsigned int NeighborhoodDistanceSphereWithWalls::
update(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  if( parameter(1) == 1.0 && t-previousTime()>parameter(2) ) {
    std::cerr << "NeighborhoodDistanceSphereWithWalls::update "
							<< " not allowed yet...\n";
    exit(0);
    return create(O,y,t);
  }
	
  return numNeighbor();
}

NeighborhoodDistanceSphereBud::
NeighborhoodDistanceSphereBud(std::vector<double> &paraValue, 
															std::vector< std::vector<size_t> > 
															&indValue ) 
{
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=3 ) {
    std::cerr << "NeighborhoodDistanceSphereBud::"
							<< "NeighborhoodDistanceSphereBud() "
							<< "Uses three parameters r_frac, upDateFlag and dt_min.\n";
    exit(0);
  }
  if( indValue.size() ) {
    std::cerr << "NeighborhoodDistanceSphereBud::"
							<< "NeighborhoodDistanceSphereBud() "
							<< "No variable index is used.\n";
    exit(0);
  }

  // Check for consistency in parameter values
  if( paraValue[0]<0.0 ) {
    std::cerr << "NeighborhoodDistanceSphereBud::"
							<< "NeighborhoodDistanceSphereBud() "
							<< "parameter(0) (r_frac) must be larger than zero.\n";
    exit(0);
  }
  if( paraValue[1]!=0.0 && paraValue[1]!=1.0 ) {
    std::cerr << "NeighborhoodDistanceSphereBud::"
							<< "NeighborhoodDistanceSphereBud() "
							<< "parameter(1) (updateFlag) must be 0 or 1.\n";
    exit(0);
  }  

  // Set the variable values
  setId("neighborhoodDistanceSphereBud");
  setParameter(paraValue);  
  setVariableIndex(indValue);

	//Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "r_frac";
  tmp[1] = "updateFlag";
  tmp[2] = "dt_min";
  setParameterId( tmp );
}

unsigned int NeighborhoodDistanceSphereBud::
create(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  if( O.numTopologyVariable() != 4 && O.numTopologyVariable() != 6 && 
      O.numTopologyVariable() !=8 ) {
    std::cerr << "NeighborhoodDistanceSphereBud::create() "
							<< "Wrong number of topology variables.\n";
    exit(-1);
  }
  size_t dimension = static_cast<size_t>( (O.numTopologyVariable()-2)/2 );
  
  size_t N = O.numCompartment();
  std::vector< std::vector<size_t> > neighbor(N);
  double d,r;
  unsigned int numNeigh=0;
  
  size_t r1Index = O.numTopologyVariable()-2;
  size_t r2Index = r1Index+1;
  size_t x1Index = 0;
  size_t x2Index = dimension;
  
  for(size_t i=0 ; i<N ; i++ )
    for(size_t j=i+1 ; j<N ; j++ ) {
      unsigned int neighborAdd=0;
      //11
      r = y[i][r1Index] + y[j][r1Index];
      d = 0.0;
      for( size_t dim=0 ; dim<dimension ; dim++ )
				d += (y[i][x1Index+dim]-y[j][x1Index+dim])*
					(y[i][x1Index+dim]-y[j][x1Index+dim]);
      d = std::sqrt( d );
      if( d<=r*parameter(0) ) {
				neighborAdd++;
      }
      else {
				//12
				r = y[i][r1Index] + y[j][r2Index];
				d = 0.0;
				for( size_t dim=0 ; dim<dimension ; dim++ )
					d += (y[i][x1Index+dim]-y[j][x2Index+dim])*
						(y[i][x1Index+dim]-y[j][x2Index+dim]);
				d = std::sqrt( d );
				if( d<=r*parameter(0) ) {
					neighborAdd++;
				}
				else {
					//21
					r = y[i][r2Index] + y[j][r1Index];
					d = 0.0;
					for( size_t dim=0 ; dim<dimension ; dim++ )
						d += (y[i][x2Index+dim]-y[j][x1Index+dim])*
							(y[i][x2Index+dim]-y[j][x1Index+dim]);
					d = std::sqrt( d );
					if( d<=r*parameter(0) ) {
						neighborAdd++;
					}
					else {
						//22
						r = y[i][r2Index] + y[j][r2Index];
						d = 0.0;
						for( size_t dim=0 ; dim<dimension ; dim++ )
							d += (y[i][x2Index+dim]-y[j][x2Index+dim])*
								(y[i][x2Index+dim]-y[j][x2Index+dim]);
						d = std::sqrt( d );
						if( d<=r*parameter(0) ) {
							neighborAdd++;
						}
					}
				}
      }
      if( neighborAdd ) {
				//Add neighbor
				neighbor[i].push_back(j);
				neighbor[j].push_back(i);
				numNeigh++;
      }
    }
  //Copy the new neighbours into the cells (Caveat: Not a good solution)
  for(size_t i=0 ; i<N ; i++ )
    O.compartment(i).setNeighbor( neighbor[i] );
  
  //Set the previous time variable to current t value
  setPreviousTime(t);
  setNumNeighbor(numNeigh);
  
  return numNeigh;
}

unsigned int NeighborhoodDistanceSphereBud::
update(Organism &O, std::vector< std::vector<double> > &y,
       double t ) {
  
  if( parameter(1) == 1.0 && t-previousTime()>parameter(2) )
    return create(O,y,t);
	
  return numNeighbor();
}

NeighborhoodDistanceEllipse::
NeighborhoodDistanceEllipse(std::vector<double> &paraValue, 
														std::vector< std::vector<size_t> > 
														&indValue ) 
{
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=3 ) {
    std::cerr << "NeighborhoodDistanceEllipse::"
							<< "NeighborhoodDistanceEllipse() "
							<< "Uses three parameters r_frac, upDateFlag and dt_min.\n";
    exit(0);
  }
  if( indValue.size() ) {
    std::cerr << "NeighborhoodDistanceEllipse::"
							<< "NeighborhoodDistanceEllipse() "
							<< "No variable index is used.\n";
    exit(0);
  }

  // Check for consistency in parameter values
  if( paraValue[0]<0.0 ) {
    std::cerr << "NeighborhoodDistanceEllipse::"
							<< "NeighborhoodDistanceEllipse() "
							<< "parameter(0) (r_frac) must be larger than zero.\n";
    exit(0);
  }
  if( paraValue[1]!=0.0 && paraValue[1]!=1.0 ) {
    std::cerr << "NeighborhoodDistanceEllipse::"
							<< "NeighborhoodDistanceEllipse() "
							<< "parameter(1) (updateFlag) must be 0 or 1.\n";
    exit(0);
  }
  
  //Set the variable values
  setId("neighborhoodDistanceEllipse");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "r_frac";
  tmp[1] = "updateFlag";
  tmp[2] = "dt_min";
  setParameterId( tmp );
}

unsigned int NeighborhoodDistanceEllipse::
create(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  if( O.numTopologyVariable() != 5 ) {
    std::cerr << "NeighborhoodDistanceEllipse::create() "
							<< "Wrong number of topology variables.\n";
    exit(-1);
  }
  size_t N = O.numCompartment(); 
  // List columns used
  size_t xaCol=0,yaCol=1,xbCol=2,ybCol=3;
  // Define additional center points and variables
  size_t dimension=2;
  double d,r,ai,aj;
  unsigned int numNeigh=0;
  std::vector< std::vector<size_t> > neighbor(N);
  std::vector<double> ci(dimension),cj(dimension);
  
  for(size_t i=0 ; i<N ; i++ ) {
    ci[0] = 0.5*( y[i][xaCol]+y[i][xbCol] );
    ci[1] = 0.5*( y[i][yaCol]+y[i][ybCol] );
    ai = 0.5*std::sqrt( (y[i][xaCol]-y[i][xbCol])*
												(y[i][xaCol]-y[i][xbCol])+
												(y[i][yaCol]-y[i][ybCol])*
												(y[i][yaCol]-y[i][ybCol]) );
    for(size_t j=i+1 ; j<N ; j++ ) {
      cj[0] = 0.5*( y[j][xaCol]+y[j][xbCol] );
      cj[1] = 0.5*( y[j][yaCol]+y[j][ybCol] );
      aj = 0.5*std::sqrt( (y[j][xaCol]-y[j][xbCol])*
													(y[j][xaCol]-y[j][xbCol])+
													(y[j][yaCol]-y[j][ybCol])*
													(y[j][yaCol]-y[j][ybCol]) );
      r = ai+aj;
      d = std::sqrt( (ci[0]-cj[0])*(ci[0]-cj[0])+ 
										 (ci[1]-cj[1])*(ci[1]-cj[1]) );
      
      if( d<=r*parameter(0) ) {
        // Add neighbor
        neighbor[i].push_back(j);
        neighbor[j].push_back(i);
        numNeigh++;
      }
    }
  }
  // Copy the new neighbours into the cells (Caveat: Not a good solution)
  for(size_t i=0 ; i<N ; i++ )
    O.compartment(i).setNeighbor( neighbor[i] );
  
  // Set the previosu time variable to current t value
  setPreviousTime(t);
	
  return numNeigh;
}

unsigned int NeighborhoodDistanceEllipse::
update(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{
  if( parameter(1) == 1.0 && t-previousTime()>parameter(2) )
    return create(O,y,t);
	
  return numNeighbor();
}

NeighborhoodDistanceCigar::
NeighborhoodDistanceCigar(std::vector<double> &paraValue, 
													std::vector< std::vector<size_t> > &indValue ) 
{ 
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=4 ) {
    std::cerr << "NeighborhoodDistanceCigar::"
							<< "NeighborhoodDistanceCigar() "
							<< "Uses four parameters "
							<< "b, b_frac, updateFlag and dt_min.\n";
    exit(0);
  }
  if( indValue.size() ) {
    std::cerr << "NeighborhoodDistanceCigar::"
							<< "NeighborhoodDistanceCigar() "
							<< "No variable index is used.\n";
    exit(0);
  }

  //Check for consistency in parameter values
  if( paraValue[0]<0.0 ) {
    std::cerr << "NeighborhoodDistanceCigar::"
							<< "NeighborhoodDistanceCigar() "
							<< "parameter(0) (b) must be larger than zero.\n";
    exit(0);
  }
  if( paraValue[1]<0.0 ) {
    std::cerr << "NeighborhoodDistanceCigar::"
							<< "NeighborhoodDistanceCigar() "
							<< "parameter(1) (b_frac) must be larger than zero.\n";
    exit(0);
  }
  if( paraValue[2]!=0.0 && paraValue[2]!=1.0 ) {
    std::cerr << "NeighborhoodDistanceCigar::"
	      << "NeighborhoodDistanceCigar() "
							<< "parameter(2) (updateFlag) must be 0 or 1.\n";
    exit(0);
  }
  
  //Set the variable values
  setId("neighborhoodDistanceCigar");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "b";
  tmp[1] = "b_frac";
  tmp[2] = "updateFlag";
  tmp[3] = "dt_min";
  setParameterId( tmp );
}

unsigned int NeighborhoodDistanceCigar::
create(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  if( O.numTopologyVariable() != 5 && 
      O.numTopologyVariable() != 7 &&
      O.numTopologyVariable() != 9 ) {
    std::cerr << "NeighborhoodDistanceCigar::create() "
							<< "Wrong number of topology variables.\n";
    exit(-1);
  }
  size_t dimension=static_cast<size_t>( (O.numTopologyVariable()-1)/2 );
  if(O.numTopologyVariable() == 9 ) 
		dimension = 2;
  size_t N = O.numCompartment(); 
  //List columns used
  std::vector<size_t> x1Col(dimension),x2Col(dimension);
  for( size_t d=0 ; d<dimension ; d++ ) {
    x1Col[d] = d;
    x2Col[d] = d+dimension;
  }
  std::vector<double> xc(dimension),nc(dimension),
    xcJ(dimension),ncJ(dimension);
  std::vector< std::vector<size_t> > neighbor(y.size());
  unsigned int numNeigh=0;
  for(size_t i=0 ; i<N ; i++ ) {
    //Set cell parameters
    for( size_t d=0 ; d<dimension ; d++ ) {
      xc[d] = 0.5*(y[i][ x1Col[d] ]+y[i][ x2Col[d] ]);
      nc[d] = xc[d]-y[i][ x1Col[d] ];
    }
    double a = 0.0;
    for( size_t d=0 ; d<dimension ; d++ )
      a += nc[d]*nc[d];
    a = std::sqrt( a );
    
    for(size_t j=i+1 ; j<N ; j++ ) {
      for( size_t d=0 ; d<dimension ; d++ ) {
				xcJ[d] = 0.5*(y[j][ x1Col[d] ]+y[j][ x2Col[d] ]);
				ncJ[d] = xcJ[d]-y[j][ x1Col[d] ];
      }
      double aJ = 0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				aJ += ncJ[d]*ncJ[d];
      aJ = std::sqrt( aJ );
			
      //Check distances from end points to neighbor line/endpoint
      std::vector<double> x1Pos(dimension),x2Pos(dimension),
				nSol(dimension);      
      double dMin;
      //x1Cell
      //////////////////////////////////////////////////
      for( size_t d=0 ; d<dimension ; d++ )
				x1Pos[d] = y[i][x1Col[d]];
      double t = 0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				t += ncJ[d]*(x1Pos[d]-xcJ[d]);
      t /= aJ*aJ;
      if( t<-1.0 ) {//outside neigh in x1 direction
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = y[j][x1Col[d]];
      }
      else if( t>1.0 ) {//outside neigh in x2 direction
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = y[j][x2Col[d]];
      }
      else {//on line
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = xcJ[d]+t*ncJ[d];
      }
      double dSol = 0.0;
      for( size_t d=0 ; d<dimension ; d++ ) {
				nSol[d] = x2Pos[d]-x1Pos[d];
				dSol += nSol[d]*nSol[d];
      }
      dSol = std::sqrt( dSol );
      dMin = dSol;
      //x2Cell
      //////////////////////////////////////////////////
      for( size_t d=0 ; d<dimension ; d++ )
				x1Pos[d] = y[i][x2Col[d]];
      t = 0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				t += ncJ[d]*(x1Pos[d]-xcJ[d]);
      t /= aJ*aJ;
      if( t<-1.0 ) {//outside neigh in x1 direction
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = y[j][x1Col[d]];
      }
      else if( t>1.0 ) {//outside neigh in x2 direction
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = y[j][x2Col[d]];
      }
      else {//on line
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = xcJ[d]+t*ncJ[d];
      }
      dSol = 0.0;
      for( size_t d=0 ; d<dimension ; d++ ) {
				nSol[d] = x2Pos[d]-x1Pos[d];
				dSol += nSol[d]*nSol[d];
      }
      dSol = std::sqrt( dSol );
      dMin = dSol<dMin ? dSol : dMin;
      //x1Neigh
      //////////////////////////////////////////////////
      for( size_t d=0 ; d<dimension ; d++ )
				x1Pos[d] = y[j][x1Col[d]];
      t = 0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				t += nc[d]*(x1Pos[d]-xc[d]);
      t /= a*a;
      if( t<-1.0 ) {//outside cell in x1 direction
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = y[i][x1Col[d]];
      }
      else if( t>1.0 ) {//outside cell in x2 direction
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = y[i][x2Col[d]];
      }
      else {//on line
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = xc[d]+t*nc[d];
      }
      dSol = 0.0;
      for( size_t d=0 ; d<dimension ; d++ ) {
				nSol[d] = x2Pos[d]-x1Pos[d];
				dSol += nSol[d]*nSol[d];
      }
      dSol = std::sqrt( dSol );
      dMin = dSol<dMin ? dSol : dMin;
      //x2Neigh
      //////////////////////////////////////////////////
      for( size_t d=0 ; d<dimension ; d++ )
				x1Pos[d] = y[j][x2Col[d]];
      t = 0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				t += nc[d]*(x1Pos[d]-xc[d]);
      t /= a*a;
      if( t<-1.0 ) {//outside cell in x1 direction
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = y[i][x1Col[d]];
      }
      else if( t>1.0 ) {//outside cell in x2 direction
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = y[i][x2Col[d]];
      }
      else {//on line
				for( size_t d=0 ; d<dimension ; d++ )
					x2Pos[d] = xc[d]+t*nc[d];
      }
      dSol = 0.0;
      for( size_t d=0 ; d<dimension ; d++ ) {
				nSol[d] = x2Pos[d]-x1Pos[d];
				dSol += nSol[d]*nSol[d];
      }
      dSol = std::sqrt( dSol );
      dMin = dSol<dMin ? dSol : dMin;
			
      // Check if close enough
      if( dMin<=2.0*parameter(0)*parameter(1) ) {
        // Add neighbor
				//std::cerr << i << " " << j << " added as neighbor pair.\n";
        neighbor[i].push_back(j);
        neighbor[j].push_back(i);
        numNeigh++;
      }
    }
  }
  // Copy the new neighbours into the cells 
	// (Caveat: Not a good solution)
  for(size_t i=0 ; i<N ; i++ )
    O.compartment(i).setNeighbor( neighbor[i] );
  
  // Set the previosu time variable to current t value
  setPreviousTime(t);
  
  return numNeigh;
}

unsigned int NeighborhoodDistanceCigar::
update(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{  
  if( parameter(2) == 1.0 && t-previousTime()>parameter(3) )
    return create(O,y,t);
  
  return numNeighbor();
}

NeighborhoodFromFileInitial::
NeighborhoodFromFileInitial(std::vector<double> &paraValue, 
			    const std::string &inFileValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=1 ) {
    std::cerr << "NeighborhoodFromFileInitial::"
	      << "NeighborhoodFromFileInitial() "
	      << "Uses one parameter areaFlag.\n";
    exit(0);
  }
  
  // Check for consistency in parameter values
  //if( paraValue[0]!=0.0 && paraValue[0]!=1.0 ) {
  // std::cerr << "NeighborhoodFromFileInitial::"
  //						<< "NeighborhoodFromFileInitial() "
  //						<< "parameter(0) (areaFlag) must be 0 or 1.\n";
  // exit(0);
  //}
  
  //Set the variable values
  setId("neighborhoodFromFileInitial");
  setParameter(paraValue);  
  setInFile(inFileValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "areaFlag";
  setParameterId( tmp );
}

unsigned int NeighborhoodFromFileInitial::
create(Organism &O, std::vector< std::vector<double> > &y,
       double t ) 
{
  int verbose=0;
  // Open the neighborhood file
  const char* fileName=inFile().c_str();  
  std::ifstream IN( fileName );
  if( !IN ) {
    std::cerr << "NeighborhoodFromFileInitial::create() - "
	      << "Cannot open file " << fileName << "\n\n\7";exit(-1);}
  
  if (verbose)
    std::cerr << "NeighborhoodFromFileInitial::create() Reading neighborhood file " 
	      << fileName << " ... ";
  
  // Create and resize temporary neighbor and area vectors
  unsigned int numNeigh=0;
  std::vector<size_t> neighbor;
  std::vector<double> neighborArea;
  std::vector< std::vector<double> > neighborVariable;
  size_t tmpInt,K,numNeighPara;
  IN >> tmpInt;//N
  if( tmpInt != O.numCompartment() ) {
    std::cerr << "NeighborhoodFromFileInitial::create() - "
	      << "Wrong number of compartments in file " << fileName 
	      << "\n\n\7";
    std::cerr << tmpInt << " " << O.numCompartment() << "\n";
    exit(-1);
  }
  
  // number of parameters (at least 1, if areas are to be read)
  IN >> numNeighPara;
  if( numNeighPara < static_cast<size_t>(parameter(0)) ) {
    std::cerr << "NeighborhoodFromFileInitial::create() - "
	      << "Wrong number of area variables in file " << fileName 
	      << "\n\n\7";
    exit(-1);
  }
  for( size_t i=0 ; i<O.numCompartment() ; i++ ) {
    IN >> tmpInt;//i
    if( tmpInt != i ) {
      std::cerr << "NeighborhoodFromFileInitial::create() - "
		<< "Strange compartment number in file " << fileName << std::endl
		<< "Read index " << tmpInt << " expected " << i << std::endl;
      exit(-1);
    }
    IN >> K;//K
    numNeigh +=K;
    neighbor.resize(K);
    //Read neighbor index
    for( size_t k=0 ; k<K ; k++ )
      IN >> neighbor[k];//neigh
    size_t paraCount=0;
    //Read neighbor areas
    if( parameter(0)!=0.0 ) {
      neighborArea.resize(K);
      for( size_t k=0 ; k<K ; k++ )
	IN >> neighborArea[k];//neighArea
      paraCount++;
    }
    
    neighborVariable.resize(numNeighPara-paraCount);
    for (size_t s = 0; s < neighborVariable.size(); ++s)
      neighborVariable[s].resize(K);
    //Read additional parameters if present (but do not use them)
    for( size_t pC=paraCount ; pC<numNeighPara ; pC++ ) {
      for( size_t k=0 ; k<K ; k++ ) {
	IN >> neighborVariable[pC-paraCount][k];//neigh
      }
    }
    
    //Set the values in the compartment
    O.compartment(i).setNeighbor( neighbor );    
    if( parameter(0)!=0.0 )
      O.compartment(i).setNeighborArea( neighborArea );    
    if( neighborVariable.size() )
      O.compartment(i).setNeighborVariable( neighborVariable );    
  }
  
  //Set the previous time variable to current t value
  setPreviousTime(t);
  if (verbose)
    std::cerr << "Done\n";
  IN.close();
  return numNeigh/2;
}

unsigned int NeighborhoodFromFileInitial::
update(Organism &O, std::vector< std::vector<double> > &y,
       double t ) {
  return numNeighbor();
}

NeighborhoodDistanceSpherePeriodicBox::
NeighborhoodDistanceSpherePeriodicBox(std::vector<double> &paraValue, 
																			std::vector< std::vector<size_t> > &indValue)
{
	// Do some checks on the parameters and variable indeces
  if (paraValue.size() < 3 || paraValue.size() > 6) {
    std::cerr << "NeighborhoodDistanceSpherePeriodicBox::"
							<< "NeighborhoodDistanceSpherePeriodicBox() "
							<< "Uses three to six parameters r_frac, upDateFlag, dt_min and "
							<< "the width/height/depth of the boundary box.\n";
    exit(EXIT_FAILURE);
  }
	if (indValue.size() > 0) {
    std::cerr << "NeighborhoodDistanceSpherePeriodicBox::"
							<< "NeighborhoodDistanceSpherePeriodicBox() "
							<< "No variable indices are used.\n";
    exit(EXIT_FAILURE);
  }

	// Check for consistency in parameter values
	if (paraValue[0] < 0.0) {
    std::cerr << "NeighborhoodDistanceSpherePeriodicBox::"
							<< "NeighborhoodDistanceSpherePeriodicBox() "
							<< "parameter(0) (r_frac) must be larger than zero.\n";
    exit(EXIT_FAILURE);
  }
	
  if (paraValue[1] != 0.0 && paraValue[1] != 1.0) {
    std::cerr << "NeighborhoodDistanceSpherePeriodicBox::"
							<< "NeighborhoodDistanceSpherePeriodicBox() "
							<< "parameter(1) (updateFlag) must be 0 or 1.\n";
    exit(EXIT_FAILURE);
  }
  
  // Set the variable values
  setId("neighborhoodDistanceSpherePeriodicBox");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp(numParameter());
  tmp.resize(numParameter());
  tmp[0] = "r_frac";
  tmp[1] = "updateFlag";
  tmp[2] = "dt_min";
  setParameterId( tmp );
}

unsigned int NeighborhoodDistanceSpherePeriodicBox::
create(Organism &O,
			 std::vector< std::vector<double> > &y, double t)
{
	size_t N = O.numCompartment();
  std::vector< std::vector<size_t> > neighbor(N);
  double d, r;
  unsigned int numNeigh = 0;
  
  size_t rIndex = O.numTopologyVariable() - 1;
	
	size_t dIndex = numParameter() - 3;
	
	std::vector<double> size(dIndex);
  for (size_t i = 3; i < numParameter(); ++i)
		size[i - 3] = parameter(i);
	
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = i + 1; j < N; ++j) {
			std::vector<int> k(dIndex + 1, -1);			
			do {
				r = y[i][rIndex] + y[j][rIndex];
				d = 0.0;
				for (size_t dim = 0; dim < rIndex; ++dim) {
					double tmp;
					if (dim < dIndex)
						tmp = (y[i][dim] - (size[dim] * k[dim] + y[j][dim]));
					else
						tmp = (y[i][dim] - y[j][dim]);
					d += (tmp * tmp);
				}
				d = std::sqrt(d);
				if (d <= r * parameter(0)) {
					neighbor[i].push_back(j);
					neighbor[j].push_back(i);
					numNeigh++;
				}						
				
				for (size_t l = 0; l < dIndex + 1; ++l) {
					++k[l];
					if (k[l] == 2)
						k[l] = -1;
					else
						break;
				}
			} while (k[dIndex] == -1);
		}
	}
	
  // Copy the new neighbours into the cells 
	// (Caveat: Not a good solution)
  for (size_t i = 0; i < N; ++i)
    O.compartment(i).setNeighbor(neighbor[i]);

  //Set the previous time variable to current t value
  setPreviousTime(t);
  setNumNeighbor(numNeigh);
	
  return numNeigh;
}

unsigned int NeighborhoodDistanceSpherePeriodicBox::
update(Organism &O,
			 std::vector< std::vector<double> > &y, double t)
{
  if (parameter(1) == 1.0 && t-previousTime() > parameter(2))
    return create(O, y, t);
	
  return numNeighbor();
}


NeighborhoodIndex::
NeighborhoodIndex(std::vector<double> &paraValue, 
									std::vector< std::vector<size_t> > &indValue)
{
	// Do some checks on the parameters and variable indeces
  if (paraValue.size()!=3) {
    std::cerr << "NeighborhoodIndex::"
							<< "NeighborhoodIndex() "
							<< "Uses three parameters periodicFlag, updateFlag and deltaTime." << std::endl;
    exit(EXIT_FAILURE);
  }
	if (indValue.size() > 0) {
    std::cerr << "NeighborhoodIndex::"
							<< "NeighborhoodIndex() "
							<< "No variable indices are used.\n";
    exit(EXIT_FAILURE);
  }

	// Check for consistency in parameter values
	if (paraValue[0]!=0.0 && paraValue[0]!=1.0) {
    std::cerr << "NeighborhoodIndex::"
							<< "NeighborhoodIndex() "
							<< "parameter(0) (periodicFlag) must be one (periodic) or zero (not).\n";
    exit(EXIT_FAILURE);
  }
	if (paraValue[1]!=0.0 && paraValue[1]!=1.0) {
    std::cerr << "NeighborhoodIndex::"
							<< "NeighborhoodIndex() "
							<< "parameter(1) (updateFlag) must be one (update) or zero (only initiate).\n";
    exit(EXIT_FAILURE);
  }
	if (paraValue[2]<0.0) {
    std::cerr << "NeighborhoodIndex::"
							<< "NeighborhoodIndex() "
							<< "parameter(2) (deltaTime) must be positive." << std::endl;
    exit(EXIT_FAILURE);
  }
	
  // Set the variable values
  setId("neighborhoodDistanceSpherePeriodicBox");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp(numParameter());
  tmp.resize(numParameter());
  tmp[0] = "periodicFlag";
  tmp[1] = "updateFlag";
  tmp[2] = "deltaTime";
  setParameterId( tmp );
}

unsigned int NeighborhoodIndex::
create(Organism &O,
			 std::vector< std::vector<double> > &y, double t)
{
	size_t N = O.numCompartment();
  std::vector< std::vector<size_t> > neighbor(N);
  unsigned int numNeigh = 0;
  
  for (size_t i = 0; i < N; ++i) {
		if (i) {
			neighbor[i].push_back(i-1);
			++numNeigh;
		}
		if (i<N-1) { 
			neighbor[i].push_back(i+1);
			++numNeigh;
		}
	}						

  // Copy the new neighbours into the cells 
	// (Caveat: Not a good solution)
  for (size_t i = 0; i < N; ++i)
    O.compartment(i).setNeighbor(neighbor[i]);

  //Set the previous time variable to current t value
  setPreviousTime(t);
  setNumNeighbor(numNeigh);
	
  return numNeigh;
}

unsigned int NeighborhoodIndex::
update(Organism &O,
			 std::vector< std::vector<double> > &y, double t)
{
  if (parameter(1) == 1.0 && t-previousTime() > parameter(2))
    return create(O, y, t);
	
  return numNeighbor();
}

NullNeighborhood::NullNeighborhood(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue)
{
	setId("nullNeighborhood");
	setParameter(paraValue);
	setVariableIndex(indValue);

	std::vector<std::string> parameterIds(0);
  	setParameterId(parameterIds);
}

unsigned int NullNeighborhood::create(Organism &O, std::vector< std::vector<double> > &y, double t)
{
	return 0;
}

unsigned int NullNeighborhood::update(Organism &O, std::vector< std::vector<double> > &y, double t)
{
	return 0;
}
