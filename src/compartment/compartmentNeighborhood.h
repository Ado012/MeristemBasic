#ifndef COMPARTMENTNEIGHBORHOOD_H
#define COMPARTMENTNEIGHBORHOOD_H
//
// Filename     : compartmentNeighborhood.h
// Description  : Classes describing compartment neighborhoods
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2005
// Revision     : $Id: compartmentNeighborhood.h 638 2015-07-08 07:30:02Z henrik $
//

#include "baseCompartmentNeighborhood.h"
#include "../organism.h"

///
/// @brief Defines neighborhood by placing each cigar-shaped cell in a lattice
/// box
///
class NeighborhoodLatticeCigarMultipleWall: public BaseCompartmentNeighborhood {
  
 public:
  	
	///
	/// @brief Main constructor
	///
	/// This is the main constructor which checks and sets the parameters and
	/// variable indices that defines the neighborhood rule.
	///
	/// @param paraValue Vector with parameters used.
	///
	/// @param indValue Vector of vectors with variable indices used in the
	/// neighborhood calculations.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
	NeighborhoodLatticeCigarMultipleWall(std::vector<double> &paraValue, 
																			 std::vector< std::vector<size_t> > 
																			 &indValue );
	
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
	unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
											double t=0.0 );

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
	unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
											double t=0.0 );  
 private:
	std::vector< std::vector< size_t > > wallNeighborLattice;
};

//!Defines neighborhood by placing each cell in a lattice cell
class NeighborhoodLatticeCigar: public BaseCompartmentNeighborhood {
  
 public:
  
	///
	/// @brief Main constructor
	///
	/// This is the main constructor which checks and sets the parameters and
	/// variable indices that defines the neighborhood rule.
	///
	/// @param paraValue Vector with parameters used.
	///
	/// @param indValue Vector of vectors with variable indices used in the
	/// neighborhood calculations.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
  NeighborhoodLatticeCigar(std::vector<double> &paraValue, 
													 std::vector< std::vector<size_t> > 
													 &indValue );
  
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
  unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
		      double t=0.0 );

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
  unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
		      double t=0.0 );  

  bool check( size_t i, size_t j, size_t dimension,  
	      std::vector< std::vector<double> > &y);
	///
	/// @brief Checks whether two bact cells are overlapping
	///
  bool check2( std::vector<double> &yi, std::vector<double> &yj );
};

namespace Sphere {
  ///
  /// @brief Defines neighborhood from distances between spheres
  ///
  /// Simple @f$O(N^2)@f$ function for finding neighboring
  /// compartments. It uses three parameters and assumes that positions
  /// are set in topology as well as the radius (which should be the final
  /// topology variable). rFactor
  /// gives a factor of the sum of the two cell radii for which a threshold is
  /// set. If the distance is larger than the threshold, the cells are not
  /// considered as neighbors. updateFlag determines if the neighborhood
  /// is to be updated (1) or not (0). t_min sets a minimal time between updates.
  ///
  /// In the model file, the neighborhood is defined as:
  /// @verbatim
  /// Sphere::NeighborhoodDistance 3 0
  /// rFactor updateFlag t_min
  /// @endverbatim
  ///
  class NeighborhoodDistance : public BaseCompartmentNeighborhood {
    
  public:
    
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the neighborhood rule.
    ///
    /// @param paraValue Vector with parameters used.
    ///
    /// @param indValue Vector of vectors with variable indices used in the
    /// neighborhood calculations.
    ///
    /// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
    /// std::vector<size_t> > &, const std::string &)
    ///
    NeighborhoodDistance(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > 
			 &indValue );
    
    ///
    /// @brief Creates the (potential) neighborhood before simulation
    ///
    /// @see BaseCompartmentNeighborhood::create()
    ///
    unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
			double t=0.0 );
    
    ///
    /// @brief Updates the (potential) neighborhood during a simulation
    ///
    /// @see BaseCompartmentNeighborhood::update()
    ///
    unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
			double t=0.0 );  
  };
} //namespace Sphere

///
/// @brief Defines periodic neighborhood by looking at distances between
/// spheres
///
/// Simple O(N^2) function for finding neighboring compartments. rFactor gives
/// a factor of the sum of the two cell radiuses for which a threshold is
/// set. If the distance is larger than the threshold, the cells are not
/// considered as neighbors. It returns the total number of neighbors
/// found. Then it connects 'borders' in the spatial directions given as
/// variable index.
///
class NeighborhoodDistanceSpherePeriodic : 
	public BaseCompartmentNeighborhood {
  
 public:
  
	///
	/// @brief Main constructor
	///
	/// This is the main constructor which checks and sets the parameters and
	/// variable indices that defines the neighborhood rule.
	///
	/// @param paraValue Vector with parameters used.
	///
	/// @param indValue Vector of vectors with variable indices used in the
	/// neighborhood calculations.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
  NeighborhoodDistanceSpherePeriodic(std::vector<double> &paraValue, 
																		 std::vector< std::vector<size_t> > 
																		 &indValue );
  
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
  unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
											double t=0.0 );

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
  unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
											double t=0.0 );  
	///
	/// @brief Makes the system periodic by adding neighbors between 'boundary'
	/// cells
	///
  unsigned int addPeriodicBoundaries(std::vector< 
																		 std::vector<double> > &y, 
																		 std::vector< std::vector<size_t> > &neighbor,
																		 size_t dIndex,size_t rIndex);
  
};

///
/// @brief Defines neighborhood for spheres and adds walls in between
///
/// Simple O(N^2) function for finding neighboring compartments. rFactor gives
/// a factor of the sum of the two cell radiuses for which a threshold is
/// set. If the distance is larger than the threshold, the cells are not
/// considered as neighbors. Then walls in between each cell-neighbor pair are
/// created, new neighbors added and the cell-cell neighborhood is moved into
/// level two of the neighborhood matrix. 
///
class NeighborhoodDistanceSphereWithWalls : 
	public BaseCompartmentNeighborhood {
  
 public:
  
	///
	/// @brief Main constructor
	///
	/// This is the main constructor which checks and sets the parameters and
	/// variable indices that defines the neighborhood rule.
	///
	/// @param paraValue Vector with parameters used.
	///
	/// @param indValue Vector of vectors with variable indices used in the
	/// neighborhood calculations.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
  NeighborhoodDistanceSphereWithWalls(std::vector<double> &paraValue, 
																			std::vector< std::vector<size_t> > 
																			&indValue );
  
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
  unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
											double t=0.0 );

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
  unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
											double t=0.0 );  
};

///
/// @brief Defines neighborhood by looking at distances between sphere-buds
///
/// Simple O(N^2) function for finding neighboring compartments. rFactor gives
/// a factor of the sum of the two cell radiuses for which a threshold is
/// set. If the distance is larger than the threshold, the cells are not
/// considered as neighbors.
///
class NeighborhoodDistanceSphereBud : public BaseCompartmentNeighborhood {
  
 public:
  
	///
	/// @brief Main constructor
	///
	/// This is the main constructor which checks and sets the parameters and
	/// variable indices that defines the neighborhood rule.
	///
	/// @param paraValue Vector with parameters used.
	///
	/// @param indValue Vector of vectors with variable indices used in the
	/// neighborhood calculations.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
  NeighborhoodDistanceSphereBud(std::vector<double> &paraValue, 
																std::vector< std::vector<size_t> > 
																&indValue );
  
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
  unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
	     double t=0.0 );

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
  unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
	     double t=0.0 );  
};

///
/// @brief Defines neighborhood by looking at distances between ellipses
///
/// Simple O(N^2) function for finding potential neighboring elliptical
/// compartments.  rFactor gives a factor of the sum of the two cell radii for
/// which a threshold is set. If the distance is larger than the threshold,
/// the cells are not considered as potential neighbors. In this case of
/// elliptical compartments, the radius is replaced with a, the distance from
/// the center to the ellipse following the major axis. It returns the total
/// number of neighbors found.
///
class NeighborhoodDistanceEllipse : public BaseCompartmentNeighborhood {
  
 public:
  
	///
	/// @brief Main constructor
	///
	/// This is the main constructor which checks and sets the parameters and
	/// variable indices that defines the neighborhood rule.
	///
	/// @param paraValue Vector with parameters used.
	///
	/// @param indValue Vector of vectors with variable indices used in the
	/// neighborhood calculations.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
  NeighborhoodDistanceEllipse(std::vector<double> &paraValue, 
															std::vector< std::vector<size_t> > 
															&indValue );
  
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
  unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
	     double t=0.0 );

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
  unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
	     double t=0.0 );  
};

///
/// @brief Defines neighborhood by looking at distances between cigars
///
/// Simple O(N^2) function for finding potential neighboring cigar
/// compartments.  rFactor gives a factor of the sum of the two cell radii for
/// which a threshold is set. If the distance is larger than the threshold,
/// the cells are not considered as potential neighbors. In this case of
/// elliptical compartments, the radius is replaced with a, the distance from
/// the center to the ellipse following the major axis. It returns the total
/// number of neighbors found.
///
class NeighborhoodDistanceCigar : public BaseCompartmentNeighborhood {
  
 public:
  
	///
	/// @brief Main constructor
	///
	/// This is the main constructor which checks and sets the parameters and
	/// variable indices that defines the neighborhood rule.
	///
	/// @param paraValue Vector with parameters used.
	///
	/// @param indValue Vector of vectors with variable indices used in the
	/// neighborhood calculations.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
  NeighborhoodDistanceCigar(std::vector<double> &paraValue, 
														std::vector< std::vector<size_t> > 
														&indValue );
  
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
  unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
	     double t=0.0 );

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
  unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
	     double t=0.0 );  
};

///
/// @brief Defines initial neighborhood from file with no update. 
///
/// This class reads a neighborhood from file at the initial phase and does
/// not update the neighborhood after that. It is a special kind of
/// compartmentNeighborhood since it uses one parameter and a file name for
/// initiation. In the model file it should be defined as:
///
/// @verbatim
/// neighborhoodFromFileInitial 1
/// areaFlag
/// fileName
/// @endverbatim
/// 
/// where the parameter give is a flag setting whether an area should be read
/// (areaFlag=1) or not (areaFlag=0).
///
/// The neighborhood file should have the format:
///
/// @verbatim
/// N_i N_A
/// i K_i n_1 ... n_K [A_1 ... A_K]
/// ...
/// @endverbatim
///
/// where @f$N_i@f$ is the number of compartments (cells), which is checked
/// against the number compartments in the model (Organism). @f$N_A@f$ is the
/// number of additional variables (at least 1 if @f$areaFlag@f$ is set to 1,
/// see above). Even if @f$N_A>1@f$ only the neighbors and neighbor areas will
/// be read. The cell index @f$i@f$ should be increasing from 0, @f$K_i@f$ is
/// the number of neighbors, and the neighbor indices are listed in
/// @f$n_k@f$. If @f$A>0@f$ and @f$areaFlag=1@f$, also the neighbor areas will
/// be read from the list @f$A_k@f$.
///
class NeighborhoodFromFileInitial : public BaseCompartmentNeighborhood {
  
private:
  
  std::string inFile_;
  
public:
  
	///
	/// @brief Main constructor, which differ from all other since it uses a
	/// file (inFileValue) to construct a neighborhood.
	///
	/// For more information see code, espesially
	/// BaseCompartmentNeighborhood::createCompartmentNeighborhood().
	///
  NeighborhoodFromFileInitial(std::vector<double> &paraValue,
															const std::string &inFileValue);
  
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
  unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
											double t=0.0 );

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// This does not do any update, but the neighborhood read in create is
	/// kept.
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
  unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
											double t=0.0 );

	///
	/// @brief Returns the neighborhood file name
	///
  inline std::string inFile() const;

	///
	/// @brief Sets the neighborhood file name
	///
  inline void setInFile(const std::string &inFileValue);
};

inline std::string NeighborhoodFromFileInitial::inFile() const {
  return inFile_;
}

inline void NeighborhoodFromFileInitial::
setInFile(const std::string &value) {
  inFile_=value;
}

class NeighborhoodDistanceSpherePeriodicBox 
	: public BaseCompartmentNeighborhood
{

 public:

	///
	/// @brief Main constructor
	///
	/// This is the main constructor which checks and sets the parameters and
	/// variable indices that defines the neighborhood rule.
	///
	/// @param paraValue Vector with parameters used.
	///
	/// @param indValue Vector of vectors with variable indices used in the
	/// neighborhood calculations.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
	NeighborhoodDistanceSpherePeriodicBox(std::vector<double> &paraValue, 
																				std::vector< std::vector<size_t> > 
																				&indValue);
	
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
  unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
											double t = 0.0);

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
  unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
											double t = 0.0);
};

///
/// @brief Sets neighborhood from indices (i-1,i+1 neighbor to i).
///
/// This class defines the neighborhood by setting adjacent indices as
/// neighbors, i. e. i-1 and i+1 are neighbors to i. If the periodic flag is
/// set it also regards i=0 and i=N-1, where N is the number of compartments
/// as neighbors.
///
/// @verbatim
/// neighborhoodIndex 1 0
/// periodicFlag
/// @endverbatim
/// 
/// where the parameter is a flag setting whether the neighbors are periodic
/// (periodicFlag=1) or not (periodicFlag=0).
///
class NeighborhoodIndex 
	: public BaseCompartmentNeighborhood
{

 public:

	///
	/// @brief Main constructor
	///
	/// This is the main constructor which checks and sets the parameters and
	/// variable indices that defines the neighborhood rule.
	///
	/// @param paraValue Vector with parameters used.
	///
	/// @param indValue Vector of vectors with variable indices used in the
	/// neighborhood calculations.
	///
	/// @see createCompartmentNeighborhood(std::vector<double> &, std::vector<
	/// std::vector<size_t> > &, const std::string &)
	///
	NeighborhoodIndex(std::vector<double> &paraValue, 
																				std::vector< std::vector<size_t> > 
																				&indValue);
	
	///
	/// @brief Creates the (potential) neighborhood before simulation
  ///
	/// @see BaseCompartmentNeighborhood::create()
	///
  unsigned int create(Organism &O, std::vector< std::vector<double> > &y,
											double t = 0.0);

	///
	/// @brief Updates the (potential) neighborhood during a simulation
  ///
	/// @see BaseCompartmentNeighborhood::update()
	///
  unsigned int update(Organism &O, std::vector< std::vector<double> > &y,
											double t = 0.0);
};

/// 
/// @brief Neighborhood that is doing nothing and has no parameters.
///
class NullNeighborhood : public BaseCompartmentNeighborhood
{
public:
	NullNeighborhood(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
	unsigned int create(Organism &O, std::vector< std::vector<double> > &y, double t = 0.0);
	unsigned int update(Organism &O, std::vector< std::vector<double> > &y, double t = 0.0);
};

#endif


