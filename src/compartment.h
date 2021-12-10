/**
 * Filename     : compartment.h
 * Description  : A class describing a compartment with neighborhood
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : October 2003
 * Revision     : $Id: compartment.h 637 2015-07-07 12:55:36Z henrik $
 */
#ifndef COMPARTMENT_H
#define COMPARTMENT_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<cassert>
#include <cstdlib>
#include "nucleus.h" //ADDITION AD 032117
class Organism;
class Topology;

//define probability matrix to hold possible events
class probabilitySegment {
public:      //putting begin and end as int messes it up for some reason?
double begin,end;
int site;
std::string action;
int chromosome; //tracks which crm is being affected
};

///
/// @brief A class describing properties of a compartment sa embedding and topology
///
/// The Compartment handles the update and information of the embedding
/// and topology. It holds a vector of variable values for all dynamical
/// variables in the compartment. It holds the neighbor information, ie
/// the neighbor indeces and if applicable the cross section areas. It
/// can also store a vector of parameters that can be updated during a
/// simulation, and used in baseReaction calculations of the
/// derivatives.
///
class Compartment {
  
 private:
  
  size_t index_;
  std::string id_;           
  std::vector<double> variable_;
  std::vector< std::vector<size_t> > neighbor_;
  std::vector<double> neighborArea_;
  std::vector< std::vector<double> > neighborVariable_;
  double volume_;
  std::vector<double> parameter_;
  const Topology* topology_;
  const Organism* organism_;
	
public:




Nucleus cellNucleus; //ADDITION AD 032117
//current activation scheme can have both crms considered together, may have to rewrite other activation systems
int geneCRM[2][5] = {{0,0,0,0,0}, {0,0,0,0,0}};//ADDITION AD 04SOMETHING17, revised 042519
int geneCRMMarker[2][5] = {{0,0,0,0,0}, {0,0,0,0,0}};//ADDITION AD 061619
//std::array needs extra braces for some bizarre reason
//RNA polytimer and polymeraseBindSite are gene gene marker marker in arrangment
int polymeraseBindSite[2][2]={{0,0}, {0,0}}; //ADDITION AD 111819
double crmTimer[2][5] = {{0,0,0,0,0}, {0,0,0,0,0}};//ADDITION AD 042519
double crmTimer2[2][5] = {{0,0,0,0,0}, {0,0,0,0,0}};//ADDITION AD 061619
double unbindTimer[4][5] = {{0,0,0,0,0}, {0,0,0,0,0},{0,0,0,0,0}, {0,0,0,0,0}};//ADDITION AD 010421
double rnaPolTimer[2][2]={{0,0}, {0,0}}; //ADDITION AD 111819

int timer=0;//ADDITION AD 110618
double stochasticStepOverflow[2] = {0,0}; //Addition 042419 Necessary to measure stochastic steps overflowing to next ODE step. Edit: 061219 made in compatible with crm/marker switch
double probabilityDeltaSum[2] = {0,0};//Addition 061219 measures the total probability of all events at a stochastic timeStep increment, needed to calculate timestep length and event. Edit 61219: made it compatible with crm/marker switch
double previousTimeStep[2]= {0,0};
double latestTimeStep[2]= {0,0};
double previousRandTimeStep[2]= {0,0};
double latestRandTimeStep[2]= {0,0};
double monFireCountDown[2][5] = {{0,0,0,0,0}, {0,0,0,0,0}};//ADDITION AD 042519
double monFireCountDown2[2][5] = {{0,0,0,0,0}, {0,0,0,0,0}};//ADDITION AD 061619
double clv3StepContribution=0;
double wusP_FromNuc=0;
double wusP_FromCyto=0;
double newAssCk=0;
double newDissCk=0;
//tracking variables
int monTracker=0;
int dimTracker=0;
int polFire=0;
int averageAge=0;
//tracking variables

int eventNum[2]={0,0};
probabilitySegment probabilityMatrix[60];//if marker is used again may need another matrix
int firstbind=0;
  Compartment();
  Compartment( const Compartment & compartmentCopy );//Change to pass a compartment_nucleated?
  Compartment( const Organism* organismVal, 
	       const Topology* topologyVal, size_t numSpecies, 
	       size_t index=static_cast<size_t>(-1), 
	       const std::string &idValue="");
  Compartment( const Organism* organismVal,
	       const Topology* topologyVal, 
	       const std::vector<double> &variableVal,
	       size_t index=static_cast<size_t>(-1), 
	       const std::string &idValue="" );
  
  ~Compartment();
  
  //Compartment & operator=( const Compartment & compartmentCopy );
  
	///
	/// @brief Returns the id string of the compartment.
	///
  inline std::string id() const;
	///
	/// @brief Returns the index of the compartment, which should be the same as
	/// the index of the stored compartment list in the Organism class.
	///
  inline size_t index() const;
	///
	/// @brief Returns the compartment topology as a pointer.
	///
	inline const Topology* topology() const;
	///
	/// @brief Returns the Organism which the compartment is part of.
	///
	inline const Organism* organism() const;
	///
	/// @brief Returns the number of variables stored by the compartment.
	///
	inline size_t numVariable() const;  
	///
	/// @brief Returns the variable value for a specified (by index) variable.
	///
  inline double variable(size_t i) const;
	///
	/// @brief Returns the complete variable vector stored by the compartment.
	///
  inline const std::vector<double> variable() const;
	///
	/// @brief Returns the number of neighbors for the compartment (stored in the first level).
	///
  inline size_t numNeighbor() const;
	///
	/// @brief Returns the index of the specified neighbor.
	///
  inline size_t neighbor(size_t i) const;
	///
	/// @brief Returns the neighbors as a vector of indices.
	///
  inline const std::vector<size_t> neighbor() const;
	///
	/// @brief Returns a reference to the neighbor compartment.
	///
	const Compartment & neighborRef(size_t i);
	///
	/// @brief Returns the index of a by index specified neighbor compartment,
	/// and -1 if not a neighbor.
	///
	size_t neighborIndex(size_t compartmentIndex) const;
	///
	/// @brief Returns the number of different neighbor types stored by the compartment.
	///
	/// Storing several 'levels' of neighbors allows for example to store
	/// nearest neighbors as well as next-nearest neighbors, e g wall neighbors
	/// and cell neighbors.
	///
  inline size_t numNeighborLevel() const;
	///
	/// @brief Returns the number of neighbors at a specified level.
	///
	/// numNeighborAtLevel(0) is identical to numNeighbor().
	///
	/// @see Compartment::numNeighborLevel()
	/// @see Compartment::numNeighbor()
	///
  inline size_t numNeighborAtLevel(size_t l) const;
	///
	/// @brief Returns a neighbor index at a specified level.
	///
	/// neighborAtLevel(0,i) is identical to neighbor(i).
	///
	/// @see Compartment::numNeighborLevel()
	/// @see Compartment::neighbor(size_t)
	///
  inline size_t neighborAtLevel(size_t l,size_t i) const;
	///
	/// @brief Returns the neighbor vector at a specified level.
	///
	/// neighborAtLevel(0) is identical to neighbor().
	///
	/// @see Compartment::numNeighborLevel()
	/// @see Compartment::neighbor()
	///
  inline const std::vector<size_t> neighborAtLevel(size_t l) const;
	///
	/// @brief Returns the number of neighborAreas stored by the compartment.
	///
	inline size_t numNeighborArea() const;
	///
	/// @brief Returns the area between this and the neighbor.
	///
  inline double neighborArea(size_t i) const;
	///
	/// @brief Returns a vector of neighbor areas.
	///
  inline const std::vector<double> neighborArea() const;
	///
	/// @brief Returns the number of variables connected to neighbors.
	///
	/// A compartment can store variables connected to the intersection between
	/// itself and its neighbors. This could for example be transport mediators
	/// sitting in the membrane between a cell and a wall.
	///
	inline size_t numNeighborVariable() const;
	///
	/// @brief returns a variable value connected to the intersection between
	/// the compartment and its neighbor.
	///
	/// @see Compartment::numNeighborVariable()
	///
	inline double neighborVariable(size_t index, size_t neigh) const;
	///
	/// @brief Returns the volume of the compartment.
	///
	/// @note This volume is not always up to date and should not be relied upon
	/// unless you know what you are doing.
	///
	inline double volume() const;
	///
	/// @brief Calculates and returns the volume of the compartment.
	///
	/// @note This volume is not always up to date and should not be relied upon
	/// unless you know what you are doing (uses Compartment.variable).
	///
	double getVolume();
	///
	/// @brief Calculates and returns the volume of the compartment.
	///
	/// The current variable values (y) are supplied such that this will be 
	/// an up to date calculation
	///
	double getVolume(std::vector<double> &y);
	///
	/// @brief Returns the start index of the topology in the variable vector.
	///
	/// In the current implementation the topology variables are stored at the
	/// beginning of the vector, and hence this function always returns 0.
	///
  inline size_t topologyStart() const;
	///
	/// @brief Returns the number of topology variables.
	///
  size_t numTopologyVariable() const;
	///
	/// @brief Returns the start index of the species (concentrations) in the
	/// variable vector.
	///
	/// In the current implementation the topology variables are stored at the
	/// beginning of the vector, and the species follows after. Hence this
	/// function always returns the number of topology variables.
	///
  inline size_t speciesStart() const;
	///
	/// @brief Returns the number of species (concentrations) stored by the
	/// compartment.
	///
	/// This is calculated as number of variables minus the number of topology
	/// variables stored by the compartment.
	///
  inline size_t numSpecies() const;
	///
	/// @brief Returns the number of dimensions from the topology.
	///
	size_t numDimension() const;
	///
	/// @brief Returns the number of parameters stored in the compartment.
	///
	/// Right now, this is obselete since this was used to store some special
	/// variables in a flux model no longer in use. This should be used in a
	/// future to store variables that are not updated, e.g. cell/wall marker.
	///
  inline size_t numParameter() const;  
	///
	/// @brief Returns a parameter value stored by the compartment.
	///
	/// @see Compartment::numParameter()
	///
  inline double parameter(size_t i) const;
	///
	/// @brief Returns the parameter value vector stored by the compartment.
	///
	/// @see Compartment::numParameter()
	///
  inline const std::vector<double> parameter() const;
	///
	/// @brief Sets the id string.
	///
  inline void setId(const std::string &value);
	///
	/// @brief Sets the compartment index.
	///
  inline void setIndex(size_t value);
	///
	/// @brief Sets the compartment topology. 
	///
  void setTopology(const Topology* value);
	///
	/// @brief Sets the compartment organism.
	///
  void setOrganism(const Organism* value);
	///
	/// @brief Sets a variable value.
	///
  inline void setVariable(size_t i,double value);
	///
	/// @brief Sets the vector of variables. 
	///
  inline void setVariable(const std::vector<double> &value);
	///
	/// @brief Sets a neighbor index.
	///
  inline void setNeighbor(size_t i,size_t value);
	///
	/// @brief Sets a neighbor vector.
	///
  inline void setNeighbor(const std::vector<size_t> &value);
	///
	/// @brief Adds a neighbor to the compartment.
	///
  inline void addNeighbor(size_t value);
	///
	/// @brief Removes a neighbor from the compartment.
	///
	/// This function checks if the supplied index of a compartment is a
	/// neighbor, and removes it if this is the case.
	///
	/// @return Returns the number of removed neighbors (one if success zero otherwise).
	///
  unsigned int removeNeighbor(size_t index);
	///
	/// @brief Changes the index of a neighbor.
	///
  unsigned int changeNeighborIndex(size_t oldIndex,size_t newIndex);
	///
	/// @brief Sets a neighbor at specified level and position.
	///
  inline void setNeighborAtLevel(size_t l,size_t i,size_t value);
	///
	/// @brief Sets a complete neighbor vector at specified level.
	///
  inline void setNeighborAtLevel(size_t l,const std::vector<size_t> &value);
	///
	/// @brief Adds a neighbor at specified level.
	///
  inline void addNeighborAtLevel(size_t l,size_t value);
	///
	/// @brief Removes a neighbor at specified level.
	///
  unsigned int removeNeighborAtLevel(size_t l,size_t index);
	///
	/// @brief Changes the neighbor at specified level.
	///
  unsigned int changeNeighborIndexAtLevel(size_t l,size_t oldIndex,size_t newIndex);
	///
	/// @brief Sets a neighbor area.
	///
  inline void setNeighborArea(size_t i,double value);
	///
	/// @brief Sets the vector with neighbor areas.
	///
  inline void setNeighborArea(const std::vector<double> &value);
	///
	/// @brief Adds a neighborarea.
	///
  inline void addNeighborArea(double value);
	///
	/// @brief Sets a complete set of neighbor variables.
	///
	inline void setNeighborVariable( size_t varIndex, size_t neighIndex, double value);
	///
	/// @brief Sets a neighbor variable value.
	///
	inline void setNeighborVariable( const std::vector< std::vector<double> > &value);
	///
	/// @brief Sets the compartment volume.
	///
  inline void setVolume(double value);
	///
	/// @brief Sets a compartment parameter.
	///
  inline void setParameter(size_t i,double value);
	///
	/// @brief Sets a complete set of compartment parameters.
	///
  inline void setParameter(const std::vector<double> &value);
};

inline std::string Compartment::id() const 
{
  return id_;
}

inline size_t Compartment::index() const 
{
  return index_;
}

inline size_t Compartment::numVariable() const 
{
  return variable_.size();
}

inline size_t Compartment::numNeighbor() const 
{
  return neighbor_.size()>0 ? neighbor_[0].size() : 0;
}

inline size_t Compartment::numNeighborLevel() const 
{
  return neighbor_.size();
}

inline size_t Compartment::numNeighborAtLevel(size_t l) const 
{
  return neighbor_.size()>l ? neighbor_[l].size() : 0;
}

inline size_t Compartment::numNeighborArea() const 
{
  return neighborArea_.size();
}

inline size_t Compartment::numNeighborVariable() const 
{
  return neighborVariable_.size();
}

inline double Compartment::neighborVariable(size_t index, size_t neigh) const 
{
	assert(index<neighborVariable_.size());
	assert(neigh<neighborVariable_[index].size());	
	return neighborVariable_[index][neigh];
}

inline double Compartment::volume() const 
{
  return volume_;
}

inline const Topology* Compartment::topology() const 
{
  return topology_;
}

inline const Organism* Compartment::organism() const 
{
  return organism_;
}

inline double Compartment::variable(size_t i) const 
{
  return variable_[i];
}

inline const std::vector<double> Compartment::variable() const
{
  return variable_;
}

inline size_t Compartment::neighbor(size_t i) const 
{
  if( neighbor_.size()<1 ) {
    std::cerr << "Comparment::neighbor(i) No neighs defined.\n";
    exit(-1);
  }
  return neighbor_[0][i];
}

inline const std::vector<size_t> Compartment::neighbor() const
{
  if( neighbor_.size()<1 ) {
    std::cerr << "vector<int> Comparment::neighbor(i) "
							<< "No neighs defined.\n";
    exit(-1);
  }
  return neighbor_[0];
}

inline size_t Compartment::neighborAtLevel(size_t l,size_t i) const 
{
  if( neighbor_.size()<l ) {
    std::cerr << "int Comparment::neighbor(l,i) "
							<< "No neighs defined at level " << l << ".\n";
    exit(-1);
  }
  return neighbor_[l][i];
}

inline const std::vector<size_t> Compartment::neighborAtLevel(size_t l) const
{
  if( neighbor_.size()<l ) {
    std::cerr << "vector<int> Comparment::neighbor(l,i) "
							<< "No neighs defined at level " << l << ".\n";
    exit(-1);
  }
  return neighbor_[l];
}

inline double Compartment::neighborArea(size_t i) const 
{
	return neighborArea_[i];
}

inline const std::vector<double> Compartment::neighborArea() const
{
  return neighborArea_;
}

inline size_t Compartment::topologyStart() const 
{
  return 0;
}

inline size_t Compartment::speciesStart() const 
{
  return numTopologyVariable();
}

inline size_t Compartment::numSpecies() const 
{
  return numVariable()-numTopologyVariable();
}

inline size_t Compartment::numParameter() const 
{
  return parameter_.size();
}

inline double Compartment::parameter(size_t i) const 
{
  return parameter_[i];
}

inline const std::vector<double> Compartment::parameter() const
{
  return parameter_;
}

inline void Compartment::setId(const std::string &value) 
{
  id_ = value;
}

inline void Compartment::setIndex( size_t value ) 
{
  index_ = value;
}

inline void Compartment::setVariable( const std::vector<double> &value) 
{
  variable_ = value;
}

inline void Compartment::setVariable(size_t i, double value) 
{
  variable_[i] = value;
}

inline void Compartment::setNeighbor( const std::vector<size_t> &value) 
{
  if( neighbor_.size()<1 )
    neighbor_.resize(1);
  neighbor_[0] = value;
}

inline void Compartment::setNeighbor(size_t i, size_t value) 
{
  if( neighbor_.size()<1 )
    neighbor_.resize(1);
  if( neighbor_[0].size()<=i )
    neighbor_[0].resize(i+1);
  neighbor_[0][i] = value;
}

inline void Compartment::addNeighbor(size_t value) 
{
  if( neighbor_.size()<1 )
    neighbor_.resize(1);
  neighbor_[0].push_back(value);
}

inline void Compartment::
setNeighborAtLevel(size_t l, const std::vector<size_t> &value) 
{
  if( neighbor_.size()<=l )
    neighbor_.resize(l+1);
  neighbor_[l] = value;
}

inline void Compartment::
setNeighborAtLevel(size_t l,size_t i, size_t value) 
{
  if( neighbor_.size()<=l )
    neighbor_.resize(l+1);
  if( neighbor_[0].size()<=i )
    neighbor_[0].resize(i+1);
  neighbor_[l][i] = value;
}

inline void Compartment::addNeighborAtLevel(size_t l,size_t value) 
{
  if( neighbor_.size()<=l )
    neighbor_.resize(l+1);
  neighbor_[l].push_back(value);
}

inline void Compartment::setNeighborArea( const std::vector<double> &value) 
{
  neighborArea_ = value;
}

inline void Compartment::setNeighborArea(size_t i, double value) 
{
  neighborArea_[i] = value;
}

inline void Compartment::addNeighborArea(double value) 
{
  neighborArea_.push_back(value);
}

inline void Compartment::setNeighborVariable( size_t varIndex, size_t neighIndex, double value)
{
	neighborVariable_[varIndex][neighIndex] = value;
}

inline void Compartment::setNeighborVariable( const std::vector< std::vector<double> > &value) 
{
	neighborVariable_ = value;
}

inline void Compartment::setVolume(double value) 
{
  volume_ = value;
}

inline void 
Compartment::setParameter( const std::vector<double> &value) {
  parameter_ = value;
}

inline void Compartment::setParameter(size_t i, double value) 
{
  parameter_[i] = value;
}

#endif
