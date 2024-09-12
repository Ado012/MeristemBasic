//
// Filename     : organism.h
// Description  : A class describing an organism 
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2003
// Revision     : $Id: organism.h 646 2015-08-23 18:35:06Z henrik $
//
#ifndef ORGANISM_H
#define ORGANISM_H

#include<cassert>
#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<cmath> //ADDITION 042521 for massaction reactions

#include"compartment/baseCompartmentNeighborhood.h"
#include"reactions/baseReaction.h"
#include"compartment.h"
#include"common/myTimes.h"
#include"species.h"
#include"topology.h"
#include"common/typedefs.h"


///
/// @brief The organism class defines the complete model
///
/// The Organism handles the update and information of the topology
/// and species (molecules) within each compartment (e.g. cell). It is
/// the 'top' class for the defined model. It includes reactions
/// updating more than one species simultaneously. It also holds the
/// neighborhood if defined.
class Organism {
	
 private:
	std::string id_;
	std::vector<Compartment> compartment_; //Modification: AD 032117
	std::vector<int> divisionIndex_;
	std::vector<double> divisionTime_;
	Topology topology_;
	unsigned int numTopology_;
	std::vector<Species> species_;
	std::vector<BaseReaction*> reaction_;
	std::vector<double*> parameter_;
	BaseCompartmentNeighborhood *neighborhood_;
	unsigned int numNeighborhood_;
	
	///
	/// @brief Copy constructor
	///
	/// Defined private for avoiding unintended use since it is not yet
	/// implemented.
	///
	Organism( const Organism & organismCopy );
	
public:


    std::string filename_; //Addition: AD 051122 Limitation: Will always output based on name of model file.

    int simPrintFlag; //Addition: AD 060722 Tracking sim printing from Organism
    double stepsize = 0;
       double endTime = 0;
       double timePoint = 0;
       int stochasticActionFlag = 0; //ADDITION: 080222 flag to activate stochastic action once per cycle

    ///
	/// @brief Empty constructor
	///
	/// This constructor only sets the number of topologies and
	/// neighborhoods to zero.
	///
	Organism();
	
	///
	/// @brief Constructor reading a model file and optionally an init file
	///
	/// This constructor reads files and the names are givn by char*'s
	///
	Organism( const char *modelFile, const char *initFile="", int verbose=0 );
	
	///
	/// @brief Constructor reading a model file and optionally an init file
	///
	/// This constructor reads files and the names are givn by strings
	///
	Organism( const std::string &modelFile, const std::string &initFile="", int verbose=0 );

	///
	/// @brief Constructor reading a model file (already opened) and optionally an init file
	///
	/// This constructor reads files and the names are givn by strings
	/// The main advantage is, for an optimization for example, to put everything into one file
	/// (model , solver, estimator and optimizer), and call successively the constructors with
	///  the same open file.
	///
	Organism( const std::ifstream& IN, const char *initFile="", int verbose=0 );

	///
	/// @brief Destructor
	///
	~Organism();
	
	///
	/// @brief Returns the identity string for the model
	///
	inline std::string id() const;
	
	///
	/// @brief Returns the number of dynamical variables in each compartment
	///
	inline size_t numVariable() const;
	
	///
	/// @brief Returns the number of compartments in the organism
	///
	inline size_t numCompartment() const;
	
	///
	/// @brief Returns the number of divisions
	///
	inline size_t numDivision() const;
	
	///
	/// @brief Returns the dimension of the embedding
	///
	inline unsigned int numDimension() const;
	
	///
	/// @brief Returns the number of different topologies 
	///
	/// As of now, at most one topology is allowed for.
	///
	inline unsigned int numTopology() const;
	
	///
	/// @brief Returns the number of variables in the topology
	///
	inline size_t numTopologyVariable() const;
	
	///
	/// @brief Returns the number of species
	///
	inline size_t numSpecies() const;
	
	///
	/// @brief Returns the number of reactions
	///
	inline size_t numReaction() const;
	
	///
	/// @brief Returns the number of neighborhood rules 
	///
	/// As of now, at most one is allowed for
	///
	inline unsigned int numNeighborhood() const;
	
	///
	/// @brief Returns the compartment with index i for const use
	///
	inline const Compartment &compartment(size_t i) const;
	
	///
	/// @brief Returns the compartment with index i
	///
	inline Compartment &compartment(size_t i); //MODIFIED AD 032117
	
	///
	/// @brief Returns the compartment index of the ith division
	///
	inline int divisionIndex(size_t i) const;
	
	///
	/// @brief Returns the time of the ith division
	///
	inline double divisionTime(size_t i) const;
	
	///
	/// @brief Returns the Topology for const use
	///
	inline const Topology &topology() const;
	
	///
	/// @brief Returns the Topology
	///
	inline Topology &topology();
	
	///
	/// @brief Returns a pointer to the topology
	///
	inline Topology* topologyPointer();
	
	///
	/// @brief Returns the Species with index i for const use
	///
	inline const Species &species(size_t i) const;
	
	///
	/// @brief Returns the Species with index i
	///
	inline Species &species(size_t i);
	
	///
	/// @brief Returns a reaction (pointer)
	///
	inline BaseReaction* reaction(size_t i) const;
	
	///
	/// @brief Returns the reaction pointers vector
	///
	inline const std::vector<BaseReaction*> reaction() const;
	
	///
	/// @brief Returns the Neighborhood for const use
	///
	inline const BaseCompartmentNeighborhood &neighborhood() const;
	
	///
	/// @brief Returns the Neighborhood
	///
	inline BaseCompartmentNeighborhood &neighborhood();
	
	///
	/// @brief Returns a pointer to the neighborhood
	///
	inline BaseCompartmentNeighborhood* neighborhoodPointer();
	
	///
	/// @brief Returns the number of parameters for the complete organism
	///
	inline size_t numParameter() const;
	
	///
	/// @brief Returns a parameter value
	///
	inline double parameter(size_t i) const;
	
	// Set and add variables 
	//
	///
	/// @brief Sets the id string
	///
	inline void setId(const std::string &value);
	
	///
	/// @brief Copy compartment value into position i
	///
	inline void setCompartment(size_t i,Compartment &value); //MODIFIED AD 032117
	
	///
	/// @brief Adds a compartment
	///
	inline void addCompartment(Compartment &value); //MODIFIED AD 032117
	
	///
	/// @brief Adds a compartment division index and time
	///
	inline void addDivision(int indexValue,double timeValue);
	
	//Sets the Topology
	//inline void setTopology(Topology &value);
	
	///
	/// @brief Sets the number of topologies
	///
	inline void setNumTopology(unsigned int value);
	
	///
	/// @brief Copy the species value into position i 
	///
	inline void setSpecies(size_t i,Species &value);
	
	///
	/// @brief Adds a Species from a class instance
	///
	inline void addSpecies(Species &value);
	
	///
	/// @brief Adds a Species from reading in an open file
	///
	inline void addSpecies(std::ifstream &IN);
	
	//Sets the neighbors for compartment with index i
	//inline void setNeighbors(size_t i,std::vector<int> &value);
	
	///
	/// @brief Sets the reaction(pointer) with index i
	///
	inline void setReaction(size_t i,BaseReaction* value);
	
	///
	/// @brief Sets all reaction(pointers) for a species
	///
	inline void setReaction(const std::vector<BaseReaction*> value);
	
	///
	/// @brief Adds a reaction to the list from an open file
	///
	inline int addReaction( std::istream &IN );
	
	// Sets a neighborhood
	//inline void setNeighborhood(BaseCompartmentNeighborhood *value);
	
	///
	/// @brief Adds the neighborhood from reading an open file
	///
	inline int addNeighborhood( std::istream &IN );
	
	///
	/// @brief Sets the number of neighborhoods
	///
	inline void setNumNeighborhood(size_t value);
	
	///
	/// @brief Sets the parameter value for parameter i
	///
	inline void setParameter(size_t i,double value);
	
	///
	/// @brief Returns the identity string for a variable of the model
	///
	std::string variableId(size_t varIndex) const;

	///
	/// @brief Reads a model and init file into the organism
	///
	/// This function converts the filenames to char* and calls that
	/// readOrganism function.
	///
	void readOrganism( const std::string &modelFile, 
			   const std::string &initFile,
			   int verbose=0 );
			   
	
	///
	/// @brief Reads a model and init file into the organism
	///
	/// This function calls the readModel and readInit functions (id
	/// init file provided). If verbose is set, it will print
	/// information to stderr.
	///
	void readOrganism( const char *modelFile, const char *initFile, int verbose=0 );

	void readOrganism( std::ifstream &IN, const char *initFile, int verbose=0 );
	///
	/// @brief Reads a model from an open file
	///
	/// This is the readModel function implementing the actual
	/// reading. First it reads a number of parameters, a topology with
	/// its reactions and compartmentchanges. Then it reads all species
	/// (molecules) with their reactions, followed by reading additional
	/// reactions. Finally it possibly reads a neighborhood.
	///
	/// Here follows a detailed discussion on the parameters in the model,
	/// divided into reasonable sub-parts.
	///
	/// ----------------------------------------<br>
	/// ModelName @f$ N_t N_s N_r N_n @f$ <br>
	/// ----------------------------------------<br>
	///
	/// The first line in the model file reads an odentity string for the model
	/// (which is only for the user), and then reads the number of different
	/// structures to be presented in the model file.
	///
	/// <ul> 
	/// 
	/// <li> @b @f$ N_t @f$ - number of topologies, which handles the positional
	/// and size variables, including cell divisions and removal. In the present
	/// version this can be either zero or one.
	///
	/// <li> @b @f$ N_s @f$ - number of molecular species included in the model.
	///
	/// <li> @b @f$ N_r @f$ - number of reaction in the model. Note that reactions
	/// also can be defined within the definition of the topology and molecular
	/// species. 
	///
	/// <li> @b @f$ N_n @f$ - number of neighborhood rules defined. Only zero or
	/// one is available at present.
	///
	/// </ul>
	///
	/// Then it is time to read the topology:
	///
	/// ----------------------------------------<br>
	/// TopologyName @f$ N_v^t N_r^t N_c^t @f$ <br>
	/// reaction1 @f$ N_p^{r1} N_{il}^{r1} N_{i1}^{r1} @f$ ... <br>
	/// @f$ p_1^{r1} p_2^{r1} ... @f$ <br>
	/// @f$ i_{11}^{r1} i_{12}^{r1} ... @f$ <br>
	/// @f$ i_{21}^{r1} i_{22}^{r1} ... @f$ <br>
	/// ... <br>
	/// reaction2 @f$ N_p^{r2} N_{il}^{r2} N_{i1}^{r2} @f$ ... <br>
	/// @f$ p_1^{r2} p_2^{r2} ... @f$ <br>
	/// @f$ i_{11}^{r2} i_{12}^{r2} ... @f$ <br>
	/// @f$ i_{21}^{r2} i_{22}^{r2} ... @f$ <br>
	/// ... <br>
	/// compartmentChange1 ... <br>
	/// ... <br>
	/// ----------------------------------------<br>
	///
	/// The fist line reads a user-defined name of the topology followed by the
	/// number of variables (@f$ N_v^t @f$) used for describing the topology of
	/// each compartment. This number is the dimension plus one for spherical
	/// cells, two times the dimension plus one for bacterial cells, etc. Then
	/// the number of reactions acting on the topology variables, @f$ N_r^t @f$
	/// and the number of compartment changes defined, @f$ N_c^t @f$, are
	/// read. The format of the reactions can be found below. The compartment
	/// changes defines rules for cell division, cell removal etc. ADD
	/// DESCRIPTION ON COMPARTMENT CHANGE FORMAT HERE. Available compartment
	/// changes are the classes inheriting BaseCompartmentChange.
	///
	/// Each species (@f$ N_s @f$) is read as:
	///
	/// ----------------------------------------<br>
	/// SpeciesName @f$ i_s N_r^t @f$ <br>
	/// reaction1 @f$ N_p^{r1} N_{il}^{r1} N_{i1}^{r1} @f$ ... <br>
	/// @f$ p_1^{r1} p_2^{r1} ... @f$ <br>
	/// @f$ i_{11}^{r1} i_{12}^{r1} ... @f$ <br>
	/// @f$ i_{21}^{r1} i_{22}^{r1} ... @f$ <br>
	/// ... <br>
	/// reaction2 @f$ N_p^{r2} N_{il}^{r2} N_{i1}^{r2} @f$ ... <br>
	/// @f$ p_1^{r2} p_2^{r2} ... @f$ <br>
	/// @f$ i_{11}^{r2} i_{12}^{r2} ... @f$ <br>
	/// @f$ i_{21}^{r2} i_{22}^{r2} ... @f$ <br>
	/// ... <br>
	/// ----------------------------------------<br>
	///
	/// The first line reads a user defined name for the species followed by the
	/// variable index of this molecule, and the number of reactions defined
	/// here for this molecule. The index should be an increasing number in the
	/// list of species starting at the number of topology variables (see
	/// above), since the topology uses the indises @f$ [0...N_v^t-1] @f$. Then
	/// each reaction is defined as described below.
	///
	/// @f$ N_r @f$ reactions are read as:
	///
	/// ----------------------------------------<br>
	/// reaction1 @f$ N_p^{r1} N_{il}^{r1} N_{i1}^{r1} @f$ ... <br>
	/// @f$ p_1^{r1} p_2^{r1} ... @f$ <br>
	/// @f$ i_{11}^{r1} i_{12}^{r1} ... @f$ <br>
	/// @f$ i_{21}^{r1} i_{22}^{r1} ... @f$ <br>
	/// ... <br>
	/// ----------------------------------------<br>
	///
	/// Possible reactions can be found among the classes inheriting
	/// BaseReaction.
	///
	/// Finally, if @f$ N_n=1 @f$ a rule for defining a neighborhood is read:
	///
	/// ----------------------------------------<br>
	/// neighborhoodName ... <br>
	/// ... <br>
	/// ----------------------------------------<br>
	///
	/// Available neighborhood definitions can be found in classes inheriting
	/// BaseCompartmentNeighborhood.
	///
	/// Comments can be included in the file by starting a row with #.
	/// 
	/// @note No complete check on the validity of the provided values are given.
	///
	void readModel(std::ifstream &IN);
	
	///
	/// @brief Reads a model by calling the ifstream version.
	///
	void readModel(const std::string &fileName);
	
	///
	/// @brief Reads a model by calling the string version which calls the ifstream version.
	///
	void readModel(const char *fileName);
	
	///
	/// @brief Reads an initial configuration from an open file
	///
	/// This function implements the reading of an init file. An init file has
	/// the format:
	///
	/// @verbatim
	/// N M
	/// x_ij ...
	/// ...
	/// @endverbatim
	///
	/// where N is the number of compartments and M is the number of
	/// variables in each compartment, i.e. topology variables (positions
	/// and sizes) and species (molecules). Then all variables are given
	/// by x_ij where i is the compartment and j is the variable.
	///
	/// When the input is read, M is checked to agree with the number of
	/// variables in the model (organism), and then N compartments are created
	/// and initiated with values from x_ij.
	///
	/// Comments can be included in the file by starting a row with #. Caveat:
	/// No check on the validity of the initial variable values provided in x_ij
	/// is done.
	/// 
	void readInit(std::ifstream &IN);
	
	///
	/// @brief Opens the file fileName and calls readInit(ifstream)
	///
	void readInit(const std::string &fileName);
	
	///
	/// @brief Opens the file fileName and calls readInit(ifstream)
	///
	void readInit(const char *fileName);
	
	///
	/// @brief Sets init values from a matrix
	///
	/// sets the init values from an matrix input. It checks that the
	/// number of variables in each compartment (number of columns) is
	/// equal with the number in the model (organism). Then old
	/// compartments are removed and new added.
	///
	void setInit(const DataMatrix &input); 
	
	///
	/// @brief Initiates the parameter vector for all parameters within the model (organism).
	///
	/// This function collects all parameters in all species>reactions
	/// and reactions into a single vector containing all the parameters
	/// for the complete organism (model).
	///
	void initiateParameter();
	
	///
	/// @brief Creates the neighborhood according to the defined neighborhood rule
	///
	/// In the member class neighborhood_ a rule for how to extract and
	/// update the neighborhood for all compartments are defined. 
	///
	/// It calls the corresponding function defined in the neighborhood.
	///
	/// @see BaseCompartmentNeighborhood::create()
	///
	void neighborhoodCreate(DataMatrix &y,double t=0.0 );
	
	///
	/// @brief Updates the neighborhood according to the defined neighborhood rule
	///
	/// This function calls the corresponding function in neighborhood.
	///
	/// @see BaseCompartmentNeighborhood::update()
	///
	void neighborhoodUpdate(DataMatrix &y,double t=0.0 );
	
	///
	/// @brief Removes a compartment and updates the neighborhood
	///
	/// This function can be used when simulating cell migration (out of
	/// the system) or cell death.
	///
	void removeCompartment(size_t i);
	
	///
	/// @brief Calculates the derivatives given variable values in y 
	///
	/// This is the main derivatives function used when integrating the
	/// system. It calls the derivs functions for all its subclasses and
	/// adds it all up.
	///
	void derivs(DataMatrix &y,DataMatrix &dydt);


	///
	/// @brief Calculates the derivatives given variable values in y and saves derivatives and absolute values
	///
	/// A derivatives function used when integrating the
	/// system that also provide the summed absolute values of the variables. This is used by solvers
	/// adding noise to the update. It calls the derivsWithAbs functions for all its subclasses and
	/// adds it all up.
	///
	/// @see HeunIto 
	///
	void derivsWithAbs(DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt);

	///
	/// @brief Updates a Jacobian matrix from the defined reactions
	///
	/// This function is right now only
	/// intended to be used by the ImplictEuler solver 
	/// simulations. It will use information from the individual reactions
	/// and hence relies upon that these have a Jacobian function defined.
	/// It is 'templated' to be able to use e.g. ublas sparse matrix representation.
	///
	/// @see ImplicitEuler
	///
    /*
	size_t Jacobian(DataMatrix &y,JacobianMatrix &A);
    */
	///
	/// @brief Calculates the present derivatives for mechanical interactions only
	///
	/// This function is used in simulations where the mechanics is
	/// equilibriated at each time step. It calls the derivs functions
	/// for all mechanical updates as defined in mechEq.Caveat: The
	/// check that mechEq only updates positional variables are left to
	/// the user!
	///
	void derivsMechanical(DataMatrix &y,DataMatrix &dydt,
			      std::vector<size_t> &mechEq,std::vector<size_t> &posVar);
	
	///
	/// @brief Calculates the derivatives for all species flagged in the simFlag vector 
	///
	/// This function is the same as the derivatives function, but it
	/// only updates a subset of species, assuming that the rest are
	/// read from a template.
	///
	void derivsTemplate(DataMatrix &y,DataMatrix &dydt,
			    std::vector<int> &simulationFlag );
	
	///
	/// @brief Initiates all rections that requires it
	///
	/// Some reactions have initiations(updates) that are done before each 
	/// simulation. This function loops over all reactions and calls the individual
	/// reaction's initiate function. Most reactions do not require an update, and
	/// the baseReaction version does not do anything.
	///
	/// @see BaseReaction::initiate()
	///
	void reactionInitiate(double t,DataMatrix &y);

	///
	/// @brief Updates all rections that requires it
	///
	/// Some reactions have updates that are done after each simulator
	/// step. This function loops over all reactions and calls the individual
	/// reaction's update function. Most reactions do not require an update, and
	/// the baseReaction version does not do anything.
	///
	/// @see BaseReaction::update()
	///
	void reactionUpdate( double h, double t,DataMatrix &y);

	///
	/// @brief Checks for all possible compartment changes as divisions and removals.
	///
	/// This function calls the topology::compartmentChange() to check
	/// whether a compartment is to be updated (removed or
	/// divided). This is done after each step in the numerical solver.
	///
	void compartmentChange(DataMatrix &y,DataMatrix &dydt,double t=0.);
	
	///
	/// @brief Generates a list of reactions included in the model.
	///
	/// Compiles a list of all reactions used in the
	/// model. Typichally to be used by a stochastic solver.
	///
	/// @note In the current implementation, only molecular
	/// reactions are included
	///
	void createReactionList(std::vector<BaseReaction*> &reactionList,
				std::vector<size_t> &speciesList);

	///
	/// @brief Scales all space variables such as positions, areas and volumes
	///
	/// The function scales positions (and lengths) with the factor,
	/// areas with the factor to the power (dimension-1), and volumes
	/// with the factor to the power dimension. It is used when a
	/// template of 'wrong' units is used.
	///
	void scaleSpace( double spaceScalingFactor);
	
	///
	/// @brief Finds peaks in a variable using gradient ascent
	///
	unsigned int findPeaksGradientAscent(const DataMatrix &y,size_t col, 
					     std::vector<size_t> &cellMax,
					     std::vector<size_t> &flag );
	

	///
	/// @brief Prints the model including parameters in the specific model file
	/// format
	///
	void printModel(std::ostream &os=std::cerr) const;

	///
	/// @brief Prints the model including parameters in Cambium format to specific model file
	/// format
	///
	void printModelCambium ( std::ostream &os=std::cerr) const;
	
	///
	/// @brief Prints all model parameters in a space delimited row. 
	///
	void printParameter(std::ostream &os=std::cerr) const;
	
	///
	/// @brief Prints the current values of the variables in an init format
	///
	void printVariable(std::ostream &os=std::cout) const;
	
	///
	/// @brief Prints the current neighborhood variables
	///
	void printNeighbor(std::ostream &os=std::cout) const;	
};

inline std::string Organism::id() const 
{
  return id_;
}

inline size_t Organism::numVariable() const 
{
  return numTopologyVariable()+numSpecies();
}

inline size_t Organism::numCompartment() const 
{
  return compartment_.size();
}

inline size_t Organism::numDivision() const 
{
  return divisionIndex_.size();
}

inline unsigned int Organism::numDimension() const 
{
  return topology_.numDimension();
}

inline unsigned int Organism::numTopology() const 
{
  return numTopology_;
}

inline size_t Organism::numTopologyVariable() const 
{
  return topology_.numVariable();
}

inline size_t Organism::numSpecies() const 
{
  return species_.size();
}

inline size_t Organism::numReaction() const 
{
  return reaction_.size();
}

inline unsigned int Organism::numNeighborhood() const 
{
  return numNeighborhood_;
}

inline size_t Organism::numParameter() const 
{
  return parameter_.size();
}

inline BaseReaction* Organism::reaction(size_t i) const 
{
  return reaction_[i];
}

inline const std::vector<BaseReaction*> Organism::reaction() const 
{
  return reaction_;
}

inline const Compartment& Organism::compartment(size_t i) const 
{
  return compartment_[i];
}

inline Compartment& Organism::compartment(size_t i) 
{
  return compartment_[i];
}

inline int Organism::divisionIndex(size_t i) const 
{
  return divisionIndex_[i];
}

inline double Organism::divisionTime(size_t i) const 
{
  return divisionTime_[i];
}

inline const Topology& Organism::topology() const 
{
  return topology_;
}

inline Topology& Organism::topology() 
{
  return topology_;
}

inline Topology* Organism::topologyPointer() 
{
  return &topology_;
}

inline const Species& Organism::species(size_t i) const 
{
  return species_[i];
}

inline Species& Organism::species(size_t i) 
{
  return species_[i];
}

inline const BaseCompartmentNeighborhood& Organism::neighborhood() const 
{
  return *(neighborhood_);
}

inline BaseCompartmentNeighborhood& Organism::neighborhood() 
{
  return *(neighborhood_);
}

inline BaseCompartmentNeighborhood* Organism::neighborhoodPointer() 
{
  return neighborhood_;
}

inline double Organism::parameter(size_t i) const 
{
  return *(parameter_[i]);
}

inline void Organism::setId(const std::string &value) 
{
  id_ = value;
}

inline void Organism::setCompartment(size_t i,Compartment &value) //MODIFIED AD 032117
{
  compartment_[i]=value; 
}

inline void Organism::addCompartment(Compartment &value) //MODIFIED AD 032117
{
  compartment_.push_back( value );
}

inline void Organism::addDivision(int indexValue,double timeValue) 
{
  divisionIndex_.push_back( indexValue );
  divisionTime_.push_back( timeValue );
}

inline void Organism::setNumTopology(unsigned int value) 
{
  numTopology_=value;
}

//inline void Organism::setTopology(Topology &value) 
//{
//std::cerr << "Organism::setTopology(Topology &value) not yet defined!\n";
//exit(0);
//}

inline void Organism::setSpecies(size_t i,Species &value) 
{
  std::cerr << "Organism::setSpecies(size_t i,Species &value) "
	    << "not yet defined!\n";
  exit(0);
}

inline void Organism::addSpecies(Species &value) 
{
  species_.push_back(value);
}

inline void Organism::addSpecies(std::ifstream &IN) 
{
  species_.push_back( Species(IN) );
}

inline void Organism::setNumNeighborhood(size_t value) 
{
  numNeighborhood_=value;
}

int Organism::addNeighborhood( std::istream &IN ) 
{
  if( !IN )
    return -1;
  neighborhood_=BaseCompartmentNeighborhood::
    createCompartmentNeighborhood(IN);
  return 0;
}

inline void Organism::setReaction(size_t i,BaseReaction* value) 
{
  reaction_[i]=value;
}

inline void Organism::setReaction( const std::vector<BaseReaction*> value) 
{
  reaction_ = value;
}

//inline void Organism::setNeighbors(size_t i,std::vector<int> &value) 
//{
//std::cerr << "Organism::setNeighbors(size_t i,std::vector<int> &value) "
//<< "not yet defined!\n";
//exit(0);
//}

inline void Organism::setParameter(size_t i,double value) 
{
  *(parameter_[i])=value;
}

#endif
