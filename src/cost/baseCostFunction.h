//
// Filename     : baseCostFunction.h
// Description  : Classes defining cost functions
// Author(s)    : Pontus Melke (pontus@thep.lu.se)
// Created      : August 2006
// Revision     : $Id: baseCostFunction.h 657 2016-01-21 18:48:05Z korsbo $
//

#ifndef BASECOSTFUNCTION_H
#define BASECOSTFUNCTION_H

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "../common/myFiles.h"

///
///@brief A base class using a factory method to create a cost function. The
///class Solver then has a pointer to a specific cost function allowing it to
///comparing simulations to experimental data.
///
/// The BaseCostFunction class is a base (factory) class used when defining
/// different types of "cost function" classes. The name of the data file
/// containing experimental data is stored in costTemplateFile_ , the list of
/// species to be used in the cost calculation is stored in costList_ and the
/// time points at which the cost is to be calculated is stored in
/// costTemplateTime_.
///
/// Each derived class must implement a funtion called getCost which defines
/// the desired cost function. The function getCost takes three arguments -
/// simValue, which is the value obtained from simulation, dataValue which is
/// the value obtained from experiment and the third parameter, scaleFactor,
/// is optional and it merely rescale the cost value and is by default set to
/// one -
///
class BaseCostFunction {
 protected:
	std::vector<size_t> costList_;
	std::string costTemplateFile_;
	double cost_;
	size_t costTemplateIndex_;
	std::vector<double> costTemplateTime_;
	std::vector< std::vector< std::vector<double> > > costTemplate_;
	std::vector< std::vector<double> > costIndividual_;

 public:
	virtual ~BaseCostFunction();

	///
	/// @brief Factory for initiating a cost function to be used by a simulation
	/// or optimization.
	///
	/// This is the factory method for generating a cost function to measure the
	/// cost/ objective function value, energy,... for a simulation which is
	/// also used when optimizing model parameters. An id string is provided and
	/// generates the initiation of a correct cost function (class). Examples are: 
	///
	/// <pre>
	/// meanSquare
	/// meanSquareRelative
	/// </pre>
	///
	/// @see MeanSquare
	/// @see MeanSquareRelative
	/// @see RK5AdaptiveTemplate::readParameterFile()
	/// @see BaseSolver::getSolver()
	///
	static BaseCostFunction *createCostFunction(const std::string &costType);
	virtual double getCost(double simValue, double dataValue, double scaleFactor = 1.0  );
	void initiateCostCalculation(double startTime, double endTime);
	virtual double addCost(const std::vector< std::vector<double> > &y, const double t);
	
	inline double endTime();
	inline size_t sizeOfCostList();
	inline size_t sizeOfCostTemplateTime();
	inline size_t costTemplateIndex();
	inline double costTemplateTime(size_t i);
	inline double cost() const;
	inline void resetCostCalculation();
	inline size_t costList(size_t i);

	inline void setCostTemplateFile( const std::string &value );
	inline void setCostTemplate(const std::vector< std::vector< std::vector<double> > > &costTemplate );
	inline void setCostTemplateTime(const std::vector<double> &inputVector);
	inline void setCostList(const std::vector<size_t> &inputVector);

	void printCost(std::ostream &os=std::cerr) const;
};

inline double BaseCostFunction::endTime()
{
	return costTemplateTime_[costTemplateTime_.size()-1];
}

inline size_t BaseCostFunction::sizeOfCostList()
{
	return costList_.size();
}

inline size_t BaseCostFunction::sizeOfCostTemplateTime()
{
	return costTemplateTime_.size();
}

inline size_t BaseCostFunction::costTemplateIndex()
{
	return costTemplateIndex_;
}

inline double BaseCostFunction::costTemplateTime(size_t i) 
{
	return costTemplateTime_[i];
}

inline double BaseCostFunction::cost() const 
{ 
	if (costTemplateFile_.size())
		return cost_;
	else { 
		return -1;
	}
}

inline void BaseCostFunction::resetCostCalculation()
{
	cost_=0;
	costTemplateIndex_=0;
}

inline size_t BaseCostFunction::costList(size_t i)
{
	return costList_[i];
}

inline void BaseCostFunction::setCostTemplateFile( const std::string &value ) 
{
  costTemplateFile_=value;
}

inline void BaseCostFunction::setCostTemplate(const std::vector< std::vector< std::vector<double> > > &input )
{
	costTemplate_=input;
}

inline void BaseCostFunction::setCostTemplateTime(const std::vector<double> &input) 
{
	costTemplateTime_=input;
}
inline void BaseCostFunction::setCostList(const std::vector<size_t> &input) 
{
	costList_=input;
}
#endif
