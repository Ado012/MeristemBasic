//
// Filename     : baseCostFunction.cc
// Description  : Classes defining cost functions
// Author(s)    : Pontus Melke (pontus@thep.lu.se)
// Created      : August 2006
// Revision     : $Id: baseCostFunction.cc 670 2016-08-03 15:06:43Z korsbo $
//

#include "baseCostFunction.h"
#include "costFunction.h"
#include <cstdlib>

BaseCostFunction::~BaseCostFunction(){
}

double BaseCostFunction::getCost(double simValue, double dataValue, double scaleFactor)
{
  std::cerr << "Warning BaseCostFunction::getCost() should be implemented in inherited classes\n";
  exit(-1);
}

BaseCostFunction* BaseCostFunction::
createCostFunction(const std::string &costType)
{
  BaseCostFunction *tmp;
  ////////////////////////////////////////////////////////////////
  if (costType == "meanSquare")
    tmp = new MeanSquare;
  else if (costType == "manhattan")
    tmp = new Manhattan;
  else if (costType == "meanSquareRelative")
    tmp = new MeanSquareRelative;
  else if (costType == "normalizedMeanSquare")
    tmp = new NormalizedMeanSquare;
  else if (costType == "maxSquare")
    tmp = new MaxSquare;
  else if (costType == "domainWeightedMeanSquare")
    tmp = new DomainWeightedMeanSquare;
  else if (costType == "cosine")
    tmp = new Cosine;
  else if (costType == "cosineSpatial")
    tmp = new CosineSpatial;
  else if (costType == "cosineCell")
    tmp = new CosineCell;
  else if (costType == "mscs")
    tmp = new mscs;
  else if (costType == "KullbackLeibler" or costType == "relativeEntropy")
    tmp = new KullbackLeibler;
  ////////////////////////////////////////////////////////////////
  else {
    std::cerr << "BaseCostFunction* BaseCostFunction::createCostFunction() - " 
	      << "Unknown solver: " << costType << std::endl;
    exit(-1);
  }
  return tmp;
}

void BaseCostFunction::initiateCostCalculation(double startTime, double endTime) 
{
  if( costList_.size()==0 ) {
    std::cerr << "Simulator::initiateCostCalculation() "
	      << "No variables defined to calculate cost from.\n";
    exit(-1);
  }
  
  cost_=0.;
  costTemplateIndex_=0;
  costTemplateTime_.resize(0);
  costTemplate_.resize(0);    
  const char* tmp = costTemplateFile_.c_str();
  std::istream *IN = myFiles::openFile(tmp);
  if( !IN ) {
    std::cerr << "BaseCostFunction::initiateCostCalculation() "
	      << "Cannot open cost templatefile " << costTemplateFile_
	      << "\n\n\7";exit(-1);}
  int count=0;
  double tmpVal;
  size_t N1,M1,Ndiv1;
  double time1;
  
  do {
    //not first time
    if( count++ ) {
      *IN >> Ndiv1;
      for (size_t i=0 ; i<Ndiv1 ; i++ ) {
        *IN >> tmpVal;
        *IN >> tmpVal;
      }
    }
    *IN >> time1;
    *IN >> N1;
    *IN >> M1;
    if( time1>=startTime && *IN ) {
      //Start reading file
      costTemplateTime_.push_back( time1 );
      
      std::vector< std::vector<double> > tmp( N1 );
      for (size_t i=0 ; i<N1 ; i++) {
        size_t k=0;
        tmp[i].resize( costList_.size() );
        for (size_t j=0 ; j<M1 ; j++ ) {
          *IN >> tmpVal;
          if( k<costList_.size() && j == size_t(costList_[k]) ) {
            tmp[i][k] = tmpVal;
	    k++;
          }
        }
      }
      costTemplate_.push_back( tmp );
      
    }
  } while( time1<=endTime && *IN );
  delete IN;
  costIndividual_.resize(costTemplate_.size());
  for (size_t c=0; c<costIndividual_.size(); ++c)
    costIndividual_[c].resize(costList_.size());
}

//!Adds cost from applicable time point
double BaseCostFunction::addCost(const std::vector< std::vector<double> > &y, const double t) 
{ 
  if( !costList_.size() ) {
    std::cerr << "BaseCostFunction::addCost() No variables marked for "
	      << "cost calculation. Nothing added to cost.\n";
    return 0.;
  }
  if( costTemplateIndex_<0 || costTemplateIndex_>=costTemplateTime_.size() ) {
    std::cerr << "BaseCostFunction::addCost() Wrong template index for "
	      << "cost calculation. Nothing added to cost.\n";
    return 0.;    
  }  
  if( costTemplate_[costTemplateIndex_].size() != y.size() ) {
    std::cerr << "BaseCostFunction::addCost() Not correct number of cells "
	      << "to compare simulation with cost template\n";
    std::cerr << "Number of cells in cost template: " << costTemplate_[costTemplateIndex_].size()
	      << ", Number of cells in simulation result: " << y.size() << "\n";
    exit(-1);
  }
  if( std::fabs( t-costTemplateTime_[costTemplateIndex_] ) > 0.001 ) {
    std::cerr << "BaseCostFunction::addCost() Warning: time difference "
	      << "between simulation and cost template (" 
	      << t << " " << costTemplateTime_[costTemplateIndex_] << "  "
	      << "costTemplateIndex: " << costTemplateIndex_ << " "
	      << std::fabs( t-costTemplateTime_[costTemplateIndex_] )
	      << ")\n";
    exit(-1);
  }
  double costAdd=0.;
  double costAddTotal=0.;
  // sum over species used in the calculation
  for (size_t k=0 ; k<costList_.size() ; ++k ) {
    int j=costList_[k];
    // sum over cells
    costAdd=0.;
    for (size_t i=0 ; i<y.size() ; ++i ) {
      costAdd += getCost(y[i][j],costTemplate_[costTemplateIndex_][i][k]);
    }
    costAddTotal += costAdd;

		// the following if-statement prevents segmentation fault of the
		// optimizer if the size of the a cost template used in the 
		// estimator file is larger than that specified in the solver file.
		if ( costIndividual_.size() > costTemplateIndex_ ) 
			costIndividual_[costTemplateIndex_][k] = (costAdd)/(y.size()); // MSE for species k at costTemplateIndex_ 
  }
  // double costAddTotalMean = costAddTotal/(y.size()*costList_.size()*costTemplate_.size());
  // costList_.size() gives the number of variables for cost calculation of the current run
  // costTemplate_.size() gives the number of time points specified in the cost file.
  double costAddTotalMean = costAddTotal/(y.size()*costList_.size());
  
  cost_ += costAddTotalMean; // cost_ saves the total cost!
  costTemplateIndex_++;
  return costAddTotal;
}

void BaseCostFunction::printCost(std::ostream &os) const
{
  os << "Cost (all template points used): " << cost() << std::endl;
  os << "time\t";
  for (size_t k=0; k<costList_.size(); ++k)
    os << "var_" << costList_[k] << "\t";
  os << "total" << std::endl;
  for (size_t t=0; t<costTemplateTime_.size(); ++t) {
    os << costTemplateTime_[t] << "\t";
    //double tSum=0.0;  // does not always make sense to just add the mean values
    for (size_t k=0; k<costList_.size(); ++k) {
      os << costIndividual_[t][k] << "\t";
      //tSum += costIndividual_[t][k];
    }
    //os << tSum << std::endl;
    os << "-" << std::endl;
  }
  os << "total\t"; 
  for (size_t k=0; k<costList_.size(); ++k) {
    //double kSum=0.0; // does not always make sense to just add the mean values
    for (size_t t=0; t<costTemplateTime_.size(); ++t) {
      //kSum += costIndividual_[t][k];
    }
    //os << kSum << "\t"; 
    os << "-" << "\t";
  }
  os << cost() << std::endl;
}
