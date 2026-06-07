
//
// Filename     : costFunction.cc
// Description  : Classes defining cost functions
// Author(s)    : Pontus Melke (pontus@thep.lu.se)
// Created      : August 2006
// Revision     : $Id: costFunction.cc 656 2016-01-21 14:33:38Z andre $
//

#include "costFunction.h"
#include<cstdlib>

double MeanSquare::getCost(double simValue, double dataValue, double scaleFactor)
{
        //std::cerr << "caclulating cost to: ";
	double costFactor= scaleFactor*(simValue-dataValue);
	//std::cerr << costFactor*costFactor << "\n";
	return costFactor*costFactor;
}

double MeanSquareRelative::getCost(double simValue, double dataValue, double scaleFactor)
{
	double costFactor= scaleFactor*(simValue-dataValue)/dataValue;	
	return costFactor*costFactor;
}

double NormalizedMeanSquare::getCost(double simValue, double dataValue, double scaleFactor)
{
  double costFactor= scaleFactor*(simValue-dataValue)/dataValue;	
  return costFactor*costFactor;
}

double NormalizedMeanSquare::
addCost(const std::vector< std::vector<double> > &y, const double t) 
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
  double costAddPrev=0.0;
  size_t xI=0,yI=1;
  for (size_t k=0 ; k<costList_.size() ; ++k ) {
    int j=costList_[k];
    for (size_t i=0 ; i<y.size() ; ++i ) {
      double normFactor=1.0;
      if (k==0) {//WUSCHEL
	if (y[i][xI]<1 && y[i][yI]<-1)
	  normFactor=78;
	else
	  normFactor=0.1;
      }
      else if (k==1) {//KANADI
      }
      else if (k==2) {//CLAVATA3
      }
      else {
	std::cerr << "NormalizedMeanSquare::addCost() Only defined for clvWusKan...\n";
	exit(-1);
      }
      costAdd += getCost(y[i][j],costTemplate_[costTemplateIndex_][i][k],normFactor);
    } 
    //		costIndividual_[costTemplateIndex_][k] = (costAdd-costAddPrev)/(2*y.size()*costList_.size()*costTemplate_.size());
    costAddPrev = costAdd;
  }
  cost_ += costAdd/(2*y.size()*costList_.size()*costTemplate_.size());
  costTemplateIndex_++;
  return costAdd;
}

double MaxSquare::getCost(double simValue, double dataValue)
{
  double costFactor= simValue-dataValue;	
  return costFactor*costFactor;
}


double MaxSquare::
addCost(const std::vector< std::vector<double> > &y, const double t) 
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
  double current_cost=0.;
  double max_cost = cost_; // assume current cost value is previous maximum
  // Loop over species k
  for (size_t k=0 ; k<costList_.size() ; ++k ) {
    int j=costList_[k];
    // Loop over #cells (?)
    for (size_t i=0 ; i<y.size() ; ++i ) {
      current_cost = getCost(y[i][j],costTemplate_[costTemplateIndex_][i][k]);
      if( current_cost > max_cost){
	max_cost = current_cost;
      }
    } 
  }
  costTemplateIndex_++;
  if(cost_ < max_cost)
    cost_ = max_cost;
  return max_cost;
}

double DomainWeightedMeanSquare::getCost(double simValue, double dataValue)
{
  double costFactor= simValue-dataValue;// /dataValue;  
  return costFactor*costFactor;
}


double DomainWeightedMeanSquare::addCost(const std::vector< std::vector<double> > &y, const double t) 
{ 
  if( !costList_.size() ) {
    std::cerr << "DomainWeightedMeanSquare::addCost() No variables marked for "
    << "cost calculation. Nothing added to cost.\n";
    return 0.;
  }
  if( costTemplateIndex_<0 || costTemplateIndex_>=costTemplateTime_.size() ) {
    std::cerr << "DomainWeightedMeanSquare::addCost() Wrong template index for "
    << "cost calculation. Nothing added to cost.\n";
    return 0.;    
  }  
  if( costTemplate_[costTemplateIndex_].size() != y.size() ) {
    std::cerr << "DomainWeightedMeanSquare::addCost() Not correct number of cells "
    << "to compare simulation with cost template\n";
    std::cerr << "Number of cells in cost template: " << costTemplate_[costTemplateIndex_].size()
    << ", Number of cells in simulation result: " << y.size() << "\n";
    exit(-1);
  }
  if( std::fabs( t-costTemplateTime_[costTemplateIndex_] ) > 0.001 ) {
    std::cerr << "DomainWeightedMeanSquare::addCost() Warning: time difference "
    << "between simulation and cost template (" 
      << t << " " << costTemplateTime_[costTemplateIndex_] << "  "
      << "costTemplateIndex: " << costTemplateIndex_ << " "
      << std::fabs( t-costTemplateTime_[costTemplateIndex_] )
      << ")\n";
exit(-1);
  } //Copied from baseCostFunction. I'm not sure about all that it does so I will not touch it.
  
  double costAddTotal=0.;
  double costAddIndividual=0;
  // The fun begins here

  // I'm somewhat following the notation of the documentation
  // M is the nuber of spicies to optimize against (relevant cloumns of the costtemplate)
  // N is the total number of cells

  size_t M = costList_.size(); // the number of columns in the costtemplate that we optimize against
  size_t N = y.size(); // The amount of cells


  for (size_t m = 0 ; m<M ; m++){ // Loop through all relevant columns of the costtemplate
    costAddIndividual =0;
    double E_domain_1=0;
    double E_domain_2=0;
    size_t N_domain_1=0;
    size_t N_domain_2=0;
    int j=costList_[m];
    for (size_t n = 0 ; n<N ; n++){ //Loop through the number of cells in the template
      if (costTemplate_[costTemplateIndex_][n][m]!= 0){ //If the cost template does not say 0: assign domain identity 1 and compute
        E_domain_1+=getCost(y[n][j],costTemplate_[costTemplateIndex_][n][m]);
        N_domain_1+=1;
      }
      else{
        E_domain_2+=y[n][j]*y[n][j]; //if the target is zero the getCost() reduces to this
        N_domain_2+=1;
      }
    }
    if (N_domain_1*N_domain_2 !=0 ){ //if none of the domains lack cells (i.e. there are two domains) 
      costAddIndividual+=double(N)/(2.*N_domain_1)*E_domain_1 + double(N)/(2.*N_domain_2)*E_domain_2;
      //That was the magic line which assigns equal weight to the two domains
    }
    else if(N_domain_1!=0){
      costAddIndividual+=E_domain_1;
    }
    else if (N_domain_2!=0){
      costAddIndividual+=E_domain_2;
    }
    else{ //this should never occurr
      std::cerr << "DomainWeightedMeanSquare::addCost ERROR - your costtemplate"
      <<"may not have been understood" << "\n";
      std::exit;
    }
    if(N!=N_domain_2+N_domain_1){ //Only for trouble-shooting.
      std::cerr<<"N!=N_domain_1+N_domain_2" <<"\n";
    }
    costIndividual_[costTemplateIndex_][m] = (costAddIndividual)/(y.size()); // MSE for species k at costTemplateIndex_ 
    costAddTotal += costAddIndividual;
  }
  
  double costAddTotalMean = costAddTotal/(y.size()*costList_.size()*costTemplate_.size());
  costTemplateIndex_++;
  cost_ += costAddTotalMean; // cost_ saves the total cost!

  return costAddTotal;
}


double Manhattan::getCost(double simValue, double dataValue, double scaleFactor)
{
        //std::cerr << "caclulating cost to: ";
	double costFactor= (simValue-dataValue);
	//std::cerr << costFactor*costFactor << "\n";
	if(costFactor > 0)
	  return costFactor;
	else
	  return -costFactor;
}

double CosineCell::getCost(std::vector< double > simValue, std::vector< double > dataValue)
{

  if(simValue.size() != dataValue.size()){
    std::cerr << "CosineCell::getCost - Warning! Different number of elements in sim result and data!\n";
    return 0;
  }

  int vector_dim = simValue.size();
  double dot_product = 0;
  for(int i=0; i< vector_dim; i++){    
    dot_product += simValue[i]*dataValue[i];
  }

  double normSim = 0;
  double normData = 0;

  for(int i=0; i< vector_dim; i++){
    normSim += simValue[i]*simValue[i];
    normData += dataValue[i]*dataValue[i];
  }


  double denominator =  std::sqrt(normSim*normData);

  if(denominator == 0) return (3.14159); // return pi, dissimilar
  else if(dot_product > denominator) return 0; // would be outside of range for acos, this can happen due to numerical errors
  else return std::acos(dot_product/denominator);
}


double CosineCell::
addCost(const std::vector< std::vector<double> > &y, const double t) 
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


  size_t num_cost_species = costList_.size();

  std::vector< double > simValue (num_cost_species);
  std::vector< double > dataValue (num_cost_species);


  double costAdd=0.;
  // Loop over cells (i)
  for (size_t i=0 ; i<y.size() ; ++i ) {

    // Loop over species (k) to get the simulation and data vector
    for (size_t k=0 ; k<num_cost_species ; ++k ) {
      int j=costList_[k];
      simValue[k] = y[i][j]; // set to elements from simulation result
      dataValue[k] = costTemplate_[costTemplateIndex_][i][k]; // set with costTemplate elements
      }

    double phi_mod = getCost(simValue, dataValue)/(3.14159); // angle between vectors, normalized
    costAdd += phi_mod/(y.size()*costTemplate_.size());

    } 
  
  costTemplateIndex_++;
  cost_ += costAdd;
  
  return costAdd;
}

double CosineSpatial::getCost(std::vector< double > simValue, std::vector< double > dataValue)
{

  if(simValue.size() != dataValue.size()){
    std::cerr << "Cosine::getCost - Warning! Different number of elements in sim result and data!\n";
    return 0;
  }

  int vector_dim = simValue.size();
  double dot_product = 0;
  for(int i=0; i< vector_dim; i++){    
    dot_product += simValue[i]*dataValue[i];
  }

  double normSim = 0;
  double normData = 0;

  for(int i=0; i< vector_dim; i++){
    normSim += simValue[i]*simValue[i];
    normData += dataValue[i]*dataValue[i];
  }

  double denominator =  std::sqrt(normSim*normData);

  if(denominator == 0) return (3.14159); // return pi, dissimilar
  else if(dot_product > denominator) return 0; // would be outside of range for acos, this can happen due to numerical errors
  else return std::acos(dot_product/denominator);
}


double CosineSpatial::
addCost(const std::vector< std::vector<double> > &y, const double t) 
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


  size_t num_cost_species = costList_.size();
  size_t num_cells = y.size();

  std::vector< double > simValue (num_cells);
  std::vector< double > dataValue (num_cells);


  double costAdd=0.;
  // Loop over species (k)
  for (size_t k=0 ; k<num_cost_species ; ++k ) {
    int j=costList_[k];
    // Loop over cells (i)
    for (size_t i=0 ; i<num_cells ; ++i ) {

      // set to elements from simulation result
      simValue[i] = y[i][j];

      // set with costTemplate elements
      dataValue[i] = costTemplate_[costTemplateIndex_][i][k];
    }
    double phi_mod = getCost(simValue, dataValue)/(3.14159); // angle between vectors, normalized
    costAdd += phi_mod/(num_cost_species*costTemplate_.size());

    costIndividual_[costTemplateIndex_][k] += phi_mod/(costTemplate_.size()); // error for species k at costTemplateIndex_ 

    } 
  

  //std::cerr << "cost add: " << costAdd << "\n";
  costTemplateIndex_++;
  cost_ += costAdd;
  return costAdd;
}


double Cosine::getCost(std::vector< double > simValue, std::vector< double > dataValue)
{

  if(simValue.size() != dataValue.size()){
    std::cerr << "Cosine::getCost - Warning! Different number of elements in sim result and data!\n";
    return 0;
  }

  int vector_dim = simValue.size();
  double dot_product = 0;
  for(int i=0; i< vector_dim; i++){    
    dot_product += simValue[i]*dataValue[i];
  }

  double normSim = 0;
  double normData = 0;

  for(int i=0; i< vector_dim; i++){
    normSim += simValue[i]*simValue[i];
    normData += dataValue[i]*dataValue[i];
  }

  double denominator =  std::sqrt(normSim*normData);

  if(denominator == 0) return (3.14159); // zero length vector, return pi, dissimilar
  else if(dot_product > denominator) return 0; // would be outside of range for acos, this can happen due to numerical errors
  else return std::acos(dot_product/denominator);
}


double Cosine::
addCost(const std::vector< std::vector<double> > &y, const double t) 
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


  size_t num_cost_species = costList_.size();
  size_t num_cells = y.size();

  std::vector< double > simValue (num_cells*num_cost_species);
  std::vector< double > dataValue (num_cells*num_cost_species);


  double costAdd=0.;
  // Loop over species (k)
  size_t vector_index=0;
  //std::cerr << "\n";
  for (size_t k=0 ; k<num_cost_species ; ++k ) {
    int cost_species_index=costList_[k];
    // Loop over cells (i)
    for (size_t i=0 ; i<num_cells ; ++i ) {

      // set to elements from simulation result
      simValue[vector_index] = y[i][cost_species_index];

      // set with costTemplate elements
      dataValue[vector_index] = costTemplate_[costTemplateIndex_][i][k];
      vector_index++;

      //std::cerr << "vi: " << vector_index<< "\n";
    }
  }

  //std::cerr << "1\n";

  double phi_mod = getCost(simValue, dataValue)/(3.14159); // angle between vectors, normalized
  costAdd += phi_mod/(costTemplate_.size());

  //std::cerr << "cost add: " << costAdd << "\n";
  costTemplateIndex_++;
  cost_ += costAdd;
  return costAdd;
}



double mscs::getCost(std::vector< double > simValue, std::vector< double > dataValue)
{

  if(simValue.size() != dataValue.size()){
    std::cerr << "Cosine::getCost - Warning! Different number of elements in sim result and data!\n";
    return 0;
  }

  int vector_dim = simValue.size();
  double dot_product = 0;
  double sse = 0;
  for(int i=0; i< vector_dim; i++){    
    dot_product += simValue[i]*dataValue[i];
    sse += ((simValue[i]-dataValue[i])*(simValue[i]-dataValue[i]))/vector_dim;
  }

  double normSim = 0;
  double normData = 0;

  for(int i=0; i< vector_dim; i++){
    normSim += simValue[i]*simValue[i];
    normData += dataValue[i]*dataValue[i];
  }

  double denominator =  std::sqrt(normSim*normData);
  double norm_angle = 0;



  if(denominator == 0) norm_angle = 1; // 1, dissimilar
  else if (dot_product > denominator) norm_angle = 0; // can happen if angle is zero, due to rounding errors ?
  else norm_angle = std::acos(dot_product/denominator)/(3.14159);


  double sse_weight = 0.5;

  return (1-sse_weight)*norm_angle+ sse_weight*sse;
}


double mscs::
addCost(const std::vector< std::vector<double> > &y, const double t) 
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


  size_t num_cost_species = costList_.size();
  size_t num_cells = y.size();

  std::vector< double > simValue (num_cells);
  std::vector< double > dataValue (num_cells);


  double costAdd=0.;
  // Loop over species (k)
  for (size_t k=0 ; k<num_cost_species ; ++k ) {
    int j=costList_[k];
    // Loop over cells (i)
    for (size_t i=0 ; i<num_cells ; ++i ) {

      // set to elements from simulation result
      simValue[i] = y[i][j];

      // set with costTemplate elements
      dataValue[i] = costTemplate_[costTemplateIndex_][i][k];
    }
    double tempCost = getCost(simValue, dataValue); // angle between vectors, normalized
    costAdd += tempCost/(num_cost_species*costTemplate_.size());

    costIndividual_[costTemplateIndex_][k] += tempCost/(costTemplate_.size()); // error for species k at costTemplateIndex_ 

    } 
  

  //std::cerr << "cost add: " << costAdd << "\n";
  costTemplateIndex_++;
  cost_ += costAdd;
  return costAdd;
}


double KullbackLeibler::getCost(double simValue, double dataValue)
{

  //std::cerr << "data: " << dataValue << " sim: " << simValue << "\n";
  double cost;
  if(dataValue == 0) cost = 0; // ln(x)*x -> 0 when x -> 0
  else cost = std::log(dataValue/simValue)*dataValue;
  //std::cerr << "cost: " << cost << "\n";
  return cost;
}


double KullbackLeibler::
addCost(const std::vector< std::vector<double> > &y, const double t) 
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


  size_t num_cost_species = costList_.size();
  size_t num_cells = y.size();

  std::vector< double > simValue (num_cells*num_cost_species);
  std::vector< double > dataValue (num_cells*num_cost_species);
  double costAdd=0.;
  double simTot = 0;
  double dataTot= 0;
  // Loop over species (k)
  size_t vector_index=0;
  //std::cerr << "\n";

  // calculate normalization factors and store all data in 1D vectors for convenience
  for (size_t k=0 ; k<num_cost_species ; ++k ) {
    int cost_species_index=costList_[k];
    for (size_t i=0 ; i<num_cells ; ++i ) {
      simValue[vector_index] = y[i][cost_species_index];

      simTot += simValue[vector_index];
      dataValue[vector_index] = costTemplate_[costTemplateIndex_][i][k];
      dataTot += dataValue[vector_index];
      vector_index++;
    }
  }


    for (size_t k=0 ; k<simValue.size() ; ++k ) {

      double normSim = simValue[k]/simTot;
      double normData = dataValue[k]/dataTot;


      if(simTot < 0.0 or simValue[k] < 0.0){ // would mean negative probabilities?
	costAdd = 1.0/(costTemplate_.size());
	break;
      }
    

      costAdd += getCost(normSim, normData)/(costTemplate_.size()*simValue.size());

      if(getCost(normSim, normData) != getCost(normSim, normData)){
      std::cerr << "simTot: " << simTot;
      std::cerr << "simValue: " << simValue[k];
      std::cerr << " dataTot: " << dataTot << "\n";
      }
    
  
    }

  //std::cerr << "cost add: " << costAdd << "\n";
  costTemplateIndex_++;
  cost_ += costAdd;
  return costAdd;
}

