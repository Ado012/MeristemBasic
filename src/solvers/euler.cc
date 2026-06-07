//
// Filename     : euler.cc
// Description  : Simple Euler numerical solvers
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : March 2009
// Revision     : $Id:$
//
#include <cmath>
#include "euler.h"

Euler::Euler(Organism *O,std::ifstream &IN)
  :BaseSolver(O,IN)
{
  readParameterFile(IN);
}

void Euler::readParameterFile(std::ifstream &IN)
{
	//Read in the needed parameters
	IN >> startTime_;
	t_=startTime_;
	IN >> endTime_;

	IN >> printFlag_;
	IN >> numPrint_;

	IN >> h_;
}

void Euler::simulate() 
{
  // Check that h1 and endTime-startTime are > 0
  //
  if( !(h_>0. && (endTime_-startTime_)>0.) ) {
    std::cerr << "Euler::simulate() Wrong time borders or time step for "
	      << "simulation. No simulation performed.\n";
    return;
  }
  std::cerr << "Simulating using a simple Euler solver\n";

  // Initiate reactions for those where it is applicable
  //
  O_->reactionInitiate(startTime_,y_);

  // Introduce the neighborhood
  //
  if( O_->numNeighborhood() )
    O_->neighborhoodCreate(y_,t_);
  if( y_.size() && y_.size() != dydt_.size() ) {
    dydt_.resize( y_.size(),y_[0]);
  }
  assert( y_.size() == O_->numCompartment() && y_.size()==dydt_.size());
  
  // Initiate print times
  //////////////////////////////////////////////////////////////////////
  double tiny = 1e-10;
  double printTime=endTime_+tiny;
  double printDeltaTime=endTime_+2.*tiny;
  if( numPrint_<=0 )//No printing
    printFlag_=0;
  else if( numPrint_==1 ) {//Print last point (default)
  }
  else if( numPrint_==2 ) {//Print first/last point
    printTime=startTime_-tiny;
  } 
  else {//Print first/last points and spread the rest uniformly
    printTime=startTime_-tiny;
    printDeltaTime=(endTime_-startTime_)/double(numPrint_-1);
  }
  
  // Go
  //////////////////////////////////////////////////////////////////////
  t_=startTime_;
  numOk_ = numBad_ = 0;
  while( t_<endTime_ ) {
    if (debugFlag()) {
      yCopy_[debugCount()] = y_;
    } 
    //Update the derivatives
    O_->derivs(y_,dydt_);
    
    //Print if applicable 
    if( printFlag_ && t_ >= printTime ) {
      printTime += printDeltaTime;
      print();
    }

    //Update
    euler();
    numOk_++;
    
    // Make updates of rections that requires it
    O_->reactionUpdate(h_, t_, y_);

    //DISABLE 121525
    // Check for cell divisions, compartmental removal etc
    //O_->compartmentChange(y_,dydt_,t_);
    
    //Update neighborhood
    if( O_->numNeighborhood() )
      O_->neighborhoodUpdate(y_,t_);
    
    //update time variable
    if( (t_+h_)==t_ ) {
      std::cerr << "Euler::simulate() Step size too small.";
      exit(-1);
    }
    t_ += h_;
  }
  if( printFlag_ ) {
    //Update the derivatives
    O_->derivs(y_,dydt_);
    print();
  }
  std::cerr << "Simulation done.\n"; 
  return;
}

void Euler::euler()
{  
  size_t n = N();
  size_t m = M();
  for(size_t i=0 ; i<n ; ++i )
    for( size_t j=0 ; j<m ; ++j )
      y_[i][j] = y_[i][j]+h_*dydt_[i][j];
}

EulerEquilibrium::EulerEquilibrium(Organism *O,std::ifstream &IN)
  :BaseSolver(O,IN)
{
  readParameterFile(IN);
}

void EulerEquilibrium::readParameterFile(std::ifstream &IN)
{
	//Read in the needed parameters
	IN >> startTime_;
	t_=startTime_;
	IN >> endTime_;

	IN >> printFlag_;
	IN >> printInterval_;
	numPrint_ = 0;

	IN >> h_;
	IN >> threshold_;
}

void EulerEquilibrium::simulate() 
{
  // Check that h1 and endTime-startTime are > 0
  //
  if( !(h_>0. && (endTime_-startTime_)>0.) ) {
    std::cerr << "EulerEquilibrium::simulate() Wrong time borders or time step for "
	      << "simulation. No simulation performed.\n";
    return;
  }
  std::cerr << "Simulating using a simple EulerEquilibrium solver\n";

  // Initiate reactions for those where it is applicable
  //
  O_->reactionInitiate(startTime_,y_);

  // Introduce the neighborhood
  //
  if( O_->numNeighborhood() )
    O_->neighborhoodCreate(y_,t_);
  if( y_.size() && y_.size() != dydt_.size() ) {
    dydt_.resize( y_.size(),y_[0]);
  }
  assert( y_.size() == O_->numCompartment() && y_.size()==dydt_.size());
  
  // Initiate print times
  //////////////////////////////////////////////////////////////////////
  double tiny = 1e-10;
  double printTime;
  double printDeltaTime=printInterval_;
  if( printInterval_ == 0 )// Print only last point
    printTime=endTime_ + tiny;
  else {//Print first/last points and spread the rest uniformly
    printTime=startTime_-tiny;
  }

  // Go
  //////////////////////////////////////////////////////////////////////
  double maxDeriv = 1000.0;
  double tempDeriv = 1000.0;
  t_=startTime_;
  numOk_ = numBad_ = 0;
  while( t_<endTime_ && maxDeriv > threshold_) {
    if (debugFlag()) {
      yCopy_[debugCount()] = y_;
    } 
    //Update the derivatives
    O_->derivs(y_,dydt_);
    
    //Print if applicable 
    if( printFlag_ && t_ >= printTime ) {
      printTime += printDeltaTime;
      print();
      numPrint_ +=1;
    }

    //Update
    eulerEquilibrium();
    numOk_++;

    maxDeriv = 0;
    for (size_t i = 0; i < N(); i++){
      for (size_t j = 0; j < M(); j++) {
	tempDeriv = fabs(dydt_[i][j]);
	if(tempDeriv > maxDeriv){
	  maxDeriv = tempDeriv;	}
      }
    }
    // Make updates of rections that requires it
    O_->reactionUpdate(h_, t_, y_);
    
    //DISABLE 121525
    // Check for cell divisions, compartmental removal etc
    //O_->compartmentChange(y_,dydt_,t_);
    
    //Update neighborhood
    if( O_->numNeighborhood() )
      O_->neighborhoodUpdate(y_,t_);
    
    //update time variable
    if( (t_+h_)==t_ ) {
      std::cerr << "EulerEquilibrium::simulate() Step size too small.";
      exit(-1);
    }
    t_ += h_;
  }
  if( printFlag_ ) {
    //Update the derivatives
    O_->derivs(y_,dydt_);
    print();
  }
  std::cerr << "Simulation done.\n"; 
  return;
}

void EulerEquilibrium::eulerEquilibrium()
{  
  size_t n = N();
  size_t m = M();
  for(size_t i=0 ; i<n ; ++i )
    for( size_t j=0 ; j<m ; ++j )
      y_[i][j] = y_[i][j]+h_*dydt_[i][j];
}


EulerTruncated::EulerTruncated(Organism *O,std::ifstream &IN)
  :BaseSolver(O,IN)
{
  readParameterFile(IN);
}

void EulerTruncated::readParameterFile(std::ifstream &IN)
{
  //Read in the needed parameters
  IN >> startTime_;
  t_=startTime_;
  IN >> endTime_;
  
  IN >> printFlag_;
  IN >> numPrint_;
  
  IN >> h_;
  IN >> maxVal_;
}

void EulerTruncated::simulate() 
{
  //Check that h1 and endTime-startTime are > 0
  //////////////////////////////////////////////////////////////////////
  if( !(h_>0. && (endTime_-startTime_)>0.) ) {
    std::cerr << "EulerTruncated::simulate() Wrong time borders or time step for "
	      << "simulation. No simulation performed.\n";
    return;
  }
  std::cerr << "Simulating using a simple EulerTruncated solver\n";
  std::cerr << "WARNING" << std::endl 
	    << "WARNING  This solver should not be used in any circumstances" << std::endl
	    << "WARNING  where the result is used. It can only be used as " << std::endl
	    << "WARNING  an initial quick testing of a simulation." << std::endl
	    << "WARNING" << std::endl;

  // Initiate reactions for those where it is applicable
  //
  O_->reactionInitiate(startTime_,y_);

  //Introduce the neighborhood
  //////////////////////////////////////////////////////////////////////
  if( O_->numNeighborhood() )
    O_->neighborhoodCreate(y_,t_);
  if( y_.size() && y_.size() != dydt_.size() ) {
    dydt_.resize( y_.size(),y_[0]);
  }
  assert( y_.size() == O_->numCompartment() && y_.size()==dydt_.size());
  
  // Initiate print times
  //////////////////////////////////////////////////////////////////////
  double tiny = 1e-10;
  double printTime=endTime_+tiny;
  double printDeltaTime=endTime_+2.*tiny;
  if( numPrint_<=0 )//No printing
    printFlag_=0;
  else if( numPrint_==1 ) {//Print last point (default)
  }
  else if( numPrint_==2 ) {//Print first/last point
    printTime=startTime_-tiny;
  } 
  else {//Print first/last points and spread the rest uniformly
    printTime=startTime_-tiny;
    printDeltaTime=(endTime_-startTime_)/double(numPrint_-1);
  }
  
  // Go
  //////////////////////////////////////////////////////////////////////
  t_=startTime_;
  numOk_ = numBad_ = 0;
  while( t_<endTime_ ) {
    if (debugFlag()) {
      yCopy_[debugCount()] = y_;
    } 
    //Update the derivatives
    O_->derivs(y_,dydt_);
    
    //Print if applicable 
    if( printFlag_ && t_ >= printTime ) {
      printTime += printDeltaTime;
      print();
    }

    //Update
    eulerTrunc();    
    numOk_++;
    if (numInvalidValues()) {
      std::cerr << "Euler::simulate() Invalid numbers found in data or derivative matrix."
		<< std::endl;
      std::cerr << "Time: " << t_ << std::endl;
      //mySignal::myExit();
      if (debugFlag_) {
	printDebug(std::cerr);
      }
      exit(EXIT_FAILURE);
    }
    
    // Make updates of rections that requires it
    O_->reactionUpdate(h_, t_, y_);
    
    //DISABLE 121525
    // Check for cell divisions, compartmental removal etc
    //O_->compartmentChange(y_,dydt_,t_);
    
    //Update neighborhood
    if( O_->numNeighborhood() )
      O_->neighborhoodUpdate(y_,t_);
    
    //update time variable
    if( (t_+h_)==t_ ) {
      std::cerr << "Rk4::simulate() Step size too small.";
      exit(-1);
    }
    t_ += h_;
  }
  if( printFlag_ ) {
    //Update the derivatives
    O_->derivs(y_,dydt_);
    print();
  }
  std::cerr << "Simulation done.\n"; 
  return;
}

void EulerTruncated::eulerTrunc()
{  
  size_t n = N();
  size_t m = M();
  for(size_t i=0 ; i<n ; ++i )
    for( size_t j=0 ; j<m ; ++j ) {
      y_[i][j] = y_[i][j]+h_*dydt_[i][j];
      if (y_[i][j]<0.0) {
	y_[i][j]=0.0;
      }
      else if (std::isnan(y_[i][j]) || std::isinf(y_[i][j])) {
	y_[i][j] = maxVal_;
	if (std::isnan(dydt_[i][j]) || std::isinf(dydt_[i][j]))
	  dydt_[i][j] = 0.0;
      }

    }
}
