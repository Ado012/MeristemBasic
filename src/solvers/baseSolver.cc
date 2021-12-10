//
// Filename     : baseSolver.cc
// Description  : Base class for solvers 
// Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
// Created      : August 2006
// Revision     : $Id: baseSolver.cc 602 2015-01-02 18:10:19Z henrik $
//
#include <cmath>
#include<ios>

#include "baseSolver.h"
#include "rungeKutta.h"
//DISABLED 042521
/*
#include "euler.h"
#include "implicit.h"
#include "gillespie.h"
#include "heunito.h"
*/

#include "../common/myFiles.h"
#include "../common/myConfig.h"

BaseSolver::BaseSolver()
{
	C_=0;
}

BaseSolver::BaseSolver(Organism *O,std::ifstream &IN)
{
  C_=0;
  setOrganism(O);
  getInit();

  //check debugging status
  std::string debugCheck = myConfig::getValue("debug_output", 0);
  if(!debugCheck.empty()) {
    std::cerr << "Performing simulation in debug-mode\n";
    debugFlag_ = true;
    if (debugFlag_) {
      int numCopy = 10;
      yCopy_.resize(numCopy);
    }
  }
  else 
    debugFlag_=false;
}

BaseSolver::~BaseSolver()
{
	
}

size_t BaseSolver::debugCount() const
{
  static size_t count = 0;
  if (debugFlag()) 
    return count++ % yCopy_.size(); 
  std::cerr << "Warning  BaseSolver::debugCount() should never be" 
	    << " called when debugFlag == " << debugFlag() << "\n"; 
  exit(-1);
}

BaseSolver *BaseSolver::getSolver(Organism *O, const std::string &file)
{
	std::istream *IN = myFiles::openFile(file);
	if (!IN) {
		std::cerr << "Solver::Solver() - "
				<< "Cannot open file " << file << std::endl;
		exit(EXIT_FAILURE);
	}
	BaseSolver* z = getSolver(O, (std::ifstream &) *IN);
	delete IN;
	return z;
	
}

BaseSolver *BaseSolver::getSolver(Organism *O, std::ifstream &IN)
{

	std::string idValue;
	IN >> idValue;

	BaseSolver *solver;
	// rungeKutta.h
    if (idValue == "RK5Adaptive")
      solver = new RK5Adaptive(O,(std::ifstream &) IN);
    else if (idValue == "RK5AdaptiveTemplate")
      solver = new RK5AdaptiveTemplate(O,(std::ifstream &) IN);
    else if (idValue == "RK5AdaptiveEquilibrium")
      solver = new RK5AdaptiveEquilibrium(O,(std::ifstream &) IN);
    else if (idValue == "RK4")
	  solver = new RK4(O,(std::ifstream &) IN);
    else if(idValue == "Dopr853")
      solver = new Dopr853(O,(std::ifstream &) IN);
//DELETION 042521
	else {
	  std::cerr << "BaseSolver::BaseSolver() (NYTT FELMEDDELANDE)- "
		    << "Unknown solver: " << idValue << std::endl;
	  delete &IN;
	  exit(EXIT_FAILURE);
	}
	return solver;
}

void BaseSolver::getInit()
{
  // Resize and define parameters
  size_t N = O_->numCompartment();
  size_t M = O_->numVariable();
  
  // Copy values from the Organism into the y matrix and resize dydt
  
  y_.resize(N);
  
  dydt_.resize(N);
  for (size_t i = 0; i < N; i++) {
    y_[i].resize(M);
    dydt_[i].resize(M);
    for (size_t j = 0; j < M; j++)
      y_[i][j] = O_->compartment(i).variable(j); 
  }
  
}

void BaseSolver::readParameterFile(std::ifstream &IN)
{
  std::cerr << "BaseSolver::readParameterFile(std::ifstream &IN) Not Defined for the subclass of BaseSolver. Empty Solver created without type, neither simulate function.\n";
  exit(-1);
}

void BaseSolver::simulate(void)
{
  std::cerr << "BaseSolver::simulate(void) Not Defined for the subclass of BaseSolver. Empty Solver created without type, neither simulate function.\n";
  exit(-1);
}

int BaseSolver::numInvalidValues()
{
  int falseCount=0;
  
  for (size_t i=0; i<N(); ++i)
    for (size_t j=0; j<M(); ++j) {
      if (std::isnan(y_[i][j]) || std::isinf(y_[i][j]))
	++falseCount;
      if (std::isnan(dydt_[i][j]) || std::isinf(dydt_[i][j]))
	++falseCount;
    }
  return falseCount;
}

void BaseSolver::divideTemplate(int i) 
{  
  int NN=N();
  int Nval=NN+1;
  int Mval=M();
  
  //Resize organism 
  O_->addCompartment( O_->compartment(i) );
  O_->compartment(NN).setIndex(NN);
  O_->addDivision(i,t_);
  
  //Resize y_, dydt_, previousValue_ and previousTime_
  y_.resize( Nval );
  dydt_.resize( Nval );
  previousValue_.resize( Nval );
  previousTime_.resize( Nval );
  y_[NN].resize( Mval );
  dydt_[NN].resize( Mval );
  previousValue_[NN].resize( Mval);
  for( int j=0 ; j<Mval ; j++ ) {
    y_[NN][j] = y_[i][j];
    
  } 
  
  //Adjust previousValue_ and previousTime_
  previousTime_[i]=previousTime_[NN]=t_;
  for( int j=0 ; j<Mval ; j++ )
    if( !simulationFlag_[j] )
      previousValue_[i][j] = previousValue_[NN][j] = y_[i][j];
  
  //Calculate new derivatives
  double deltaTime=futureTime_-t_;
  if( deltaTime>0. )
    for( int j=0 ; j<Mval ; j++ )
      if( !simulationFlag_[j] ) {
				dydt_[i][j] = (futureValue_[i][j] - y_[i][j])/deltaTime;
				dydt_[NN][j] = (futureValue_[NN][j] - y_[NN][j])/deltaTime;;
      }  
}

void BaseSolver::initiateValuesTemplate(std::ifstream &IN,
	std::vector<int> &divIndex,
	std::vector<double> &divTime)
{
  size_t N1,M1,N2,M2,numDiv;
  double time1,time2;
  
  IN >> time1;
  IN >> N1;
  IN >> M1;
  
  if( time1>t_ ) {
    //Create a bogus point at time before start time
    double prevTimeTmp=t_-1.;//Anything less than t would do
    time2=time1;
    N2=N1;
    M2=M1;
    previousTime_.resize(N1);
    for ( size_t i=0 ; i<N1 ; i++ )
      previousTime_[i]=prevTimeTmp;
    futureTime_=time2;
    previousValue_.resize(N1);
    futureValue_.resize(N2);
    for (size_t i=0 ; i<N1 ; i++ ) {
      previousValue_[i].resize(M1);
      futureValue_[i].resize(M2);
      for (size_t j=0 ; j<M1 ; j++ ) {
        IN >> previousValue_[i][j];
        futureValue_[i][j]= previousValue_[i][j];
      }
    }
  }
  else {
    //Complete first time point (previousValue)
    previousTime_.resize(N1);
    for ( size_t i=0 ; i<N1 ; i++ )
      previousTime_[i]=time1;
    previousValue_.resize(N1);
    for (size_t i=0 ; i<N1 ; i++ ) {
      previousValue_[i].resize(M1);
      for (size_t j=0 ; j<M1 ; j++ ) {
        IN >> previousValue_[i][j];
      }
    }
    int count=0;
    do {
      if( count ) {
        //Copy the future time point into the previous one
        previousValue_.resize( futureValue_.size() );
        previousTime_.resize( futureValue_.size() );
        for( size_t i=0 ; i<previousValue_.size() ; i++ ) {
          previousTime_[i]=futureTime_;
          previousValue_[i].resize( futureValue_[i].size() );
          for( size_t j=0 ; j<previousValue_[i].size() ; j++ )
            previousValue_[i][j] = futureValue_[i][j];
        }
        count++;
      }
      IN >> numDiv;
      divIndex.resize( numDiv );
      divTime.resize( numDiv );
      for (size_t i=0 ; i<numDiv ; i++ ) {
        IN >> divIndex[i];
        IN >> divTime[i];
      }
      //Read next time point (futureValue)
      IN >> time2;
      IN >> N2;
      IN >> M2;
      futureTime_=time2;
      futureValue_.resize(N2);
      for (size_t i=0 ; i<N2 ; i++ ) {
        futureValue_[i].resize(M2);
        for (size_t j=0 ; j<M2 ; j++ ) {
          IN >> futureValue_[i][j];
        }
      }
    } while( time2<=t_ && IN );
    
    if( time2<=t_ ) {
      //All time points before "start time" t (end of file reached)
      //Create a bogus point after max time to be able to finish simulation
      
      previousTime_.resize(N1);
      for (size_t i=0 ; i<N1 ; i++ )
        previousTime_[i]=time1;
      previousValue_.resize(N1);
      for (size_t i=0 ; i<N1 ; i++ ) {
        previousValue_[i].resize(M1);
        for (size_t j=0 ; j<M1 ; j++ ) {
          IN >> previousValue_[i][j];
        }
      }
    }
  }
  //Do some checks for consistancy
  if( N1 != N() ) {
    std::cerr << "Simulator::initiateValuesFile : Wrong number of "
							<< "data points in input file " << inputTemplateFile() << "\n";
    exit(-1);
  }
  if( M1 != M() || M2 != M() ) {
    std::cerr << "Simulator::initiateValuesFile : Wrong number of "
							<< "data columns in input file " << inputTemplateFile() << "\n";
    exit(-1);
  }
  if( N1+numDiv != N2 ) {
    std::cerr << "Simulator::initiateValuesFile : Wrong number of "
							<< "divisions in input file " << inputTemplateFile() << "("
							<< N1+numDiv << " " << N2 << ")\n";
    exit(-1);
  }
  //Calculate derivatives and present values for columns that should be read
  for (size_t j=0 ; j<M() ; ++j)
    if( !simulationFlag(j) )
      for (size_t i=0 ; i<N() ; ++i) {
        dydt_[i][j] = ( futureValue_[i][j] - previousValue_[i][j] )
          / (futureTime_ - previousTime_[i]);
        y_[i][j] = previousValue_[i][j] + (t_-previousTime_[i])*dydt_[i][j];
      }
  //std::cerr << "Initiated from file:\n";
  //std::cerr << previousValue_.size() << " " << M1 << " " << previousTime_[0] 
  //    <<"\n";
  //std::cerr << futureValue_.size() << " " << M2 << " " << futureTime_ 
  //<< "\n";
}

void BaseSolver::updateValuesTemplate( std::ifstream &IN,double tmax,
	std::vector<int> &divIndex,
	std::vector<double> &divTime) 
{
  static int updateFlag=1;
	
  if( !updateFlag ) {
    //Don't need to update since end of file was reached last update.
    return;
  }
  // Copy the old futureValues into previousValues
  //////////////////////////////////////////////////////////////////////
  if( previousValue_.size() != futureValue_.size() ) {
    //Following all division these should have the same size right now.
    std::cerr << "Simulator::updateValuesFile() : Inconsistancy in "
							<< "number of cells at shift.\n";
    exit(-1);
  }
  for (size_t i=0 ; i<previousValue_.size() ; i++ ) {
    previousTime_[i]=futureTime_;
    for (size_t j=0 ; j<previousValue_[i].size() ; j++ )
      previousValue_[i][j] = futureValue_[i][j];
  }
  size_t N2,M2,numDiv;
  double time2;
  // Try to read new division data
  //////////////////////////////////////////////////////////////////////
  IN >> numDiv;
  if( IN ) {//More data to be read (is this fool proof?)
    divIndex.resize( numDiv );
    divTime.resize( numDiv );
    for (size_t i=0 ; i<numDiv ; ++i) {
      IN >> divIndex[i];
      IN >> divTime[i];
    }
    //Read new values into futureValues
    IN >> time2;
    IN >> N2;
    IN >> M2;
    futureTime_=time2;
    futureValue_.resize(N2);
    for (size_t i=0 ; i<N2 ; i++ ) {
      futureValue_[i].resize(M2);
      for (size_t j=0 ; j<M2 ; j++ ) {
        IN >> futureValue_[i][j];
      }
    }
  }
  else { //No more values to read. Assume constant variable values
    numDiv=0;
    divIndex.resize( numDiv );
    divTime.resize( numDiv );
    //set future time to time above timemax and deflag updateFlag
    futureTime_=tmax+1.;//Anything above tmax would work
    updateFlag=0;
  }
  size_t N1 = previousValue_.size();
  //Do some consistency checks
  if( N1 != N() ) {
    std::cerr << "Simulator::updateValuesFile : Wrong number of "
							<< "data points in input file " << inputTemplateFile() << "\n";
    exit(-1);
  }
  if( N1+numDiv != N2 ) {
    std::cerr << "Simulator::updateValuesFile : Wrong number of "
							<< "divisions in input file " << inputTemplateFile() << "("
							<< N1+numDiv << " " << N2 << ")\n";
    exit(-1);
  }
  //Calculate derivatives and present values for columns that should be read
  for (size_t j=0 ; j<M() ; j++ )
    if( !simulationFlag(j) )
      for (size_t i=0 ; i<N() ; i++ ) {
        dydt_[i][j] = ( futureValue_[i][j] - previousValue_[i][j] )
          / (futureTime_ - previousTime_[i]);
        y_[i][j] = previousValue_[i][j] + (t_-previousTime_[i])*dydt_[i][j];
      }
  //size_t M1=M();
  //std::cerr << "Updated from file:\n";
  //std::cerr << N1 << " " << M1 << " " << previousTime_[0] << "\n";
  //std::cerr << N2 << " " << M2 << " " << futureTime_ << "\n";
}

void BaseSolver::allValuesTemplate() 
{
  for (size_t j=0 ; j<M() ; ++j)
    if( !simulationFlag(j) )
      for (size_t i=0 ; i<N() ; ++i)
        y_[i][j] = previousValue_[i][j] + (t_-previousTime_[i])*dydt_[i][j];
}

double BaseSolver::valueTemplate( int i, int j) 
{
  if( !simulationFlag(j) )
    return previousValue_[i][j] + (t_-previousTime_[i])*dydt_[i][j];
  std::cerr << "Simulator::valueFile : Asking for value which "
						<< "should be simulated.\n";
  exit(-1);
}

void BaseSolver::setSimulationFlag(int verbose) 
{
	simulationFlag_.resize( M() );
	
	for (size_t i=0 ; i<M() ; ++i) {
		simulationFlag_[i] = 1;
	}
	for (size_t i=0; i<simulationIndex_.size(); ++i)
		simulationFlag_[simulationIndex_[i]] = 0;
	if (verbose) {
		std::cerr << "Simulation flag for the " << M() << " variables.\n";
		for (size_t i=0 ; i<M() ; ++i) {
			std::cerr << simulationFlag_[i] << " ";
		}
		std::cerr << "\n";
	}
}

void BaseSolver::print(std::ostream &os) 
{  
  static int tCount=0;
  static int NOld=0,okOld=0,badOld=0;
  static double tOld=0.0;
	double time=myTimes::getDiffTime();
	if (!(tCount%40)) {
		std::cerr << "#Count time deltaT NumCompart deltaNc NumVar ";
		std::cerr << "NumSteps deltaNS NumFail deltaNF ";
		//printStatisticsHeader();
		std::cerr << "userTime" << std::endl;
	}
	std::cerr << tCount << "\t" << t_ << "\t" << t_-tOld << "\t" 
		  << N() << "\t" << static_cast<int>(N())-static_cast<int>(NOld) << "\t"
		  << M() << "\t";
	std::cerr << numOk_ << "\t" << numOk_-okOld << "\t" << numBad_ << "\t" 						 
		  << "\t" << numBad_-badOld << "\t";
	//printStatistics();
	std::cerr << time << std::endl;
	
	tOld = t_;
	NOld = N();
	okOld = numOk_;
	badOld = numBad_;
	
	// For most home-made plotting
	// 10 is kept for old time sake...
	if( printFlag_==1 || printFlag_==10 ) {
	  int overlapFlag=0;
	  if( tCount==0 )
	    os << numPrint_ << "\n";
	  os << N() << " " << M()+1+overlapFlag << " " << t_ << "\n";
	  for( size_t i=0 ; i<N() ; i++ ) {
	    for( size_t j=0 ; j<M() ; j++ ) {
	      os << y_[i][j] << " ";
	    }
	    os << O_->compartment(i).numNeighbor();
	    if( overlapFlag ) {
	      //Print summed overlap larger than k_fac
	      double k_fac = 0.75, overlap = 0.0;
	      for( size_t k=0 ; k<O_->compartment(i).numNeighbor() ; k++ ) {
		size_t ii = O_->compartment(i).neighbor(k);
		double r1 = y_[i][O_->compartment(i).numDimension()]; 
		double r2 = y_[ii][O_->compartment(i).numDimension()]; 
		double dist=0.0;
		for( size_t d=0 ; d<O_->compartment(i).numDimension() ; d++ )
		  dist += (y_[i][d]-y_[ii][d])*(y_[i][d]-y_[ii][d]);
		dist = std::sqrt( dist);
		if( dist < k_fac*(r1+r2) )
		  overlap += k_fac*(r1+r2)-dist;
	      }
	      os << " " << overlap/O_->compartment(i).numNeighbor();
	    }
	    os << "\n";
	  }
	  os << "\n";
	}
	else if (printFlag_ == 2) { // For gnuplot to be able to follow individual cells including derivatives
	  for (size_t i = 0; i < N(); ++i) {
	    os << tCount << " " << t_ << " " << i << " ";
	    os << O_->compartment(i).numNeighbor() << " ";	  
	    for (size_t j = 0; j < M(); ++j)
	      os << y_[i][j] << " ";      
	    for (size_t j = 0; j < M(); ++j)
	      os << dydt_[i][j] << " ";      
	    os << std::endl;
	  }
	  os << std::endl;
	}
	else if( printFlag_==3 ) {// Cost template style
	  static size_t numDivOld=0;
	  if( tCount ) {
	    os << O_->numDivision()-numDivOld << "\n";;
	    for( size_t i=numDivOld ; i<O_->numDivision() ; i++ )
				os << O_->divisionIndex(i) << " " << O_->divisionTime(i) << "\n";
      os << "\n";
      numDivOld = O_->numDivision();
    }    
    os << t_ << " " << N() << " " << M() << "\n";
    for( size_t i=0 ; i<N() ; i++ ) {
      for( size_t j=0 ; j<M() ; j++ ) {
				os << y_[i][j] << " ";
      }
      os << "\n";
    }
    os << "\n";
  }
  else if (printFlag_==4) {// Init style
    printInit(os);
  }
  else if (printFlag_==5) { //  Output formats 1 and 2, directed to different files
    int overlapFlag=0;
    if( tCount==0 )
      os << numPrint_ << "\n";
    os << N() << " " << M()+1+overlapFlag << " " << t_ << "\n";
    for( size_t i=0 ; i<N() ; i++ ) {
      for( size_t j=0 ; j<M() ; j++ ) {
	os << y_[i][j] << " ";
      }
      os << O_->compartment(i).numNeighbor();
      if( overlapFlag ) {
	//Print summed overlap larger than k_fac
	double k_fac = 0.75, overlap = 0.0;
	for( size_t k=0 ; k<O_->compartment(i).numNeighbor() ; k++ ) {
	  size_t ii = O_->compartment(i).neighbor(k);
	  double r1 = y_[i][O_->compartment(i).numDimension()]; 
	  double r2 = y_[ii][O_->compartment(i).numDimension()]; 
	  double dist=0.0;
	  for( size_t d=0 ; d<O_->compartment(i).numDimension() ; d++ )
	    dist += (y_[i][d]-y_[ii][d])*(y_[i][d]-y_[ii][d]);
	  dist = std::sqrt( dist);
	  if( dist < k_fac*(r1+r2) )
	    overlap += k_fac*(r1+r2)-dist;
	}
	os << " " << overlap/O_->compartment(i).numNeighbor();
      }
      os << "\n";
    }
    os << "\n";
    // Printing gnuplot style in file organism.gdata
    std::ofstream of;
    if (tCount==0) {
      of.open("organism.gdata");
    }
    else {
      of.open("organism.gdata",std::ios_base::app | std::ios_base::out);
    }
    for (size_t i = 0; i < N(); ++i) {
      of << tCount << " " << t_ << " " << i << " ";
      of << O_->compartment(i).numNeighbor() << " ";	  
      for (size_t j = 0; j < M(); ++j)
	of << y_[i][j] << " ";      
      for (size_t j = 0; j < M(); ++j)
	of << dydt_[i][j] << " ";      
      of << std::endl;
    }
    of << std::endl;
  }
  else if (printFlag_==101) {// For plotting neighborvariables
    for (size_t i=0; i<N(); ++i) {
			for (size_t k=0; k<O_->compartment(i).numNeighbor(); ++k) {
				size_t j = O_->compartment(i).neighbor(k);
				assert(O_->compartment(i).numNeighborVariable()>2);
				double f = O_->compartment(i).neighborVariable(2,k);
				std::vector<double> d(2);
				d[0] = y_[j][0]-y_[i][0]; 
				d[1] = y_[j][1]-y_[i][1]; 
				double norm = 1.0/(d[0]*d[0]+d[1]*d[1]);
				d[0] *= norm;
				d[1] *= norm;
				os << tCount << " 0 " << i << " " << j << " " 
					 << y_[i][0] << " " << y_[i][1] << std::endl
					 << tCount << " 1 " << i << " " << j << " " 
					 << y_[i][0]+f*d[0] << " " << y_[i][1]+f*d[1] << std::endl << std::endl << std::endl;
			}
		}
	}
	else if (printFlag_==102) {// For plotting parameters and specific data variables for simplified root
	  size_t lIndex = 1;
	  size_t gradIndex=11;//9
	  size_t allIndex=12; //10
	  size_t c1 = 1;
	  size_t c2 = 7;//4
	  size_t c3 = 10;//7
	  for (size_t i=0; i<O_->numParameter(); ++i)
	    os << O_->parameter(i) << " ";
	  os << y_[c1][lIndex] << " " << y_[c2][lIndex] << " " << y_[c3][lIndex] << " "
	     << y_[c1][gradIndex] << " " << y_[c2][gradIndex] << " " << y_[c3][gradIndex] << " "
	     << y_[c1][allIndex] << " " 
	     << y_[c1][gradIndex]-y_[c2][gradIndex] << " " 
	     << (y_[c1][gradIndex]-y_[c2][gradIndex])/(y_[c1][gradIndex]+y_[c2][gradIndex]) << " "
	     << cost();
	  os << std::endl;
	}	
	//
	// For vtk format (only works for 3D spheres and a single time point)
	else if (printFlag_==103) {
		os << "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\">" << std::endl
		<< "<UnstructuredGrid>" << std::endl
		<< "<Piece  NumberOfPoints=\"" << N() << "\" NumberOfCells=\"1\">" << std::endl
		<< "<Points>" << std::endl
		<< "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
		for (size_t i=0; i<N(); ++i) {
			for (size_t j=0; j<3; ++j) {
				os << y_[i][j] << " ";
			}
			os << std::endl;
		}
		os << "</DataArray>" << std::endl << "</Points>" << std::endl;

		os << "<Cells>" << std::endl
		<< "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
		for (size_t i=0; i<N(); ++i) {
			os << i << " ";
		}
		os << std::endl;
		os << "</DataArray>" << std::endl
		<< "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl
		<< N() << std::endl << "</DataArray>" << std::endl;
		os << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl
		<< "2" << std::endl << "</DataArray>" << std::endl << "</Cells>" << std::endl;

		os << "<PointData>" << std::endl
		<< "<DataArray  type=\"Float64\"  Name=\"radius\"  format=\"ascii\">" << std::endl;
		for (size_t i=0; i<N(); ++i) {
			os << y_[i][3] << " ";
		}
		os << "</DataArray>" << std::endl;

		for (size_t j=4; j<M(); ++j) {
			os << "<DataArray  type=\"Float64\"  Name=\"X" << j << "\" format=\"ascii\">" 
				<< std::endl;
			for (size_t i=0; i<N(); ++i) {
				os << y_[i][j] << " ";
			}
			os << "</DataArray>" << std::endl;
		}
		os << "</PointData>" << std::endl << "</Piece>" << std::endl 
		<< "</UnstructuredGrid>" << std::endl << "</VTKFile>" << std::endl;

		std::cerr << "PrintFlag=103, vtk printing, works only for 3D spheres and a single time point."
			<< std::endl << "Exiting after first print point." << std::endl;
		exit(0);    
	}  
	else if (printFlag_==111) { //Printing added variable values and number of cells (for averages)
		std::vector<double> sum(M());
		for (size_t i=0; i<N(); ++i) {
			for (size_t j=0; j<M(); ++j) {
				sum[j] += y_[i][j];
			}
		}
		for (size_t j=0; j<M(); ++j)
			os << sum[j] << " ";
		os << N() << " " << t_ << std::endl;
	}
	else if( printFlag_==121 ) {// For most home-made plotting
	  int overlapFlag=0;
	  if( tCount==0 )
	    os << numPrint_ << "\n";
	  os << N() << " " << M()+2+overlapFlag << " " << t_ << "\n";
	  for( size_t i=0 ; i<N() ; i++ ) {
	      os << y_[i][0] << " 0.0 ";
	    for( size_t j=1 ; j<M() ; j++ ) {
	      os << y_[i][j] << " ";
	    }
	    os << O_->compartment(i).numNeighbor();
	    if( overlapFlag ) {
	      //Print summed overlap larger than k_fac
	      double k_fac = 0.75, overlap = 0.0;
	      for( size_t k=0 ; k<O_->compartment(i).numNeighbor() ; k++ ) {
		size_t ii = O_->compartment(i).neighbor(k);
		double r1 = y_[i][O_->compartment(i).numDimension()]; 
		double r2 = y_[ii][O_->compartment(i).numDimension()]; 
		double dist=0.0;
		for( size_t d=0 ; d<O_->compartment(i).numDimension() ; d++ )
		  dist += (y_[i][d]-y_[ii][d])*(y_[i][d]-y_[ii][d]);
		dist = std::sqrt( dist);
		if( dist < k_fac*(r1+r2) )
		  overlap += k_fac*(r1+r2)-dist;
	      }
	      os << " " << overlap/O_->compartment(i).numNeighbor();
	    }
	    os << "\n";
	  }
	  os << "\n";
	}
	else {
	  std::cerr << "BaseSolver::print() Wrong printFlag value (" << printFlag_
		    << "). Allowed values are:\n";
	  std::cerr << "1 - For most home-made plotting, Newman included." << std::endl;
	  std::cerr << "2 - For gnuplot including time and cell indices and derivatives on each row." << std::endl;
	  std::cerr << "3 - For cost template format." << std::endl;
	  std::cerr << "4 - For init format." << std::endl;
	  std::cerr << std::endl << "Note that the numbers were changed and obselete versions removed"
		    << " at revision 398 (2008-11-20) and that more temporary versions do exist."
		    << " (1,2,3,4) correspond to (10,7,3,4) before the change."
		    << std::endl;
	}
  tCount++;
}

void BaseSolver::printInit(std::ostream &os) const
{
  os << N() << " " << M() << "\n";
  for (size_t i=0; i<N(); i++) {
    for (size_t j=0; j<M(); j++) {
      os << y_[i][j] << " ";
    }
    os << "\n";
  }
  os << "\n";
}

void BaseSolver::printDebug(std::ostream &os) const
{
  os << yCopy_.size() << "\n";
  size_t startElement = debugCount();
  for (size_t c=startElement; c<startElement+yCopy_.size(); ++c) {
    size_t n = c%yCopy_.size();
    os << yCopy_[n].size() << " " << M() << " " << c-startElement << "\n";
    for (size_t i=0; i<yCopy_[n].size(); i++) {
      for (size_t j=0; j<M(); j++) {
	os << yCopy_[n][i][j] << " ";
      }
      os << "\n";
    }
    os << "\n\n";
  }
}

double BaseSolver::maxDerivative() {  
  double max=0.0,val;
  for( size_t i=0 ; i<N() ; i++ ) {
    for( size_t j=0 ; j<M() ; j++ ) {
      if( y_[i][j]!=0 )
	val = fabs( dydt_[i][j]/y_[i][j] );
      else
	val = fabs( dydt_[i][j] );
      
      if( val>max ) 
	max = val;
    }
  }
  return max;
}

