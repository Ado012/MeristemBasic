#include <vector>
#include <iostream>
#include "../organism.h"
#include "../solvers/baseSolver.h"
#include "../common/myConfig.h"
#include "../common/mySignal.h"
#include "../common/myRandom.h"
#include "../common/myFiles.h" //added in reorganization 042121

int main(int argc,char *argv[]) 
{
	myConfig::registerOption("init_output", 1);
	myConfig::registerOption("debug_output", 1);
	myConfig::registerOption("parameter_input", 1);
	myConfig::registerOption("verbose", 1);
	myConfig::registerOption("help", 0);

	std::string configFile(getenv("HOME"));
	configFile.append("/.organism");
	myConfig::initConfig(argc, argv, configFile);
	
	if (myConfig::getBooleanValue("help")) {
	  std::cerr << std::endl 
		    << "Usage: " << argv[0] << " modelFile initFile "
		    << "simulatorParaFile." << std::endl 
		    << std::endl;
	  std::cerr << "Possible additional flags are:" << std::endl;
	  std::cerr << "-init_output file - Set filename for output of"
		    << " final state in init file format." << std::endl;
	  std::cerr << "-debug_output file - Saves the last ten variable"
		    << " states before exiting." << std::endl;
	  std::cerr << "-parameter_input file - Read parameter values from file"
		    << " overriding those read in the model file." << std::endl;
	  std::cerr << "-verbose flag - Set flag for verbose (flag=1) or "
		    << "silent (0) output mode to stderr." << std::endl; 
	  exit(EXIT_FAILURE);
	} else if (myConfig::argc() != 4) {
	  std::cerr << "Wrong number of arguments given to " << argv[0] << std::endl
		    << "Type '" << argv[0] << " -help' for usage." << std::endl;
	  exit(EXIT_FAILURE);
	}
	std::string modelFile = myConfig::argv(1);
	std::string initFile = myConfig::argv(2);
	std::string simPara = myConfig::argv(3);
	
	int verboseFlag=0;
	std::string verboseString;
	verboseString = myConfig::getValue("verbose", 0);
	if( !verboseString.empty() ) {
	  verboseFlag = atoi( verboseString.c_str() );
	  if( verboseFlag != 0 || verboseFlag !=1 ) {
	    verboseFlag=0;
	    std::cerr << "Flag given to -verbose not recognized (0, 1 allowed)."
		      << " Setting it to zero (silent)." << std::endl;
	  }
	}
	//Randomize
	myRandom::Randomize();
	// Get current time (at start of program)
	myTimes::getTime();
	
	// Define the organism (model)
	Organism O(modelFile,initFile,verboseFlag);
	
	// Initiate parameter vector
	O.initiateParameter();
	
	// Overrride parameter values in model file if applicable 
	// Caveat: currently it assumes correct number of parameters and 
	// do not allow for comments in the file.
	std::string parameterValueFile;
	parameterValueFile = myConfig::getValue("parameter_input", 0);
	if( !parameterValueFile.empty() ) {
		std::istream *IN = myFiles::openFile(parameterValueFile);
		if (!IN) {
			std::cerr << "Warning: main() -"
								<< "Cannot open file for parameter reading ("
								<< parameterValueFile << ")." << std::endl;
		} 
		else {
			std::cerr << "Overriding model parameters with values from file "
								<< parameterValueFile << "." << std::endl; 
			double tmpDouble;
			for( size_t i=0; i<O.numParameter(); ++i ) {
				*IN >> tmpDouble;			
				O.setParameter(i,tmpDouble);
			}
			delete IN;
		}
	}
	
  // Get a solver from parameter file.
  BaseSolver *S = BaseSolver::getSolver(&O, simPara);
	
	// Print the model to standard error if applicable
	if( verboseFlag )
		O.printModel();
	
  // Add solver to signal handler.
  mySignal::addSolver(S);
	
  // Simulate with updates of the neighborhood
  std::cerr << "Start simulation.\n";
  S->getInit();
	S->simulate();
	if ( S->cost() != -1 ) {
		S->printCost();
	}
	// Print init if applicable
	std::string fileName;
	fileName = myConfig::getValue("init_output", 0);
	if(!fileName.empty()) {
		std::ofstream OUT(fileName.c_str());
		if (!OUT) {
			std::cerr << "Warning: main() -"
								<< "Cant open file for init output.\n";
		} else {
			S->printInit(OUT);
			OUT.close();
		}
	}

	//if( verboseFlag ) {
		O.printModel();
		for (size_t i=0; i<O.numParameter(); ++i)
			std::cerr << O.parameter(i) << " ";
		std::cerr << std::endl;
		//}
	// Print debug information if applicable
	fileName = myConfig::getValue("debug_output", 0);
	if(!fileName.empty()) {
		std::ofstream OUT(fileName.c_str());
		if (!OUT) {
			std::cerr << "Warning: main() -"
								<< "Cant open file for debug output.\n";
		} else {
			S->printDebug(OUT);
		}
	}

  // Delete the solver.
  delete S;
}
