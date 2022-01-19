
// Filename     : extendedMeristemReactions.cc
// Description  : Reactions for extended meristem sims
// Author(s)    : Al Do (ado012@ucr.edu)
// Created      : June 2017
// Revision     :
//

#include"baseReaction.h"
#include"extendedMeristemReactions.h"
#include"extendedMeristemReactionsHelperFunctions.h"
#include"../organism.h"
#include <random> //ADDITION 051517 Currently unusable due to problems with C11
#include <fstream>//ADDITION 051517
#include <stdlib.h>//ADDITION 051517
#include <string> //ADDITION 100218


Spatial::Spatial(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=7 ) {
		std::cerr << "Spatial::Spatial() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 3 ) {
		std::cerr << "Spatial::Spatial() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("spatial");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	  tmp[0] = "Co";
  tmp[1] = "xo";
  tmp[2] = "yo";
   tmp[3] = "zo";
  tmp[4] = "sig1";
  tmp[5] = "sig2";
  tmp[6] = "sig3";
	setParameterId( tmp );
}

void Spatial::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  


  double Co=parameter(0);
  double xo=parameter(1);
  double yo=parameter(2);
  double zo=parameter(3);
  double sig1=parameter(4);
  double sig2=parameter(5); 
  double sig3=parameter(6); 
  double xvar=y[compartment.index()][variableIndex(0,0)];
  double yvar=y[compartment.index()][variableIndex(0,1)];
  double zvar=y[compartment.index()][variableIndex(0,2)];
    double offsetval=0;
   //std::cerr <<"parameterlist"<<Co<<" "<< xo<< " " << yo<< " " << sig1<< " " <<sig2<< " " << xvar << " " <<yvar <<std::endl; 
    double b = -1*((std::pow(xvar-xo,2)/sig1)+(std::pow(yvar-yo,2)/sig2)+(std::pow(zvar-zo,2)/sig3));
   
   
   double c =  std::pow(2.71828,b); //seems to cause errors need to look into solving this for accuracy AD121817
  // double c = b*2.71828;//temp standin for spatial
       
//std::cerr  << "c" << c << std::endl;
   

   //check to ensure proper parentheses test stepper
    // double contribution = 1*std::pow(c,b);
     double contribution=c;
   


if (((xvar*xvar)+(yvar*yvar)+std::pow(zvar-zo,2)) <=36)
 {
 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
 // std::cerr  << "contribution " << contribution << std::endl;
  //std::cerr  << "wRNA " << dydt[compartment.index()][varIndex]  << std::endl;
  }
  }
  
  
  
}


////////////////////////////////////////////////////

Hill_Weitao::Hill_Weitao(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=9 ) {
		std::cerr << "Hill_Weitao::Hill_Weitao() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}//messed with error messages...check later AD102017
	if( indValue.size() != 1 || indValue[0].size() != 4 ) {
		std::cerr << "Hill_Weitao::Hill_Weitao() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("hill_Weitao");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	  tmp[0] = "Wp";
  tmp[1] = "xo";
  tmp[2] = "yo";
  tmp[3] = "zo";
  tmp[4] = "sig1";
  tmp[5] = "sig2";
   tmp[6] = "sig3";
  tmp[7] = "Kck";
  tmp[8] = "n";

	setParameterId( tmp );
	
}

void Hill_Weitao::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  


    double offsetval=0;
  double Wp=parameter(0);
  double xo=parameter(1);
  double yo=parameter(2);
  double zo=parameter(3);
  double sig1=parameter(4);
  double sig2=parameter(5); 
   double sig3=parameter(6); 
     double khill=parameter(7);
  double nhill=parameter(8);
  double xvar=y[compartment.index()][variableIndex(0,0)];
  
  double yvar=y[compartment.index()][variableIndex(0,1)];
  double zvar=y[compartment.index()][variableIndex(0,2)];
 
  double hillvar=y[compartment.index()][variableIndex(0,3)];
    
   //std::cerr <<"parameterlist"<<Co<<" "<< xo<< " " << yo<< " " << sig1<< " " <<sig2<< " " << xvar << " " <<yvar <<std::endl; 
    //stepper test
             //-10*((std::pow(xvar-xo,2))/sig1)+((std::pow(yvar-yo,2))/sig2)+((std::pow(zvar-zo,2))/sig3);
  double b = -1*((std::pow(xvar-xo,2)/sig1)+(std::pow(yvar-yo,2)/sig2)+(std::pow(zvar-zo,2)/sig3));
     double c = Wp*std::pow(2.71828,b);

 
 //std::cerr  << "test" << std::endl;
 
  double contribution = c
      / ( 1+std::pow((hillvar/khill),nhill));

 
 

  
if (((xvar*xvar)+(yvar*yvar)+std::pow(zvar-zo,2)) <=36)
 {
 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
  }
  }
  
  //std::cerr  << "wRNA " << y[compartment.index()][varIndex]  << std::endl;

  
}



////////////////////////////////////////////////////

Hill_Weitao2::Hill_Weitao2(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=11 ) {
		std::cerr << "Hill_Weitao2::Hill_Weitao2() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 6 ) {
		std::cerr << "Hill_Weitao2::Hill_Weitao2() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("hill_Weitao2");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	  tmp[0] = "dparam";
  tmp[1] = "khill";
  tmp[2] = "n2";
  tmp[3] = "khill2";
  tmp[4] = "n7";
  tmp[5] = "khill3";
  tmp[6] = "n6";
  tmp[7] = "WB";
  tmp[8] = "EL";
  tmp[9] = "n8";
  tmp[10]= "Clv3stabparam";
  
  
	setParameterId( tmp );
}

void Hill_Weitao2::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
    double offsetval=0;
  double dparam=parameter(0);
  double khill=parameter(1);
   double n2=parameter(2);
    double khill2=parameter(3);
 double n7=parameter(4);
  double khill3=parameter(5);
   double n6=parameter(6);
 double WB=parameter(7);
 double EL=parameter(8);
  double n8=parameter(9);
 double clv3stabparam=parameter(10);
 double distance;
 
  double hillvar=y[compartment.index()][variableIndex(0,0)];
    double hillvar2=y[compartment.index()][variableIndex(0,1)];
double hillvar3=y[compartment.index()][variableIndex(0,2)];
    
      double xvar=y[compartment.index()][variableIndex(0,3)];
  
  double yvar=y[compartment.index()][variableIndex(0,4)];
  double zvar=y[compartment.index()][variableIndex(0,5)];
    
    
   //std::cerr <<"parameterlist"<<dparam<<" "<< khill<< " " << nhill<< " " << nhill<< " " <<khill2<< " " << nhill2 << " " <<
   //" "<< WB <<" "<< EL <<" "<< hillvar << " " << hillvar2<< std::endl; 
    double c = (0.5*dparam+(0.5*dparam))/(1+std::pow(WB/0.5,-1*n8));

distance=sqrt(xvar*xvar+yvar*yvar+zvar*zvar);

if (distance<6)
distance=1;


double d= ( 1+std::pow((hillvar*WB/khill),n2));
double e= ( 1+std::pow((hillvar2/khill2),n7));
double f= ( 1+std::pow((hillvar3/khill3),n6));
double g= distance*clv3stabparam*y[compartment.index()][varIndex];



double contribution = -1*(c/ (d*e*f))*y[compartment.index()][varIndex];
//first constant temp was adjusted from 1 to 40




 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
  }

}


////////////////////////////////////////////////////

Hill_Weitao3::Hill_Weitao3(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=10 ) {
		std::cerr << "Hill_Weitao3::Hill_Weitao3() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 3 ) {
		std::cerr << "Hill_Weitao3::Hill_Weitao3() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("hill_Weitao3");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	  tmp[0] = "dparam";
  tmp[1] = "khill";
  tmp[2] = "nhill";
  tmp[3] = "khill2";
  tmp[4] = "nhill2";
  tmp[5] = "khill3";
  tmp[6] = "nhill3";
  tmp[7] = "WB";
  tmp[8] = "EL";
  tmp[9] = "nhill4";
	setParameterId( tmp );
}

void Hill_Weitao3::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
    double offsetval=0;
  double dparam=parameter(0);
  double khill=parameter(1);
   double nhill=parameter(2);
    double khill2=parameter(3);
 double nhill2=parameter(4);
     double khill3=parameter(5);
 double nhill3=parameter(6);
 double WB=parameter(7);
 double EL=parameter(8);
 double nhill4=parameter(9);
  double hillvar=y[compartment.index()][variableIndex(0,0)];
    double hillvar2=y[compartment.index()][variableIndex(0,1)];
    double hillvar3=y[compartment.index()][variableIndex(0,2)];
    
   //std::cerr <<"parameterlist"<<dparam<<" "<< khill<< " " << nhill<< " " << nhill<< " " <<khill2<< " " << nhill2 << " " <<
   //" "<< WB <<" "<< EL <<" "<< hillvar << " " << hillvar2<< std::endl; 
    double c = 0.05*dparam+(0.1*dparam)/(1+std::pow(EL/0.5,-1*nhill4));
double d= ( 1+std::pow((hillvar*WB/khill),nhill));
double e= ( 1+std::pow((hillvar2/khill2),-1*nhill2));
double f= ( 1+std::pow((hillvar3/khill3),nhill3));


  double contribution = -1*(0.1*dparam+(c / (d*e*f)))*y[compartment.index()][varIndex];

 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
  }
}


///////////////////////////////////////////////



Hill_Gated::Hill_Gated(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=4 ) {
		std::cerr << "Hill_Gated::Hill_Gated() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 4 ) {
		std::cerr << "Hill_Gated::Hill_Gated() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("hill_Gated");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
  tmp[0] = "wo1";
  tmp[1] = "wo2";
  tmp[2] = "nparam";
  tmp[3] = "klow";

	setParameterId( tmp );
}

void Hill_Gated::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
    double offsetval=0;
  double wo1=parameter(0);
  double wo2=parameter(1);
   double nparam=parameter(2);
    double klow=parameter(3);

  double hillvar=y[compartment.index()][variableIndex(0,0)];
      double xvar=y[compartment.index()][variableIndex(0,1)];
    double yvar=y[compartment.index()][variableIndex(0,2)];
    double zvar=y[compartment.index()][variableIndex(0,3)];
  double contributionprelim=0;
  double contribution=0;
  double Co=10;
    
    //play around with CLV3 gradient
    double b= -1*((std::pow(xvar,2)+std::pow(yvar,2)+std::pow(zvar-10,2))/3);
    Co=10*std::pow(2.71828,b); //stepper
    //Co=10*2.71828*b;



     if (hillvar>wo1)//what about equals?
    contributionprelim = 2*Co/ (1+std::pow(((hillvar-wo1)/(wo2-wo1)),nparam));
    
    else if(hillvar<wo1)
        contributionprelim = 2*Co/ (1+std::pow(((wo1-hillvar)/klow),nparam));
        
        contribution=hillvar*contributionprelim;
     

//Check to make control statements don't flow into each other



 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
  }
  
  
  
}

//////////////////////////////////////////////////////////////////


MutatableNuclearExport::MutatableNuclearExport(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=6 ) {
		std::cerr << "MutatableNuclearExport::MutatableNuclearExport() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 1 ) {
		std::cerr << "MutatableNuclearExport::MutatableNuclearExport() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("mutatableNuclearExport");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	  tmp[0] = "rnc";
  tmp[1] = "ro";
  tmp[2] = "Ro2";
  tmp[3] = "WB";
  tmp[4] = "n5";
  tmp[5] = "EL";

	setParameterId( tmp );
}

void MutatableNuclearExport::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
    double offsetval=0;
   double rnc = parameter(0);
  double ro = parameter(1);
  double Ro2=parameter(2);
    double WB=parameter(3);
    double n5=parameter(4);
    double EL=parameter(5);
    double varmne=y[compartment.index()][variableIndex(0,0)];
 
  double contribution=0;
    
  
  double c = (Ro2-ro)/(1+(std::pow(2*WB,n5)));
   
   
     contribution = varmne*rnc*(2*(ro + c )/(1 + std::pow(EL,-1*n5)));

 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
  }


}

/////////////////////////////////////////////////////////////////////////


Multiplicative::Multiplicative(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=1 ) {
		std::cerr << "Multiplicative::Multiplicative() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 5 ) {
		std::cerr << "Multiplicative::Multiplicative() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("multiplicative");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	  tmp[0] = "k";

	setParameterId( tmp );
}

void Multiplicative::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
    double offsetval=0;
double kon=parameter(0);
double firstvar=y[compartment.index()][ variableIndex(0,0) ];
double secondvar=y[compartment.index()][ variableIndex(0,1) ];
double xvar=y[compartment.index()][variableIndex(0,2)];
double yvar=y[compartment.index()][variableIndex(0,3)];
double zvar=y[compartment.index()][variableIndex(0,4)];


   double contribution = kon*firstvar*secondvar;

 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
  }

}


/////////////////////////////////////////////////////////////////////////

CreationTest::CreationTest(std::vector<double> &paraValue, 
												 std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indices
	if( paraValue.size()!=1 ) {
		std::cerr << "CreationTest::CreationTest() "
			  << "Uses only one parameter k_c\n";
		exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 4  ) {
		std::cerr << "CreationTest::CreationTest() "
							<< "One variable index needed.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("creationTest");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	tmp[0] = "k_c";
	setParameterId( tmp );
}

void CreationTest::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{



	dydt[compartment.index()][varIndex] += parameter(0)*
		y[compartment.index()][ variableIndex(0,0) ];
		
		//std::cerr << "parameter" << parameter(0) << std::endl;
		//std::cerr << "variable" << y[compartment.index()][ variableIndex(0,0) ] << std::endl;

	//y[compartment.index()][varIndex]=100; 
	
/*	
	
		if ((varIndex)==6)
{

y[compartment.index()][varIndex-1]=100; 


std::cerr << "derivative" << dydt[compartment.index()][varIndex-1] << std::endl;


	std::cerr << "parameter" << parameter(0) << std::endl;
	std::cerr << "WusMRNA" << y[compartment.index()][ variableIndex(0,0) ] << std::endl;
	std::cerr << "contribution" << parameter(0)*y[compartment.index()][ variableIndex(0,0) ] << std::endl;
	
		}
*/		
		//subtracted 1 from varIndex
		
}

/////////////////////////////////////////////////////////////////////////

DegradationTest::DegradationTest(std::vector<double> &paraValue, 
												 std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indices
	if( paraValue.size()!=1 ) {
		std::cerr << "DegradationTest::DegradationTest() "
			  << "Uses only one parameter k_c\n";
		exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 1  ) {
		std::cerr << "DegradationTest::DegradationTest() "
							<< "One variable index needed.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("degradationTest");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	tmp[0] = "k_c";
	setParameterId( tmp );
}

void DegradationTest::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{



	dydt[compartment.index()][varIndex] += parameter(0)*
		y[compartment.index()][ variableIndex(0,0) ];
		
		//std::cerr << "parameter" << parameter(0) << std::endl;
		//std::cerr << "variable" << y[compartment.index()][ variableIndex(0,0) ] << std::endl;

	//y[compartment.index()][varIndex]=100; 
	
		
}

/////////////////////////////////////////////////////////////////////////














CreationLimited::CreationLimited(std::vector<double> &paraValue, 
												 std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indices
	if( paraValue.size()!=1 ) {
		std::cerr << "CreationLimited::CreationLimited() "
			  << "Uses only one parameter k_c\n";
		exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 3  ) {
		std::cerr << "CreationLimited::CreationLimited() "
							<< "One variable index needed.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("creationLimited");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	tmp[0] = "k_c";
	setParameterId( tmp );
}

void CreationLimited::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{
    double offsetval=0;
	double xvar=y[compartment.index()][variableIndex(0,0)];
  double yvar=y[compartment.index()][variableIndex(0,1)];
  double zvar=y[compartment.index()][variableIndex(0,2)];



 	double contribution= parameter(0);
 

if (((xvar*xvar)+(yvar*yvar)+std::pow(zvar,2)) <=36)
 {
 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
 // std::cerr  << "contribution " << contribution << std::endl;
  //std::cerr  << "wRNA " << dydt[compartment.index()][varIndex]  << std::endl;
  }
  }
	
	
	
		
}



///////////////////////////////////////////////////////////////////////

DegradationOne_alt::DegradationOne_alt(std::vector<double> &paraValue, 
			       std::vector< std::vector<size_t> > 
			       &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=1 ) {
		std::cerr << "DegradationOne_alt::DegradationOne_alt() "
			  << "Uses only one parameter k_deg (lambda)\n";
		exit(0);
	}
	if( indValue.size() ) {
		std::cerr << "DegradationOne_alt::DegradationOne_alt() "
			  << "No variable index are used.\n";
		exit(0);
	}

	// Set the variable values
	setId("degradationOne_alt");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	tmp[0] = "k_deg";
	setParameterId( tmp );
}

void DegradationOne_alt::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
	dydt[compartment.index()][varIndex-1] -= parameter(0)*
		y[compartment.index()][varIndex-1];
		//subtracted 1 from varIndex
}




////////////////////////////////////////////////////////////

DiffusionSimple_alt::DiffusionSimple_alt(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > &indValue ) 
{  
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "DiffusionSimple_alt::DiffusionSimple_alt() "
	      << "Uses only one parameter D.\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "DiffusionSimple_alt::DiffusionSimple_alt() "
	      << "No variable index used.\n";
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("diffusionSimple_alt");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "D";
  setParameterId( tmp );
}

void DiffusionSimple_alt::
derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt) 
{    
   //if (compartment.index()==300)
 //std::cerr << " cell " << compartment.index() << " neighbors " << compartment.numNeighbor() << " Conc " << y[compartment.index()][species] << " c.\n";

  size_t i=compartment.index();
  for( size_t n=0 ; n<compartment.numNeighbor() ; n++ ) {
    size_t j=compartment.neighbor(n);
    //Only update if compartments index less than neighbors
      //std::cerr << " i "
	     // << i << " neighbors " << compartment.numNeighbor() << " j " << j << " Diffuse.\n";
    
    //110818 hack to prevent diffusion from turning concentrations negative. Should review later
    if( i<j ) {
    
    
      double diff = parameter(0)*(y[j][species]-y[i][species]);
      
         if (j==300)
     std::cerr << " Wus Cyto Concentration before diffusion " << y[j][species] << " Conc Change "<< diff<< " c.\n";
      
      
      if (y[j][species]-diff < 0)
      {
        
      if (y[j][species] > 0)
      diff =y[j][species];
      
      else 
      diff=0;
      
      }
      
            if (y[i][species]+diff < 0)
      {
      
    if (y[i][species] > 0)
      diff =y[i][species];
      
      else 
      diff=0;
      
      }
      
      
      
      dydt[i][species] += diff;
      dydt[j][species] -= diff;
   
   
     if (j==300)
   std::cerr << " 300 Conc Dy/dt " << dydt[j][species] << " c.\n";
   
    }
    
    
  }
}


///////////////////////////////////////////////////////////////////
CreationLinear_alt::CreationLinear_alt(std::vector<double> &paraValue, 
			       std::vector< std::vector<size_t> > 
			       &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=1 ) {
		std::cerr << "CreationLinear_alt::CreationLinear_alt() "
			  << "Uses only one parameter k_c\n";
		exit(0);
	}
	if( indValue.size() != 0 && indValue.size() != 1 ) {
		std::cerr << "CreationLinear_alt::CreationLinear_alt() "
			  << "At most one level of variables are allowed.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("creationLinear");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	tmp[0] = "k_c";
	setParameterId( tmp );
}

void CreationLinear_alt::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
	double rate = parameter(0);
	if( numVariableIndexLevel() )
	  for( size_t k=0; k<numVariableIndex(0); ++k )
	    rate *= y[compartment.index()][ variableIndex(0,k) ];
	dydt[compartment.index()][varIndex-1] += rate;
	//subtracted 1 from varIndex
}

////////////////////////////////////////////////


///////////////////////////////////////////////



CLV3_Dynamics::CLV3_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue )
{
	// Do some checks on the parameters and variable indeces
    if( paraValue.size()!=43 ) {
		std::cerr << "CLV3_Dynamics::CLV3_Dynamics() "
                  << "Uses only one parameter k_deg\n"
                  << "parameter(0)" << " 1st parameter "
                  << "parameter(1)" << " 2nd parameter "
                  << paraValue.size() << " parameter size\n";
        exit(0);
	}
    if( indValue.size() != 1 || indValue[0].size() != 5 ) {
		std::cerr << "CLV3_Dynamics::CLV3_Dynamics() "
                  << "Two variable indices are used.\n";
		exit(0);
	}

	// Set the variable values
	setId("cLV3_Dynamics");
    setParameter(paraValue);
	setVariableIndex(indValue);

	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
    tmp[0] = "clv3P";
    tmp[1] = "wusMonomerCoefficient";
    tmp[2] = "wusDimerCoefficient";
    tmp[3] = "clv3Flux";
    tmp[4] = "clv3DecayParameter";
    tmp[5] = "crmActivityCoefficient";
    tmp[6] = "clv3SourceWidth";
    tmp[7] = "timeStep";
    tmp[8] = "satpoint";
    tmp[9] = "cooptMonEffect";
    tmp[10] = "cooptDimEffect";

    tmp[11] = "bind1";
    tmp[12] = "bind2";
    tmp[13] = "bind3";
    tmp[14] = "bind4";
    tmp[15] = "bind5";
    tmp[16] = "unbind1";
    tmp[17] = "unbind2";
    tmp[18] = "unbind3";
    tmp[19] = "unbind4";
    tmp[20] = "unbind5";
    tmp[21] = "activationOnlyFlag";
    tmp[22] = "clv3Barrier";
    tmp[23] = "activation";
    tmp[24] = "suppressDecay";
    tmp[25] = "FrozenWus";
    tmp[26] = "CRMorMarkerSwitch";
    tmp[27] = "crmTimerLength";
    tmp[28] = "m4Flag";
    tmp[29] = "970BonusCoopt";
    tmp[30] = "neighborOnlyCoopt";
    tmp[31] = "dimerBindP";
    tmp[32] = "polBaseBindAffinity";
    tmp[33] = "polTimeLimit";
    tmp[34] = "monFireLimit";
    tmp[35] = "dimerUnbindP1";
    tmp[36] = "dimerUnbindP2";
    tmp[37] = "dimerUnbindP3";
    tmp[38] = "dimerUnbindP4";
    tmp[39] = "dimerUnbindP5";
    tmp[40] = "L1nodimer";
    tmp[41] = "bonusL1MonCoopt";
    tmp[42] = "unbindLimit";
	setParameterId( tmp );
}

void CLV3_Dynamics::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt)
{
//parameters
    double clv3P=parameter(0);
    double wusMonomerCoefficient=parameter(1);
    double wusDimerCoefficient=parameter(2);
    double clv3Flux=parameter(3);
    double clv3DecayParameter=parameter(4);
    double crmActivityCoefficient=parameter(5);
    double clv3SourceWidth=parameter(6);
    double timeStep=parameter(7);
    double satpoint=parameter(8);
    double cooptMonEffect= parameter(9);
    double cooptDimEffect= parameter(10);
    double activationOnlyFlag= parameter(21);
    double clv3Barrier=parameter(22);
    double activation=parameter(23);
    double suppressDecay=parameter(24);
    double frozenWUS=parameter(25);
    int crmOrMarkerSwitch=parameter(26);
    int crmTimerLength=parameter(27);
    int m4Flag=parameter(28);
    int HABonusCoopt=parameter(29);
    int neighborOnlyCoopt=parameter(30);
    double dimerBindP=parameter(31);
    double polBaseBindAffinity=parameter(32);
    double polTimeLimit=parameter(33);
    int monFireLimit=parameter(34);
    double dimerUnbindP [5] ={parameter(35), parameter(36), parameter(37), parameter(38), parameter(39)};//modify later
    double L1nodimer =parameter(40);
    int bonusL1MonCoopt = parameter(41); //added 122020
    int unbindLimit = parameter(42);

//variables
    double hillvar=y[compartment.index()][variableIndex(0,0)];
    double xvar=y[compartment.index()][variableIndex(0,1)];
    double yvar=y[compartment.index()][variableIndex(0,2)];
    double zvar=y[compartment.index()][variableIndex(0,3)];



//CRM variables
    int crmBindingSites=4;
    double geneCRMSiteBindMaxBaseChance [5] ={parameter(11), parameter(12), parameter(13), parameter(14), parameter(15)};//modify later
    double geneCRMSiteChance_Unbind [5] ={parameter(16), parameter(17), parameter(18), parameter(19), parameter(20)};//modify later. If site binding goes to zero what of unbinding?
    int wusMonomer=0;
    int wusDimer=0;
    int wusMonomer2=0;
    int wusDimer2=0;
    int chromoCycle;
    int emptyCRM=0;
    int stochTimeSkip=0;

    double clv3Creation=0;
    double clv3Decay;
    double timeStepIncrement, previousTimeStepIncrement, timeStepRemainder,timeStepOverflow, timeStepInteger;
    timeStepIncrement=0;

//loop iterators
    int i, j, l, n , p, q;

    int eventFlag=0;
    //int neighbors=0;
    //double cooperativeEffect=0;
    double contribution=0;
    double concModifier=0;
    double timeStepRemain=timeStep;
    std::ofstream outfile;//file testing
    //std::ofstream occupancyfile;//file testing

    //cylinder definition for expression domain
    double centralAxis[3]={0,0,zvar};
    double distanceFromCentralAxis=sqrt(pow(centralAxis[0]-xvar,2)+pow(centralAxis[1]-yvar,2)+pow(centralAxis[2]-zvar,2));

    //inner layer clv3 effectiveness
    double distanceFromCentralBase=sqrt(pow(xvar,2)+pow(yvar,2)+pow(zvar,2));
//rand variables
    double randvalue1, randvalue2;
    double eventChance, eventChance2;

    //seed pseudorandom num.
    //select from Mer Twist

    std::random_device rd;
    std::mt19937 mt(rd());

//do you need to adjust distribution to current timestep size?
    std::uniform_real_distribution<double> dist(0,1);


    //Check to make control statements don't flow into each other



    int stochasticLoopTrack=0;

    //1. OVERFLOW HANDLING: Handle overflow from previous timestemp

    OverFlowResults overFlowResults1;

    timeStepOverflow=compartment.stochasticStepOverflow[crmOrMarkerSwitch];

    //1b: Advance crm timers, generate activation, and remaining timestep based on overflow
    TimeStepOverflowHandler(timeStep, polTimeLimit, timeStepOverflow, crmOrMarkerSwitch,
                            clv3P, timeStepRemain, compartment, overFlowResults1);

    timeStepRemain= overFlowResults1.timeStepRemain;//gets remaining timestep from timestepoverflowhandler
    clv3Creation += overFlowResults1.clv3;//if 0 overflow than full timestep is given back.

    //2: PreStochastic Loop Safety Checks:
// check crm occupancy if it is empty
    //removed marker crm check as its no longer used 010520
    if (std::all_of(compartment.geneCRM[0],compartment.geneCRM[0]+5, [](int x){
        return x == 0;
    })
        && std::all_of(compartment.geneCRM[1],compartment.geneCRM[1]+5, [](int x){
        return x == 0;
    }))
    {emptyCRM=1;}

//2b: no going through stochastic loop if low CLV3 concentration and empty CRM
    //this prevents a weakness in the Gillespie algorithm of timesteps spiking to high values due to not much going on in the beginning
    if (hillvar < 1 && emptyCRM==1)
    {
        timeStepRemain = 0;
        compartment.previousTimeStep[crmOrMarkerSwitch]=compartment.latestTimeStep[crmOrMarkerSwitch];
        compartment.latestTimeStep[crmOrMarkerSwitch]=1;


    }


    //3: STOCHASTIC TIMESTEP LOOP: Take stochastic timesteps within ODE timestep
    //Consists of five main steps: CRMEventPicker, CRMProbabilityGenerator, TimeStepGenerator, CRMSummer, Clavata3ActivationMechanisms
    while (timeStepRemain>0)
    {

//4. Loop Preparations:
        //concentration affects CRM up to saturation point
        //scale to prevent fractional multiplication and large step sizes
        concModifier=hillvar/satpoint;

        if (concModifier>1)
            concModifier=1;



        //4b: Cleanup: done with eventNum, can zero
        compartment.eventNum[crmOrMarkerSwitch]=0;
//integrate in matrix erasure for safety later
//ResetProbabilityMatrix(compartment);
        //4c: zero deltasum (probability matrix magnitude) in preparation for generating a new one
        compartment.probabilityDeltaSum[crmOrMarkerSwitch]=0;



        //5: CRMProbability Generator: generate an event set and associated probabilities

        //ultra low probabilities will automatically be zero due to limits of double
        //also there is some small fluctuation of the precise values shown in the debugger for some reason.


        if (crmOrMarkerSwitch==0)
        {
            for (n=0; n <= 1; n++)
            {
                for (i=0; i <= crmBindingSites; i++)
                {

                    chromoCycle=n;
                    //int samplearray[2][5] = {{0,1,1,0,0}, {0,0,0,11,13}};

                    CRMProbabilityGenerator(m4Flag, crmOrMarkerSwitch, chromoCycle, compartment.probabilityMatrix, compartment.geneCRM[n],
                                            i, crmActivityCoefficient, cooptMonEffect, cooptDimEffect,
                                            geneCRMSiteBindMaxBaseChance[i], geneCRMSiteChance_Unbind[i], concModifier, compartment.eventNum[crmOrMarkerSwitch],
                                            compartment.probabilityDeltaSum[crmOrMarkerSwitch],HABonusCoopt,neighborOnlyCoopt, dimerBindP, polBaseBindAffinity,
                                            compartment, dimerUnbindP[i], L1nodimer, bonusL1MonCoopt, distanceFromCentralBase);

                }
            }
        }

        if (crmOrMarkerSwitch==1)
        {
            for (n=0; n <= 1; n++)
            {
                for (i=0; i <= crmBindingSites; i++)//fixed bug missing = sign
                {

                    chromoCycle=n;
                    //should the deltasum for the two chromosomes be seperate or together? Currently they are together
                    CRMProbabilityGenerator(m4Flag, crmOrMarkerSwitch, chromoCycle, compartment.probabilityMatrix, compartment.geneCRMMarker[n],
                                            i, crmActivityCoefficient, cooptMonEffect, cooptDimEffect,
                                            geneCRMSiteBindMaxBaseChance[i], geneCRMSiteChance_Unbind[i], concModifier, compartment.eventNum[crmOrMarkerSwitch],
                                            compartment.probabilityDeltaSum[crmOrMarkerSwitch],HABonusCoopt,neighborOnlyCoopt, dimerBindP, polBaseBindAffinity,
                                            compartment, dimerUnbindP[i], L1nodimer, bonusL1MonCoopt, distanceFromCentralBase);


                }
            }
        }



        //6: pick a random value to determine the event that takes place from the proability matrix
        std::uniform_real_distribution<double> dist2(0, compartment.probabilityDeltaSum[crmOrMarkerSwitch]);
        randvalue2 = dist2(mt);//need a new seed?





        //if this is not included than a monomer binding will automatically occur on the
        //first step for cells which have a certain threshold of WUSNuc.

        //7: Step 1: Event Picker: Select Event to happen from probability matrix
        //Event Picker is first because...?
        //After first timestep, update CRM, determine which event has happened. Go through each event to see if random value falls within it. Check if k goes all the way
        if(compartment.timer > 1)
        {//only take into account events calculated for this timestep. -1 eventnum
            for(l=0; l< compartment.eventNum[crmOrMarkerSwitch]; l++)
            {

                chromoCycle=compartment.probabilityMatrix[l].chromosome;



                if (crmOrMarkerSwitch==0)
                {



                    eventFlag=CRMEventPicker(compartment.probabilityMatrix[l].site,eventFlag, compartment.probabilityMatrix[l].action,
                                             compartment.probabilityMatrix[l].begin, compartment.probabilityMatrix[l].end,y[compartment.index()][variableIndex(0,0)],randvalue2,
                                             compartment.eventNum[crmOrMarkerSwitch], crmTimerLength, compartment.crmTimer[crmOrMarkerSwitch][compartment.probabilityMatrix[l].site],
                                             compartment.geneCRM[chromoCycle][compartment.probabilityMatrix[l].site], compartment, chromoCycle, monFireLimit, crmOrMarkerSwitch, polTimeLimit, unbindLimit);
                }


                else if(crmOrMarkerSwitch==1)
                {

                    eventFlag=CRMEventPicker(compartment.probabilityMatrix[l].site,eventFlag, compartment.probabilityMatrix[l].action,
                                             compartment.probabilityMatrix[l].begin, compartment.probabilityMatrix[l].end,y[compartment.index()][variableIndex(0,0)],randvalue2,
                                             compartment.eventNum[crmOrMarkerSwitch], crmTimerLength, compartment.crmTimer2[crmOrMarkerSwitch][compartment.probabilityMatrix[l].site],
                                             compartment.geneCRMMarker[chromoCycle][compartment.probabilityMatrix[l].site], compartment, chromoCycle, monFireLimit, crmOrMarkerSwitch, polTimeLimit, unbindLimit);
                }


                //5b: if an event has be chosen, break out
                if (eventFlag==1)
                {
                    eventFlag=0;
                    break;
                }


            }





        }



        //8: TimeStepGenerator: Generate new timestep
        //Generate a random value to use for timestep
        randvalue1=dist(mt);

        //7b: the very first step, probabilityDeltaSum will be zero. How to handle?
        //in cases where total probabilities equal zero, or a negative number. This will break out of the loop for that.
        //Check if right later
        if (compartment.probabilityDeltaSum[crmOrMarkerSwitch] <= 0)
        {
            break;
        }
        //otherwise generate timestep
        else if (compartment.probabilityDeltaSum[crmOrMarkerSwitch] > 0)
        {
            timeStepIncrement = TimeStepGenerator(crmOrMarkerSwitch, randvalue1, polTimeLimit, compartment);


        }

        //If the stepsize is larger than the ODE timestep it should be changed to 1. A stochastic timestep that skips over
        //several ODE timesteps does not make sense since there should be feedback from the ODE components. This is also superior to simply passing
        //to the next timestep since it won't be inadvertantly shorter than any other timestep. Check with Weitao.

        //7c: Subtract stochastic timestep from remaining timestep
        timeStepRemain = timeStepRemain-timeStepIncrement;



        //9: CRMSummer evaluate CRM for dimers and monomers bound
        //modulate activator and repressor contribution by the amount of monomer and dimer present
        for (j=0; j<crmBindingSites; j++)//Correct later
        {
            if (crmOrMarkerSwitch==0)
            {
                CRMSummer(compartment.geneCRM[0][j], j, compartment.crmTimer[0][j], wusMonomer, wusDimer, timeStepIncrement);
                CRMSummer(compartment.geneCRM[1][j], j, compartment.crmTimer[1][j], wusMonomer, wusDimer, timeStepIncrement);
            }

            else if (crmOrMarkerSwitch==1)
            {
                CRMSummer(compartment.geneCRMMarker[0][j], j, compartment.crmTimer2[0][j], wusMonomer2, wusDimer2, timeStepIncrement);
                CRMSummer(compartment.geneCRMMarker[1][j], j, compartment.crmTimer2[1][j], wusMonomer2, wusDimer2, timeStepIncrement);
            }
        }




        //10: Clavata3ActivationMechanisms: Calculate Activation for current CRM Configuration

        //if within proper range clv3 activation takes place
        if (distanceFromCentralAxis<=clv3SourceWidth && zvar>=clv3Barrier)
        {


            if (crmOrMarkerSwitch==0)
            {
                chromoCycle=0;//unused
                clv3Creation +=Clavata3ActivationMechanisms(activation, clv3Creation, clv3P, wusMonomer, wusDimer, wusMonomerCoefficient,
                                                            wusDimerCoefficient, chromoCycle, compartment, polTimeLimit, crmOrMarkerSwitch);
            }

            else if (crmOrMarkerSwitch==1)
            {
                chromoCycle=0;//unused
                clv3Creation += Clavata3ActivationMechanisms(activation, clv3Creation, clv3P, wusMonomer2, wusDimer2, wusMonomerCoefficient,
                                                             wusDimerCoefficient, chromoCycle, compartment, polTimeLimit, crmOrMarkerSwitch);
            }

        }

        else //no activation
            clv3Creation=0;



//11: Per Cycle CleanUp
        //zero out monomers and dimers after use
        wusMonomer=0;
        wusDimer=0;
        wusMonomer2=0;
        wusDimer2=0;

        contribution += (clv3Creation);
        clv3Creation=0;
        //scale activation by time step increment, implement mean later?

        stochasticLoopTrack++;


    }


//12: Outside Loop Cleanup
    //finished tracking loop
    stochasticLoopTrack=0;

    //Track times in clv3 dynamics: Do not let binding within first timestep?
    compartment.timer=compartment.timer+1;

//after binding loop track if stochastic timestep overflows to next dynamics timestep
    if (timeStepRemain<0)
        compartment.stochasticStepOverflow[crmOrMarkerSwitch]=-1*timeStepRemain;


//if theres exactly zero time remaining, no overflow (possibly redundant)



    //13: Final Output Calculations
//Depending on sim conditions, have a decay rate for CLV3 or not
    if (suppressDecay==0)
        clv3Decay=y[compartment.index()][varIndex]*clv3DecayParameter;

    else
        clv3Decay=0;


    compartment.clv3StepContribution=contribution;
    //With activation flag variable is set to output only to immediate clv3 generated value generated this cycle
    if (activationOnlyFlag==1)
    {
        //y[compartment.index()][varIndex] = (contribution/timeStep); //scale by timeStep for some reason Weitao wants this
        y[compartment.index()][varIndex] = (contribution);
    }

    else if (activationOnlyFlag==0)//normal clv3 output
    {
        //contribution=(contribution*timeStep) -clv3Decay;//times contribution by timestep?
        contribution=(contribution) -clv3Decay;

        //do not let concentration go to negative, include timestep modulation?
        if (y[compartment.index()][varIndex]+contribution < 0 )
        {
            y[compartment.index()][varIndex]=0;

        }

        else
        {//add contribution to derivative

            dydt[compartment.index()][varIndex] += contribution;
            //y[compartment.index()][varIndex] += contribution;

        }
    }

}


///////////////////////////////////////////////////////////////////
CrmupdateM::CrmupdateM(std::vector<double> &paraValue,
                       std::vector< std::vector<size_t> >
                       &indValue )
{
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=1 ) {
		std::cerr << "CrmupdateM::CrmupdateM() "
                  << "Uses only zero parameter k_c\n";
		exit(0);
	}
	if( indValue.size() != 0 && indValue.size() != 1 ) {
		std::cerr << "CrmupdateM::CrmupdateM() "
                  << "At most one level of variables are allowed.\n";
		exit(0);
	}

	// Set the variable values
	setId("crmupdatem");
    setParameter(paraValue);
	setVariableIndex(indValue);

	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
    tmp[0] = "crmsetting";
	setParameterId( tmp );
}

void CrmupdateM::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt)
{
    double wusIsoformValue=y[compartment.index()][variableIndex(0,0)];
    double wusIsoform=0;
    int crmOrMarkerSwitch=parameter(0);
    //dydt[compartment.index()][varIndex] += -1*wusIsoformValue;


    //evaluate CRM for dimers and monomers bound


    for (int j=0; j<=4; j++)//Correct later
	{

        if (crmOrMarkerSwitch==0)
        {
            if (compartment.geneCRM[0][j]==1)
                wusIsoform +=1;
            if (compartment.geneCRM[1][j]==1)
                wusIsoform +=1;
        }

        if (crmOrMarkerSwitch==1)//will not work properly since both gene and marker go to same output?
        {
            if (compartment.geneCRMMarker[0][j]==1)
                wusIsoform +=1;
            if (compartment.geneCRMMarker[1][j]==1)
                wusIsoform +=1;
        }


	}


    // std::cerr << "variableIndex" << " " << varIndex << std::endl;

    y[compartment.index()][varIndex] = wusIsoform;




}


///////////////////////////////////////////////////////////////////
CrmupdateD::CrmupdateD(std::vector<double> &paraValue, 
			       std::vector< std::vector<size_t> > 
			       &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
	if( paraValue.size()!=1 ) {
		std::cerr << "CrmupdateD::CrmupdateD() "
			  << "Uses only one parameter k_c\n";
		exit(0);
	}
	if( indValue.size() != 0 && indValue.size() != 1 ) {
		std::cerr << "CrmupdateD::CrmupdateD() "
			  << "At most one level of variables are allowed.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("crmupdated");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
    tmp[0] = "crmsetting";
	setParameterId( tmp );
}

void CrmupdateD::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt)
{
    double wusIsoformValue=y[compartment.index()][variableIndex(0,0)];

    //dydt[compartment.index()][varIndex] += -1*wusIsoformValue;
    double wusIsoform=0;
    int crmOrMarkerSwitch=parameter(0);

    //evaluate CRM for dimers and monomers bound

    //evaluate CRM for dimers and monomers bound
    for (int j=0; j<=4; j++)       //Correct later
    {

        if (crmOrMarkerSwitch==0)
        {
            if (compartment.geneCRM[0][j]==2)
                wusIsoform +=1;
            if (compartment.geneCRM[1][j]==2)
                wusIsoform +=1;
        }

        if (crmOrMarkerSwitch==1)   //will not work properly since both gene and marker go to same output?
        {
            if (compartment.geneCRMMarker[0][j]==2)
                wusIsoform +=1;
            if (compartment.geneCRMMarker[1][j]==2)
                wusIsoform +=1;
        }


    }



    // if (compartment.index()==5)
    //std::cerr << "CRM state of cell 5" <<  compartment.geneCRM[0] << compartment.geneCRM[1] << compartment.geneCRM[2] << compartment.geneCRM[3] << compartment.geneCRM[4] << std::endl;



    y[compartment.index()][varIndex] = wusIsoform;

}


////////////////////////////////////////////////////

WUSNuc_Dynamics::WUSNuc_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
    if( paraValue.size()!=10 ) {
		std::cerr << "WUSNuc_Dynamics::WUSNuc_Dynamics() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}
    if( indValue.size() != 1 || indValue[0].size() != 7 ) {
		std::cerr << "WUSNuc_Dynamics::WUSNuc_Dynamics() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("wUSNuc_Dynamics");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	  tmp[0] = "wusRNAtransP";
  tmp[1] = "wusNucDegP";
  tmp[2] = "clv3_wusNucP";
  tmp[3] = "clv3_wusNucNexp";
  tmp[4] = "wusCytoInP";
  tmp[5] = "wusNucOutP";
  tmp[6] = "freeze";
  tmp[7] = "wusNuc_CkP";
  tmp[8] = "wusNuc_OuterLayerExportP";
  tmp[9] = "wusNuc_StabP";
  
	setParameterId( tmp );
}

void WUSNuc_Dynamics::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
    double offsetval=0;
 
 
 double wusRNAtransP=parameter(0);
 double wusNucDegP=parameter(1);
 double clv3_wusNucP=parameter(2);
 double clv3_wusNucNexp=parameter(3);
 double wusCyto_NucP=parameter(4);
 double wusNuc_CytoP=parameter(5);
 double freeze=parameter(6);
 double wusNuc_CkP=parameter(7);
double wusNuc_OuterLayerExportP=parameter(8);
double wusNuc_StabP=parameter(9);
 

 double xvar=y[compartment.index()][variableIndex(0,0)];
 double yvar=y[compartment.index()][variableIndex(0,1)];
 double zvar=y[compartment.index()][variableIndex(0,2)];
  double wusRNAvar=y[compartment.index()][variableIndex(0,3)];
    double clv3var=y[compartment.index()][variableIndex(0,4)];
    double wusCytovar=y[compartment.index()][variableIndex(0,5)];
      double ckvar=y[compartment.index()][variableIndex(0,6)];
double wusNucVar=y[compartment.index()][varIndex];

      double contribution=0;
      double outerLayerExportMod=0;

      //inner layer clv3 effectiveness addition
      double distanceFromCentralBase=sqrt(pow(xvar,2)+pow(yvar,2)+pow(zvar,2));

      double wn_Production, wn_Degradation, wn_NtoC, wn_CtoN;



   //std::cerr <<"parameterlist"<<dparam<<" "<< khill<< " " << nhill<< " " << nhill<< " " <<khill2<< " " << nhill2 << " " <<
   //" "<< WB <<" "<< EL <<" "<< hillvar << " " << hillvar2<< std::endl;

//Production-Degradation+Import-Export
//double contribution =((wusNucVar*wusNuc_CytoP)/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp)))/((ckvar*wusNuc_CkP)+1)
//wusRNAvar*wusRNAtransP-(wusNucDegP*y[compartment.index()][varIndex])+compartment.wusP_FromCyto -wusNuc_CytoP/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp))*y[compartment.index()][varIndex];

      //wusRNAvar*wusRNAtransP-(wusNucDegP*y[compartment.index()][varIndex])+compartment.wusP_FromCyto -(wusNuc_CytoP/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp))*y[compartment.index()][varIndex])/((ckvar*wusNuc_CkP)+1);



      if (distanceFromCentralBase > 7.7)
      outerLayerExportMod=wusNuc_OuterLayerExportP;



    //Production-Degradation+Import-Export
    //contribution =
    //wusRNAvar*wusRNAtransP-(wusNucDegP*wusNucVar)+compartment.wusP_FromCyto -outerLayerExportMod*((wusNucVar*wusNuc_CytoP)/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp)))/((ckvar*wusNuc_CkP)+1);

     wn_Production=wusRNAvar*wusRNAtransP;
      wn_Degradation=(wusNucDegP*wusNucVar)/(1+pow((wusNucVar/wusNuc_StabP),2));
      wn_NtoC=((wusNucVar*wusNuc_CytoP)/(1+outerLayerExportMod*pow(clv3var/clv3_wusNucP,clv3_wusNucNexp)))/((ckvar*wusNuc_CkP)+1);
       wn_CtoN=compartment.wusP_FromCyto;

contribution= wn_Production - wn_NtoC + wn_CtoN -wn_Degradation;
    //contribution = wusRNAvar*wusRNAtransP-(wusNucDegP*wusNucVar)/(1+wusNucVar/wusNuc_StabP)+compartment.wusP_FromCyto - ((wusNucVar*wusNuc_CytoP)/(1+outerLayerExportMod*pow(clv3var/clv3_wusNucP,clv3_wusNucNexp)))/((ckvar*wusNuc_CkP)+1);

    //(ckvar/threshold+1)
//double wusNProd=wusRNAvar*wusRNAtransP;
//double wusNDeg=(wusNucDegP*y[compartment.index()][varIndex]);
//double wusNImport=compartment.wusP_FromCyto;
//double wusNExportOld=wusNuc_CytoP/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp))*y[compartment.index()][varIndex];
//double wusNExport=(wusNuc_CytoP/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp))*y[compartment.index()][varIndex]/((ckvar*wusNuc_CkP)+1));

                       //outerLayerExportMod*((wusNucVar*wusNuc_CytoP)/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp)))/((ckvar*wusNuc_CkP)+1)
    //keep track of exported WUS
//compartment.wusP_FromNuc=outerLayerExportMod*((wusNucVar*wusNuc_CytoP)/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp)))/((ckvar*wusNuc_CkP)+1);
//compartment.wusP_FromNuc=((wusNucVar*wusNuc_CytoP)/(1+outerLayerExportMod*pow(clv3var/clv3_wusNucP,clv3_wusNucNexp)))/((ckvar*wusNuc_CkP)+1);
//exported wus is now useless
compartment.wusP_FromNuc=0;

//old export to cyto term
//-wusCyto_NucP*y[compartment.index()][varIndex]
//old equation
//wusRNAvar*wusRNAtransP-(wusNucDegP*y[compartment.index()][varIndex])+wusCyto_NucP*wusCytovar -wusNuc_CytoP/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp))*y[compartment.index()][varIndex];
//export term from wusnuc
//-wusNuc_CytoP/(1+pow(clv3var/clv3_wusNucP,clv3_wusNucNexp))*y[compartment.index()][varIndex]
//import term from wuscyto
//wusNuc_CytoP/(1+pow(clv3var/clv3_wusCytoP,clv3_wusCytoNexp))*wusNucvar
//first constant temp was adjusted from 1 to 40
// std::cerr << "In Nuclear Code " << compartment.index() << std::endl;

if (freeze==0)
{
 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
  }
}


}



////////////////////////////////////////////////////

WUSRNA_Dynamics::WUSRNA_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
    if( paraValue.size()!=9 ) {
		std::cerr << "WUSRNA_Dynamics::WUSRNA_Dynamics() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}//messed with error messages...check later AD102017
	if( indValue.size() != 1 || indValue[0].size() != 4 ) {
		std::cerr << "WUSRNA_Dynamics::WUSRNA_Dynamics() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("wUSRNA_Dynamics");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	tmp[0] = "wusRnaSource";
	  tmp[1] = "wusRnaSourceModP";
  tmp[2] = "wusRnaDegP";
  tmp[3] = "clv3RnaP";
  tmp[4] = "clv3RnaNexp";
  tmp[5] = "wusRnaSourceWidth";
  tmp[6] = "underL1Thickness";
  tmp[7] = "WUSRNAbarrier";
  tmp[8] = "freeze";
	setParameterId( tmp );
	
}

void WUSRNA_Dynamics::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  


    double offsetval=0;
    double wusRnaSource=parameter(0);
  double wusRnaSourceModP=parameter(1);
   double wusRnaDegP=parameter(2);
     double clv3_wusRnaP=parameter(3);
  double clv3_wusRnaNexp=parameter(4);
  double wusRnaSourceWidth=parameter(5);
  double underL1Thickness=parameter(6);
  double WUSRNAbarrier=parameter(7);
  int freeze=parameter(8);

  double xvar=y[compartment.index()][variableIndex(0,0)];
  double yvar=y[compartment.index()][variableIndex(0,1)];
  double zvar=y[compartment.index()][variableIndex(0,2)];
double clv3var=y[compartment.index()][variableIndex(0,3)];
double distanceFromCentralAxis;
double distanceFromCentralBase;

    double centralAxis[3]={0,0,zvar};
 
  
  double contribution=0;
  double wusRNAvar=y[compartment.index()][varIndex];

distanceFromCentralAxis=sqrt(pow(centralAxis[0]-xvar,2)+pow(centralAxis[1]-yvar,2)+pow(centralAxis[2]-zvar,2));
distanceFromCentralBase=sqrt(pow(xvar,2)+pow(yvar,2)+pow(zvar,2));

    //stepper test
             //Create WUSRNA along a central cylinder defined by wusRNASourceWidth except in the L1 layer as defined by underL1Thickness
             if (distanceFromCentralAxis <=wusRnaSourceWidth && distanceFromCentralBase < underL1Thickness && distanceFromCentralBase > WUSRNAbarrier )
 contribution = (wusRnaSourceModP*wusRnaSource)/ ( 1+std::pow((clv3var/clv3_wusRnaP),clv3_wusRnaNexp));
 
 //WUSRNA degradation
 contribution -= wusRnaDegP*y[compartment.index()][varIndex];
      

 if (freeze==0)
 {
      //if the change will put the species amount below zero
 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
//rezero the amount
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 //add the contribution
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
  }
  }

}

////////////////////////////////////////////////////

WUSCyto_Dynamics::WUSCyto_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
    if( paraValue.size()!=8 ) {
		std::cerr << "WUSCyto_Dynamics::WUSCyto_Dynamics() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}
    if( indValue.size() != 1 || indValue[0].size() != 4 ) {
		std::cerr << "WUSCyto_Dynamics::WUSCyto_Dynamics() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
	setId("wUSCyto_Dynamics");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	  tmp[0] = "wusRNAtransP";
  tmp[1] = "wusNucDegP";
  tmp[2] = "wusCytoOutP";
  tmp[3] = "wusNucInP";
  tmp[4] = "wusCytoP";
  tmp[5] = "wusCytoNexp";
  tmp[6] = "freeze";
  tmp[7] = "wusCyto_CkP";
  
	setParameterId( tmp );
}

void WUSCyto_Dynamics::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
    double offsetval=0;
 
 
 double wusRNAtransP=parameter(0);
 double wusCytoDegP=parameter(1);
 double wusNuc_CytoP=parameter(2);
 double wusCyto_NucP=parameter(3);
 //Addition 092519
 double clv3_wusCytoP=parameter(4);
 double clv3_wusCytoNexp=parameter(5);




 double freeze=parameter(6);
 double wusCyto_CkP=parameter(7);
 
  double wusRNAvar=y[compartment.index()][variableIndex(0,0)];
    double wusNucvar=y[compartment.index()][variableIndex(0,1)];
    double clv3var=y[compartment.index()][variableIndex(0,2)];//addition 092519
    double ckvar=y[compartment.index()][variableIndex(0,3)];
    double wusC_Prod, wusC_Deg, wusC_NtoC, wusC_CtoN;
    // y[compartment.index()][varIndex] =8;


   //std::cerr <<"y "<<yvar<<" cyto conc "<< y[compartment.index()][varIndex] << std::endl;

//Creation-Degradation-Export+Import
    //updated import and export by tying it with compartment to better ensure consistency

    wusC_Prod=wusRNAvar*wusRNAtransP;
    //wusC_Deg=((wusCytoDegP*y[compartment.index()][varIndex])/(1+pow(clv3var/clv3_wusCytoP,clv3_wusCytoNexp)))/((ckvar*wusCyto_CkP)+1);

    wusC_Deg=(wusCytoDegP*y[compartment.index()][varIndex])/((ckvar*wusCyto_CkP)+1);
    wusC_NtoC=compartment.wusP_FromNuc;
    wusC_CtoN=wusCyto_NucP*y[compartment.index()][varIndex];

    double contribution = wusC_Prod-wusC_Deg-wusC_CtoN+wusC_NtoC;

    compartment.wusP_FromCyto=wusC_CtoN;
    //double contribution =  wusRNAvar*wusRNAtransP-(((wusCytoDegP*y[compartment.index()][varIndex])/(1+pow(clv3var/clv3_wusCytoP,clv3_wusCytoNexp)))/((ckvar*wusCyto_CkP)+1))-wusCyto_NucP*y[compartment.index()][varIndex] + compartment.wusP_FromNuc;
//keep track of WUS protein imported to nucleus
//compartment.wusP_FromCyto=wusCyto_NucP*y[compartment.index()][varIndex];

//old export term
  //wusNuc_CytoP/(1+pow(clv3var/clv3_wusCytoP,clv3_wusCytoNexp))*wusNucvar

        //wusRNAvar*wusRNAtransP-((wusCytoDegP*y[compartment.index()][varIndex])/(1+pow(clv3var/clv3_wusCytoP,clv3_wusCytoNexp)))-wusCyto_NucP*y[compartment.index()][varIndex] +wusNuc_CytoP/(1+pow(clv3var/clv3_wusCytoP,clv3_wusCytoNexp))*wusNucvar;

//Old equation: Before 092519
//double contribution =
//wusRNAvar*wusRNAtransP-(wusCytoDegP*y[compartment.index()][varIndex])-wusCyto_NucP*y[compartment.index()][varIndex] +wusNuc_CytoP*wusNucvar;

//if (compartment.index() ==300)
//std::cerr <<" wuscytoConc "<<y[compartment.index()][varIndex] <<" wus creation "<< wusRNAvar*wusRNAtransP << " wus degradation "<< (wusCytoDegP*y[compartment.index()][varIndex]) << 
//" wus to nuc "<< -wusCyto_NucP*y[compartment.index()][varIndex] << " wus from nuc " << wusNuc_CytoP*wusNucvar << std::endl;


if (freeze==0)
{
 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }
 
 else
 {
  dydt[compartment.index()][varIndex] += contribution; 
  }
}


}

///////////////////////////////////////////////////

CLV3PEPTIDE_Dynamics::CLV3PEPTIDE_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue )
{
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=5 ) {
        std::cerr << "WUSRNA_Dynamics::WUSRNA_Dynamics() "
                            << "Uses only one parameter k_deg\n"
          << "parameter(0)" << " 1st parameter "
          << "parameter(1)" << " 2nd parameter "
          << paraValue.size() << " parameter size\n";
          exit(0);
    }//messed with error messages...check later AD102017
    if( indValue.size() != 1 || indValue[0].size() != 4 ) {
        std::cerr << "CLV3Peptide_Dynamics::CLV3Peptide_Dynamics() "
                            << "Two variable indices are used.\n";
        exit(0);
    }

    // Set the variable values
    setId("CLV3PEPTIDE_Dynamics");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
      tmp[0] = "clv3PeptideProdP";
  tmp[1] = "clv3PeptideDegP";
  tmp[2] = "freeze";
  tmp[3] = "activationOnlyFlag";
  tmp[4] = "innerLayerEffectiveness";
    setParameterId( tmp );

}

void CLV3PEPTIDE_Dynamics::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt)
{


    double offsetval=0;

  double clv3PeptideProdP=parameter(0);
  double clv3PeptideDegP=parameter(1);
int freeze=parameter(2);
  double activationOnlyFlag=parameter(3);
  double innerLayerEffectiveness=parameter(4);

double clv3RnaSource=y[compartment.index()][variableIndex(0,0)];
double xvar=y[compartment.index()][variableIndex(0,1)];
double yvar=y[compartment.index()][variableIndex(0,2)];
double zvar=y[compartment.index()][variableIndex(0,3)];

  double contribution=0;

  //inner layer clv3 effectiveness addition
  double distanceFromCentralBase=sqrt(pow(xvar,2)+pow(yvar,2)+pow(zvar,2));


 contribution = (clv3RnaSource*clv3PeptideProdP);



//clv3 effectiveness addition
if (distanceFromCentralBase<8.5)
contribution=contribution*innerLayerEffectiveness;


 //CLV3Peptide degradation
 contribution -= clv3PeptideDegP*y[compartment.index()][varIndex];


 if (freeze==0)
 {

if (activationOnlyFlag==0)
{
      //if the change will put the species amount below zero
 if (y[compartment.index()][varIndex]+contribution < 0 && contribution < 0)
{
//rezero the amount
offsetval=y[compartment.index()][varIndex]+contribution;
dydt[compartment.index()][varIndex] +=contribution-offsetval;
 }

 //add the contribution
 else
 {
  dydt[compartment.index()][varIndex] += contribution;
  }
  }

else if (activationOnlyFlag==1)// if activationOnlyFlag is set, CLV3Peptide is equivalent to CLV3RNA
 y[compartment.index()][varIndex] = clv3RnaSource;
 }


}



////////////////////////////////////////////////////

CLV3_Tracker::CLV3_Tracker(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue ) 
{  
	// Do some checks on the parameters and variable indeces
    if( paraValue.size()!=1 ) {
		std::cerr << "WUSNuc_Dynamics::WUSNuc_Dynamics() "
							<< "Uses only one parameter k_deg\n"
	      << "parameter(0)" << " 1st parameter " 
	      << "parameter(1)" << " 2nd parameter "
	      << paraValue.size() << " parameter size\n";
	      exit(0);
	}
	if( indValue.size() != 1 || indValue[0].size() != 1 ) {
		std::cerr << "WUSNuc_Dynamics::WUSNuc_Dynamics() "
							<< "Two variable indices are used.\n";
		exit(0);
	}
	
	// Set the variable values
    setId("cLV3_Tracker");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	// Set the parameter identities
	std::vector<std::string> tmp( numParameter() );
	tmp.resize( numParameter() );
	  tmp[0] = "conversionFactor";


  
  
	setParameterId( tmp );
}

void CLV3_Tracker::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt) 
{  
 
 
 double conversionFactor=parameter(0);
 
  //double CLV3Activation=y[compartment.index()][variableIndex(0,0)];
    double CLV3Activation=compartment.clv3StepContribution;
    

   //std::cerr <<"parameterlist"<<dparam<<" "<< khill<< " " << nhill<< " " << nhill<< " " <<khill2<< " " << nhill2 << " " <<
   //" "<< WB <<" "<< EL <<" "<< hillvar << " " << hillvar2<< std::endl;

//Production-Degradation+Import-Export
double contribution = CLV3Activation;

//first constant temp was adjusted from 1 to 40
// std::cerr << "In Nuclear Code " << compartment.index() << std::endl;


//add the contribution

 dydt[compartment.index()][varIndex] += contribution;
 



}


////////////////////////////////////////////////////


WUSP_ExportRate::WUSP_ExportRate(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue )
{
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=1 ) {
        std::cerr << "WUSP_ExportRate::WUSP_ExportRate() "
                            << "Uses only one parameter k_deg\n"
          << "parameter(0)" << " 1st parameter "
          << "parameter(1)" << " 2nd parameter "
          << paraValue.size() << " parameter size\n";
          exit(0);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
        std::cerr << "WUSP_ExportRate::WUSP_ExportRate() "
                            << "Two variable indices are used.\n";
        exit(0);
    }

    // Set the variable values
    setId("wUSP_ExportRate");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
      tmp[0] = "conversionFactor";




    setParameterId( tmp );
}

void WUSP_ExportRate::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt)
{


 double conversionFactor=parameter(0);

  //double CLV3Activation=y[compartment.index()][variableIndex(0,0)];


   //std::cerr <<"parameterlist"<<dparam<<" "<< khill<< " " << nhill<< " " << nhill<< " " <<khill2<< " " << nhill2 << " " <<
   //" "<< WB <<" "<< EL <<" "<< hillvar << " " << hillvar2<< std::endl;

//Production-Degradation+Import-Export
double contribution = compartment.wusP_FromNuc;

//first constant temp was adjusted from 1 to 40
// std::cerr << "In Nuclear Code " << compartment.index() << std::endl;


//add the contribution

 y[compartment.index()][varIndex] = contribution;




}


////////////////////////////////////////////////////


CkLigand_Dynamics::CkLigand_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue )
{
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=5 ) {
        std::cerr << "CkLigand_Dynamics::CkLigand_Dynamics() "
                            << "Uses only one parameter k_deg\n"
          << "parameter(0)" << " 1st parameter "
          << "parameter(1)" << " 2nd parameter "
          << paraValue.size() << " parameter size\n";
          exit(0);
    }//messed with error messages...check later AD102017
    if( indValue.size() != 1 || indValue[0].size() != 3 ) {
        std::cerr << "CkLigand_Dynamics::CkLigand_Dynamics() "
                            << "Two variable indices are used.\n";
        exit(0);
    }

    // Set the variable values
    setId("ckLigand_Dynamics");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "ckL_SourceP";
      tmp[1] = "ckL_ProdP";
  tmp[2] = "ckL_DegP";
  tmp[3] = "ckL_ReceptorAssP";
  tmp[4] = "ckL_ComplexDissP";
    setParameterId( tmp );

}

void CkLigand_Dynamics::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt)
{


    double ckL_SourceP=parameter(0);
  double ckL_ProdP=parameter(1);
  double ckL_DegP=parameter(2);
     double ckL_GradientCenterHeight=parameter(3);
     double ckL_ProdDiameter=parameter(4);


  double xvar=y[compartment.index()][variableIndex(0,0)];
  double yvar=y[compartment.index()][variableIndex(0,1)];
  double zvar=y[compartment.index()][variableIndex(0,2)];
//double ckL_ckReceptorVar=y[compartment.index()][variableIndex(0,3)];
//double ckL_ckComplexVar=y[compartment.index()][variableIndex(0,3)];
double distanceFromGradientCenter;
double ckL_SpatialMod;

  double contribution=0;

distanceFromGradientCenter=sqrt(pow(xvar,2)+pow(yvar,2)+pow(zvar-ckL_GradientCenterHeight,2));
ckL_SpatialMod=(-1*(xvar*xvar)/(3*3))-((yvar*yvar)/(3*3))-((zvar*zvar)/(6*6));
//Production-Degradation-cytokininonreceptor+cytokininoffreceptor
    //stepper test
             //Create cytokinin ligand along a central cylinder defined by wusRNASourceWidth except in the L1 layer as defined by underL1Thickness

if (distanceFromGradientCenter < ckL_ProdDiameter)
contribution = (ckL_ProdP*ckL_SourceP);

    //contribution = (ckL_ProdP*ckL_SourceP)*pow(2.718,ckL_SpatialMod) ;


                 //contribution = (wusRnaSourceModP*wusRnaSource)/ ( 1+std::pow((clv3var/clv3_wusRnaP),clv3_wusRnaNexp));

 //ckLigand degradation and complex dynamics. Should be located after the calculation of newAssComplex and newDissComplex
 contribution += -1*(ckL_DegP*y[compartment.index()][varIndex]) -(compartment.newAssCk)+(compartment.newDissCk);




 //add the contribution

  dydt[compartment.index()][varIndex] += contribution;



}

////////////////////////////////////////////////////


CkReceptor_Dynamics::CkReceptor_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue )
{
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=10 ) {
        std::cerr << "CkReceptor_Dynamics::CkReceptor_Dynamics() "
                            << "Uses only one parameter k_deg\n"
          << "parameter(0)" << " 1st parameter "
          << "parameter(1)" << " 2nd parameter "
          << paraValue.size() << " parameter size\n";
          exit(0);
    }//messed with error messages...check later AD102017
    if( indValue.size() != 1 || indValue[0].size() != 5 ) {
        std::cerr << "CkReceptor_Dynamics::CkReceptor_Dynamics() "
                            << "Two variable indices are used.\n";
        exit(0);
    }

    // Set the variable values
    setId("ckReceptor_Dynamics");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "ckR_SourceP";
     tmp[1] = "ckR_SourceDiameter";
      tmp[2] = "ckR_GradientHeight";
  tmp[3] = "ckR_ProdP";
  tmp[4] = "ckR_DegP";
  tmp[5] = "ckR_ComplexAssP";
  tmp[6] = "ckR_ComplexDissP";
  tmp[7] = "ckR_GradXaxe";
  tmp[8] = "ckR_GradYaxe";
  tmp[9] = "ckR_GradZaxe";
    setParameterId( tmp );

}

void CkReceptor_Dynamics::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt)
{

    double ckR_SourceP=parameter(0);
    double ckR_SourceDiameter=parameter(1);
    double ckR_GradientHeight=parameter(2);
  double ckR_ProdP=parameter(3);
  double ckR_DegP=parameter(4);
   double ckR_ComplexAssP=parameter(5);
   double ckR_ComplexDissP=parameter(6);
   double ckR_GradXaxe=parameter(7);
   double ckR_GradYaxe=parameter(8);
   double ckR_GradZaxe=parameter(9);



  double xvar=y[compartment.index()][variableIndex(0,0)];
  double yvar=y[compartment.index()][variableIndex(0,1)];
  double zvar=y[compartment.index()][variableIndex(0,2)];
  double ckR_LigandVar=y[compartment.index()][variableIndex(0,3)];
double ckR_ckComplexVar=y[compartment.index()][variableIndex(0,4)];
 double ovalValueFromGradientCenter;
 double newDissCk=0;
 double newAssCk=0;

    double centralAxis[3]={0,0,zvar};


  double contribution=0;

//distanceFromGradientCenter=sqrt(pow(xvar,2)+pow(yvar,2)+pow(zvar-ckR_GradientHeight,2));
  ovalValueFromGradientCenter=((xvar*xvar)/(ckR_GradXaxe*ckR_GradXaxe))+((yvar*yvar)/(ckR_GradYaxe*ckR_GradYaxe))+(((zvar-ckR_GradientHeight)*(zvar-ckR_GradientHeight))/(ckR_GradZaxe*ckR_GradZaxe));



    //stepper test
             //Maintain ckReceptor in a spherical domain under the L2. Dynamics should be before the ck_Ligand
             if (ovalValueFromGradientCenter < 1 )
             {

                 newAssCk=ckR_ComplexAssP*y[compartment.index()][varIndex]*ckR_LigandVar;
                 newDissCk=ckR_ComplexDissP*ckR_ckComplexVar;


                 //if calculated Association is higher than available ligands or receptors
                 if (newAssCk > ckR_LigandVar || newAssCk > y[compartment.index()][varIndex])
                 {
                     if (ckR_LigandVar < y[compartment.index()][varIndex])
                         newAssCk=ckR_LigandVar;

                     else
                         newAssCk=y[compartment.index()][varIndex];
                 }

                 if (newDissCk > ckR_ckComplexVar)
                 newDissCk = ckR_ckComplexVar;



                 contribution = ckR_ProdP+newDissCk-newAssCk-ckR_DegP*y[compartment.index()][varIndex];

             }

  compartment.newAssCk=newAssCk;
  compartment.newDissCk=newDissCk;


 //Should receptors have degradation?
// contribution -= ckR_DegP*y[compartment.index()][varIndex];



  dydt[compartment.index()][varIndex] += contribution;


}

////////////////////////////////////////////////////




CkComplex_Dynamics::CkComplex_Dynamics(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue )
{
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=1 ) {
        std::cerr << "CkComplex_Dynamics::CkComplex_Dynamics() "
                            << "Uses only one parameter k_deg\n"
          << "parameter(0)" << " 1st parameter "
          << "parameter(1)" << " 2nd parameter "
          << paraValue.size() << " parameter size\n";
          exit(0);
    }//messed with error messages...check later AD102017
    if( indValue.size() != 1 || indValue[0].size() != 3 ) {
        std::cerr << "CkComplex_Dynamics::CkComplex_Dynamics() "
                            << "Two variable indices are used.\n";
        exit(0);
    }

    // Set the variable values
    setId("ckComplex_Dynamics");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "ckC_DegP";
    setParameterId( tmp );

}

void CkComplex_Dynamics::
derivs(Compartment &compartment,size_t varIndex,DataMatrix &y,DataMatrix &dydt)
{


  double ckC_DegP=parameter(0);
  // double ckC_ReceptorAssP=parameter(1);
    // double ckC_ComplexDissP=parameter(2);


  double xvar=y[compartment.index()][variableIndex(0,0)];
  double yvar=y[compartment.index()][variableIndex(0,1)];
  double zvar=y[compartment.index()][variableIndex(0,2)];
//double ckC_ckLigandVar=y[compartment.index()][variableIndex(0,3)];
//double ckC_ckReceptorVar=y[compartment.index()][variableIndex(0,3)];


    double centralAxis[3]={0,0,zvar};


  double contribution=0;



    //stepper test
            //Association-Disassociation-Degradation
                 contribution = (compartment.newAssCk)-(compartment.newDissCk)-ckC_DegP*y[compartment.index()][varIndex];

                 //should be last in the ck series
                 compartment.newDissCk=0;
                 compartment.newAssCk=0;
 // degradation
 //contribution -= ckC_DegP*y[compartment.index()][varIndex];



  dydt[compartment.index()][varIndex] += contribution;


}

////////////////////////////////////////////////////


DiffusionSimple::DiffusionSimple(std::vector<double> &paraValue,
                 std::vector< std::vector<size_t> > &indValue )
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "DiffusionSimple::DiffusionSimple() "
          << "Uses only one parameter D.\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "DiffusionSimple::DiffusionSimple() "
          << "No variable index used.\n";
    exit(0);
  }
  //
  // Set the variable values
  //
  setId("diffusionSimple");
  setParameter(paraValue);
  setVariableIndex(indValue);
  //
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "D";
  setParameterId( tmp );
}

void DiffusionSimple::
derivs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt)
{
  size_t i=compartment.index();
  for( size_t n=0 ; n<compartment.numNeighbor() ; n++ ) {
    size_t j=compartment.neighbor(n);
    //Only update if compartments index less than neighbors
    if( i<j ) {
      double diff = parameter(0)*(y[j][species]-y[i][species]);
      dydt[i][species] += diff;
      dydt[j][species] -= diff;
    }
  }
}

void DiffusionSimple::
derivsWithAbs(Compartment &compartment,size_t species,DataMatrix &y,DataMatrix &dydt,DataMatrix &sdydt)
{
  size_t i=compartment.index();
  for( size_t n=0 ; n<compartment.numNeighbor() ; n++ ) {
    size_t j=compartment.neighbor(n);
    //Only update if compartments index less than neighbors
    if( i<j ) {
      double diff = parameter(0)*(y[j][species]-y[i][species]);
      dydt[i][species] += diff;
      dydt[j][species] -= diff;
      sdydt[i][species] += std::fabs(diff);
      sdydt[j][species] += std::fabs(diff);
    }
  }
}

//DISABLED 042521
/*
size_t DiffusionSimple::
Jacobian(Compartment &compartment,size_t species,DataMatrix &y,JacobianMatrix &A)
{
  size_t ci=compartment.index();
  for (size_t n=0; n<compartment.numNeighbor(); ++n) {
    size_t cj=compartment.neighbor(n);
    //Only update if compartments index less than neighbors
    if( ci<cj ) {
      size_t i = ci*y[0].size() + species;
      size_t j = ci*y[0].size() + species;
      A(i,i) -= parameter(0);
      A(j,j) -= parameter(0);
      A(i,j) += parameter(0);
      A(j,i) += parameter(0);
    }
  }
  return 1;
}

double DiffusionSimple::
propensity(Compartment &compartment,size_t species,DataMatrix &y)
{
  return parameter(0)*y[compartment.index()][species]*compartment.numNeighbor();
}

void DiffusionSimple::
discreteUpdate(Compartment &compartment,size_t species,DataMatrix &y)
{
    size_t neigh = static_cast<size_t> (myRandom::ran3()*compartment.numNeighbor());
    y[compartment.index()][species] -= 1.0;
    y[compartment.neighbor(neigh)][species] += 1.0;
}

*/

void DiffusionSimple::printCambium( std::ostream &os, size_t varIndex ) const
{
  std::string varName=organism()->variableId(varIndex);
  os << "Arrow[Organism[" << id() << "],{"
     << varName << "[i]},{},{"
     << varName << "[j]},";
  os << "Parameters["
     << parameter(0);
  os << "], ParameterNames[" << parameterId(0) << "], VarIndices[],";
  os << "Solve{" << varName << "[i]\' =  - p_0 ("
     << varName << "[i] - "
     << varName << "[j])"
     <<", ";
  os << varName << "[j]\' =   p_0 ("
     << varName << "[i] - "
     << varName << "[j])"
     <<", "
     << "Using[{nbr,"
     << organism()->topology().id() << "::Cell, i,"
     << organism()->topology().id() << "::Cell, j}"
     << "}";
  os << "]" << std::endl;
}

