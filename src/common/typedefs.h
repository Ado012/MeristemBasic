#ifndef TYPEDEFS_H
#define TYPEDEFS_H
//
// Filename     : typedefs.h
// Description  : Contains all typedefs used within the organism code
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Modded : AD
// Created      : February 2010
// Revision     : $Id:$
//
//#include "la/SparseMatrix.h"
//#include "la/DenseVector.h"
#include <vector>

using namespace std;

typedef std::vector< std::vector<double> > DataMatrix;
//typedef SparseMatrix JacobianMatrix;
//typedef DenseVector ConstantVector;

#endif
