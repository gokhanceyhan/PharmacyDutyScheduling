

#ifndef _BS_H_
#define _BS_H_

#include <math.h>
#include <time.h> 
#include <ilcplex/ilocplex.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib> 
#include <fstream>
#include <iomanip>
#define M 1000000000
#define Imax 61
#define Jmax 91
#define Tmax 16
#define Kmax 10
#define InsMax 1
#define RC_EPS 1.0e-6
using namespace std;

ILOSTLBEGIN

/* The following 3 declarations are for use of the random-number generator
   lcgrand and the associated functions lcgrandst and lcgrandgt for seed
   management.  This file (named lcgrand.h) should be included in any program
   using these functions by executing
       #include "lcgrand.h"
   before referencing the functions. */

float lcgrand(int stream);
void  lcgrandst(long zset, int stream);
long  lcgrandgt(int stream);

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumArray> IloNumArrayArray; 
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix>   NumVar3Matrix;


typedef IloArray<IloBoolVarArray> IloBoolVarArray2;
typedef IloArray<IloBoolArray> IloBoolArrayArray;


#endif //_BS_H_
