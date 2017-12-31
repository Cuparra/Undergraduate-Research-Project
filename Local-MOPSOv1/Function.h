
#ifndef FUNCTION_H
#define FUNCTION_H

#include <stdio.h>
#include <math.h>
#include "Swarm.h"

# define CONSTANT   0.18257418583
# define LIMIT      0.00001
# define A1         0.87364856231
# define A2         2.74857244327

//  0<= x <= 10, 6 variables,  2 functions, 6 constraints
void OsyczkaFunction(Particle *);

//  0<= x <= 1, 30 variables, 4 functions
void ZDT1(Particle *);

//  0<= x <= 1, 30 variables, 4 functions
void ZDT2(Particle *P);

//  0<= x <= 1, 30 variables, 4 functions
void ZDT3(Particle *);
//  -5 <= x <= 5, 10 variables, 4 functions
void ZDT4(Particle *);

//  0 <= x <= 1, 10 variables, 4 functions
void ZDT6(Particle *);

// -4 <= x <= 4, 30 variables, 2 functions
void FonFleFunction(Particle *);

// -7 <= x <= 4, 2 variables, 2 functions, 3 constraints
void TestFunction4(Particle *);

// 0 <= x <= 1, 2 variables, 2 functions, 2 constraints
void CTP1Function(Particle *);

// 0 <= x <= 5, 2 variables, 2 functions, 2 constraints
void BihnKornFunction(Particle *);

// -M_PI <= x <= M_PI, 2 variable, 2 functions
void PoloniFunction(Particle *);


double ErrorRation(Swarm *, FILE*);

#endif
