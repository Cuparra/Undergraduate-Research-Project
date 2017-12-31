
#ifndef FUNCTION_H
#define FUNCTION_H

#include <stdio.h>
#include <math.h>
#include "Swarm.h"

# define CONSTANT   0.18257418583
# define LIMIT      0.00001
# define A1         0.87364856231
# define A2         2.74857244327

#define SIZE_D      30

//  0<= x <= 1, 30 variables, 4 functions
void ZDT1(Particle *);

//  0<= x <= 1, 30 variables, 4 functions
void ZDT3(Particle *);
//  -5 <= x <= 5, 10 variables, 4 functions
void ZDT4(Particle *);

//  0 <= x <= 1, 10 variables, 4 functions
void ZDT6(Particle *);

double ErrorZDT3(Swarm *, FILE*);
double ErrorZDT6(Swarm *, FILE*);
#endif
