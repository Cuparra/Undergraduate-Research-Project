#ifndef POPULATION_H
#define POPULATION_H

#include "Particle.h"
#include "ParetoArchive.h"

#define SIZE_ARCHIVE 10

#define MINIMUM     0.4
#define MAXIMUM     0.5

typedef struct SWARM{
    Archive *ParetoFront;
    Particle *Particles;
    int SizeOfNeighbor;
    int SizeOfFitness;      /*It includes objectives and constraints*/
    int SizeOfObjective;    /*Includes objectives only*/
    int SizeOfArchive;
    int SizeOfDimention;
    int SizeOfSwarm;
    double accelaration1;
    double accelaration2;
    int MaxInteraction;
    double weight;
    double HighestNumber;
    double LowestNumber;
    double Result;
}Swarm;

void UpdateArchive(Swarm *, Archive *);
int RouletteWheel(Archive *);

#endif
