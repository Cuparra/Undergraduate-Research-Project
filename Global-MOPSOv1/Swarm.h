#ifndef POPULATION_H
#define POPULATION_H

#include <stdlib.h>

typedef struct PARETOPARTICLE{
    double *Position;
    double *Fitness;
    double Density;
}ParetoParticle;

typedef struct ARCHIVE{
    int CurrentSize;
    int **VectorLeader;
    int SizeOfArchive;
    double TotalDensity;
    ParetoParticle *Leader;
}Archive;

typedef struct PARTICLE{
    int *Neighbor;
    double *Velocity;
    double *CurrentFitness;
    double *ParticleBestFitness;
    double *CurrentPosition;
    double *ParticleBestPosition;
}Particle;

typedef struct SWARM{
    Archive *ParetoFront;
    Particle *Particles;
    Particle *BestParticle;
    int (*Compare)(double,double);
    int SizeOfFitness;
    int SizeOfDimention;
    int SizeOfSwarm;
    double accelaration1;
    double accelaration2;
    int MaxInteraction;
    double weight;
    double Range;
    double HighestNumber;
    double LowestNumber;
    double sigma;
}Swarm;

void UpdateArchive(Swarm *);
int RouletteWheel(Archive *);

#endif
