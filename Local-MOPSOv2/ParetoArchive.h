#ifndef PARETOARCHIVE_H
#define PARETOARCHIVE_H

#include "Particle.h"

typedef struct PARETOPARTICLE{
    double *Position;
    double *Fitness;
    double Density;
}ParetoParticle;

typedef struct ARCHIVE{
    int CurrentSize;
    int **VectorLeader;
    double TotalDensity;
    Particle *Particles;
    ParetoParticle *Leader;
}Archive;

int TournamentSelection(Archive *, int);
double CalculateDensity(Archive *, double *, int );

#endif
