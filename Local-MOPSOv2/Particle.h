#ifndef PARTICLE
#define PARTICLE

#define SIZE_D      32

typedef struct PARTICLE{
    int *Neighbor;
    double *Velocity;

    double *CurrentFitness;
    double *CurrentPosition;

    double *ParticleBestFitness;
    double *ParticleBestPosition;

}Particle;


#endif // PARTICLE
