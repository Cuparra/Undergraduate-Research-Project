#include "Inicialization.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

static double *GeneratePosition(int SizeOfDimention, double HighestNumber, double LowestNumber){

    double *Position = (double*)calloc(SizeOfDimention, sizeof(double));

    while(SizeOfDimention--)
        Position[SizeOfDimention]   = ( (double) rand() / RAND_MAX)*(HighestNumber-LowestNumber) + LowestNumber;

    return Position;

}

static void ConstructMesh(Swarm *S){

    int i,j;
    Archive *A          = S->ParetoFront;
    Particle *Particles = S->Particles;
    int SizeOfSwarm     = S->SizeOfSwarm;

    int SizeOfNeighbor  = sqrt(SizeOfSwarm);

    S->SizeOfNeighbor = 4;

    for(i=0; i<SizeOfSwarm;i++)
        Particles[i].Neighbor = (int*)malloc(4*sizeof(int));

    for(i=0; i<SizeOfNeighbor; i++){

        for(j=0; j<SizeOfNeighbor; j++){

            int pos = i*SizeOfNeighbor + j;

            Particles[pos].Neighbor[0] = (i-1)*SizeOfNeighbor + j;
            Particles[pos].Neighbor[1] = (i+1)*SizeOfNeighbor + j;
            Particles[pos].Neighbor[2] = pos - 1;
            Particles[pos].Neighbor[3] = pos + 1;

            if(i - 1 < 0)
               Particles[pos].Neighbor[0] = (SizeOfNeighbor - 1)*SizeOfNeighbor + j;

            if(i + 1 > SizeOfNeighbor - 1)
               Particles[pos].Neighbor[1] = j;

            if(pos - 1 < i*SizeOfNeighbor)
                Particles[pos].Neighbor[2] = (i+1)*SizeOfNeighbor - 1;

            if(pos + 1 > (i+1)*SizeOfNeighbor - 1)
                Particles[pos].Neighbor[3] = i*SizeOfNeighbor;

        }
    }

    for(i=0;i<SizeOfSwarm;i++){

        A[i].Particles = (Particle*)malloc(4*sizeof(Particle));

        for(j=0;j<4;j++){
            int idx           = Particles[i].Neighbor[j];
            A[i].Particles[j] = Particles[idx];
        }
    }

}

Swarm *InicializeSwarm(int SizeOfSwarm, int SizeOfDimention, int SizeOfFitness, int SizeOfObjective, double accelaration1, double accelaration2, double HighestNumber, double LowestNumber,  int MaxInteraction){

    int i,j,k;
    Swarm *S                = (Swarm*)malloc(sizeof(Swarm));
    S->SizeOfSwarm          = SizeOfSwarm;
    S->MaxInteraction       = MaxInteraction;
    S->SizeOfDimention      = SizeOfDimention;
    S->SizeOfFitness        = SizeOfFitness;
    S->accelaration1        = accelaration1;
    S->accelaration2        = accelaration2;
    S->weight               = MAXIMUM + MINIMUM;
    S->HighestNumber        = HighestNumber;
    S->LowestNumber         = LowestNumber;
    S->SizeOfObjective      = SizeOfObjective;
    S->SizeOfArchive        = 10;
    S->Particles            = (Particle*)calloc(SizeOfSwarm,sizeof(Particle));
    Particle *Particles     = S->Particles;

    for(i=0;i<SizeOfSwarm;i++){

        Particles[i].ParticleBestFitness     = (double*)malloc((SizeOfFitness)*sizeof(double));
        Particles[i].CurrentFitness          = (double*)malloc((SizeOfFitness)*sizeof(double));
        Particles[i].Velocity                = (double*)calloc(SizeOfDimention, sizeof(double));
        Particles[i].CurrentPosition         = GeneratePosition(SizeOfDimention,HighestNumber,LowestNumber);
        Particles[i].ParticleBestPosition    = (double*)calloc(SizeOfDimention, sizeof(double));

    }

    for(i=0;i<SizeOfSwarm;i++){

        for(j=0;j<SizeOfFitness;j++){

            Particles[i].ParticleBestFitness[j]     =  INT_MAX;
            Particles[i].CurrentFitness[j]          =  INT_MAX;
        }
    }

    S->ParetoFront = (Archive*)malloc(SizeOfSwarm*sizeof(Archive));

    for(i=0;i<SizeOfSwarm;i++){

        S->ParetoFront[i].CurrentSize           = SIZE_ARCHIVE;
        S->ParetoFront[i].Leader                = (ParetoParticle*)malloc(SIZE_ARCHIVE*sizeof(ParetoParticle));

        for(j=0; j<SIZE_ARCHIVE ;j++){

            S->ParetoFront[i].Leader[j].Fitness     = (double*)malloc(SizeOfFitness*sizeof(double));
            S->ParetoFront[i].Leader[j].Position    = (double*)malloc(SizeOfDimention*sizeof(double));

            for(k=0;k<SizeOfDimention;k++)
                S->ParetoFront[i].Leader[j].Position[k]  = ( (double) rand() / RAND_MAX)*(HighestNumber-LowestNumber) + LowestNumber;

            for(k=0;k<SizeOfFitness;k++)
                S->ParetoFront[i].Leader[j].Fitness[k]  = INT_MAX;
        }

        for(j=0;j<SIZE_ARCHIVE;j++)
            S->ParetoFront[i].Leader[j].Density = 0;

        S->ParetoFront[i].VectorLeader = (int**)malloc(SizeOfObjective*sizeof(int*));

        for(j=0;j<SizeOfObjective;j++)
            S->ParetoFront[i].VectorLeader[j] = (int*)malloc(SIZE_ARCHIVE*sizeof(int));

        for(j=0;j<SizeOfObjective;j++){
            for(k=0;k<SIZE_ARCHIVE;k++)
                S->ParetoFront[i].VectorLeader[j][k] = k;
        }

    }

    ConstructMesh(S);

    return S;

}

