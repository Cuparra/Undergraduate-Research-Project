#include "ParticleSwarmOptimization.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

static inline void CopyPosition(double *Pos1, double *Pos2, int Size){

    while(Size--)
        Pos1[Size] = Pos2[Size];

}

/* Update the velocity and position*/
static void UpdateParticles(Swarm *S){

    int i,j;

    double w                = S->weight;
    Particle *particle      = S->Particles;
    Archive *A              = S->ParetoFront;
    double c1               = S->accelaration1;
    double c2               = S->accelaration2;
    double HighestNumber    = S->HighestNumber;
    double LowestNumber     = S->LowestNumber;
    int SizeOfDimention     = S->SizeOfDimention;
    int SizeOfSwarm         = S->SizeOfSwarm;

    for(i=0;i<SizeOfSwarm;i++){

        int LeaderIdx       = TournamentSelection(&A[i], S->SizeOfFitness);
        double *velocity    = particle[i].Velocity;
        double *gb          = A[i].Leader[LeaderIdx].Position;
        double *x           = particle[i].CurrentPosition;
        double *pb          = particle[i].ParticleBestPosition;

        for(j=0;j<SizeOfDimention;j++){

            double r1       = ( (double) rand()/RAND_MAX );
            double r2       = ( (double) rand()/RAND_MAX );

            velocity[j]     = w*velocity[j] + c1*r1*(pb[j] - x[j]) + c2*r2*(gb[j] - x[j]);

            x[j]            = x[j] + velocity[j];

            if(x[j] > HighestNumber)
                x[j] = HighestNumber;
            if(x[j] < LowestNumber)
                x[j] = LowestNumber;
        }
    }

}

static void SetParticleBestPosition(Swarm *S){

    int i,j;
    Archive *A              = S->ParetoFront;
    Particle *particle      = S->Particles;
    int SizeOfObjective     = S->SizeOfObjective;
    int SizeOfFitness       = S->SizeOfFitness;
    int SizeOfDimention     = S->SizeOfDimention;
    int SizeOfSwarm         = S->SizeOfSwarm;

    for(i=0;i<SizeOfSwarm;i++){

        int cCurr           = 0;
        int cBest           = 0;
        double *CurrFitness = particle[i].CurrentFitness;
        double *BestFitness = particle[i].ParticleBestFitness;

        for(j=0;j<SizeOfFitness;j++){

            cCurr += (CurrFitness[j] < BestFitness[j]);
            cBest += (BestFitness[j] < CurrFitness[j]);

        }

        double densityCurr = CalculateDensity(&A[i],particle[i].CurrentFitness,SizeOfObjective);
        double densityBest = CalculateDensity(&A[i],particle[i].ParticleBestFitness,SizeOfObjective);

        if(cCurr > 0 && cBest >= 0 && densityCurr > densityBest ){

            for(j=0;j<SizeOfFitness;j++)
                particle[i].ParticleBestFitness[j] = particle[i].CurrentFitness[j];

            CopyPosition(particle[i].ParticleBestPosition, particle[i].CurrentPosition,SizeOfDimention);
        }
    }
}

static void SetNeighborhoodBestPosition(Swarm *S){

    int i;
    Archive *A      = S->ParetoFront;
    int SizeOfSwarm = S->SizeOfSwarm;

    for(i=0;i<SizeOfSwarm;i++)
        UpdateArchive(S,&A[i]);

}

void SwarmOptimization(Swarm *S, void (*EvaluateParticles)(Particle *), int (*TerminationCriteria)(Swarm *)){
	
    int Cycle		    = 0;
    int MaxInteraction  = S->MaxInteraction;
    Particle *Particles = S->Particles;

    do{

        /*Set the fitness of all particles*/
        EvaluateParticles(Particles);

        /*Set particle's best position*/
        SetParticleBestPosition(S);

        /*Set the neighborhood(local) best position */
        SetNeighborhoodBestPosition(S);

        /* Update the velocity and position of all particles*/
        UpdateParticles(S);

        /*Update weight and sigma*/
        ++Cycle;
        S->weight  = MAXIMUM*( 1.0*(MaxInteraction - Cycle)/MaxInteraction ) + MINIMUM;

    }while(TerminationCriteria(S));

}
