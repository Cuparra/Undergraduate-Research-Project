#include "ParticleSwarmOptimization.h"

#define  ALFA       0.975

#define MINIMUM     0.4
#define MAXIMUM     0.9

#define MAX_SIGMA   0.92
#define MIN_SIGMA   0.02

static inline int CompareMin(double num1, double num2){

    return (num1 < num2) ? 1 : 0;

}

static inline int CompareMax(double num1, double num2){

    return (num1 > num2) ? 1 : 0;

}

static double *GeneratePosition(int SizeOfDimention, double HighestNumber, double LowestNumber){

    double *Position = (double*)calloc(SizeOfDimention, sizeof(double));

    while(SizeOfDimention--)
        Position[SizeOfDimention]   = ( (double) rand() / RAND_MAX)*(HighestNumber-LowestNumber) + LowestNumber;

    return Position;

}

Swarm *InicializeSwarm(int SizeOfSwarm, int SizeOfDimention, int SizeOfFitness, double accelaration1, double accelaration2, double HighestNumber, double LowestNumber,  int MaxInteraction, int type){

    int i,j;
    Swarm *S                = (Swarm*)malloc(sizeof(Swarm));
    S->SizeOfSwarm          = SizeOfSwarm;
    S->MaxInteraction       = MaxInteraction;
    S->SizeOfDimention      = SizeOfDimention;
    S->SizeOfFitness        = SizeOfFitness;
    S->accelaration1        = accelaration1;
    S->accelaration2        = accelaration2;
    S->weight               = MAXIMUM;
    S->sigma                = MAX_SIGMA;
    S->HighestNumber        = HighestNumber;
    S->LowestNumber         = LowestNumber;

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

            Particles[i].ParticleBestFitness[j]     =  (type == MINIMIZATION) ? INT_MAX : INT_MIN;
            Particles[i].CurrentFitness[j]          =  (type == MINIMIZATION) ? INT_MAX : INT_MIN;
        }
    }

    S->ParetoFront                  = (Archive*)malloc(sizeof(Archive));
    S->ParetoFront->SizeOfArchive   = 50;
    S->ParetoFront->CurrentSize     = 50;
    S->ParetoFront->Leader          = (ParetoParticle*)malloc(S->ParetoFront->SizeOfArchive*sizeof(ParetoParticle));

    for(i=0; i<S->ParetoFront->SizeOfArchive ;i++){

        S->ParetoFront->Leader[i].Fitness   = (double*)malloc(SizeOfFitness*sizeof(double));
        S->ParetoFront->Leader[i].Position  = (double*)malloc(SizeOfDimention*sizeof(double));

        for(j=0;j<SizeOfFitness;j++)
            S->ParetoFront->Leader[i].Fitness[j]   = Particles[i].ParticleBestFitness[j];
        for(j=0;j<SizeOfDimention;j++)
            S->ParetoFront->Leader[i].Position[j]  = Particles[i].ParticleBestPosition[j];

    }

    S->ParetoFront->VectorLeader = (int**)malloc(SizeOfFitness*sizeof(int*));
    for(i=0;i<SizeOfFitness;i++)
        S->ParetoFront->VectorLeader[i] = (int*)malloc(S->ParetoFront->SizeOfArchive*sizeof(int));

    for(i=0;i<SizeOfFitness;i++){
        for(j=0;j<S->ParetoFront->SizeOfArchive;j++)
            S->ParetoFront->VectorLeader[i][j] = j;
    }

    if (type == MINIMIZATION)
        S->Compare = &CompareMin;
    else
        S->Compare = &CompareMax;

    return S;

}

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

        int LeaderIdx       = RouletteWheel(A);

        double *velocity    = particle[i].Velocity;
        double *gb          = A->Leader[LeaderIdx].Position;
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
    Particle *particle              = S->Particles;
    int SizeOfFitness               = S->SizeOfFitness;
    int SizeOfDimention             = S->SizeOfDimention;
    int SizeOfSwarm                 = S->SizeOfSwarm;
    int (*Compare)(double,double)   = S->Compare;

    for(i=0;i<SizeOfSwarm;i++){

        int cCurr           = 0;
        int cBest           = 0;
        double *CurrFitness = particle[i].CurrentFitness;
        double *BestFitness = particle[i].ParticleBestFitness;

        for(j=0;j<SizeOfFitness;j++){

            cCurr += Compare(CurrFitness[j],BestFitness[j]) ? 1 : 0;
            cBest += Compare(BestFitness[j],CurrFitness[j]) ? 1 : 0;

        }

        if(cBest == 0 && cCurr > 0){

            for(j=0;j<SizeOfFitness;j++)
                particle[i].ParticleBestFitness[j] = particle[i].CurrentFitness[j];

            CopyPosition(particle[i].ParticleBestPosition, particle[i].CurrentPosition,SizeOfDimention);
        }

    }
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

        UpdateArchive(S);

        /* Update the velocity and position of all particles*/
        UpdateParticles(S);

        /*Update weight and sigma*/
        ++Cycle;
        S->weight  = (MAXIMUM - MINIMUM)*( 1.0*(MaxInteraction - Cycle)/MaxInteraction ) + MINIMUM;

    }while(TerminationCriteria(S));

}
