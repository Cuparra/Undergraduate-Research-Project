#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "ParticleSwarmOptimization.h"

#define SIZE_P      100
#define SIZE_D      30
#define SIZE_O      2
#define SIZE_F      2
#define N_HI        1
#define N_LO        0
#define SESSION     100

// Comparar densidade no best position

Swarm *S;

int session = 0;

void ZDT1(Particle *P){

    int i;
    double *x = P->CurrentPosition;
    double f1 = 0, f2 = 0, h = 0, g = 0;

    // F1 function
    f1 = x[0];

    // G function
    for(i=1;i<SIZE_D;i++)
        g += x[i];
    g = 1 + (9.0/SIZE_D-1)*g;

    // H function
    h = 1 - sqrt(f1/g);

    // F2 function
    f2 = g*h;

    P->CurrentFitness[0] = f1;
    P->CurrentFitness[1] = f2;
}

void EvaluateParticles(Particle * particles){

    int i;

    for(i=0;i<SIZE_P;i++)
        ZDT1(&particles[i]);

}

int TerminationCriteria(Swarm *S){
    int i,j;

    Archive *A                      = S->ParetoFront;
    ParetoParticle *Leader          = A->Leader;

    printf("%d\n",++session);

    if(session < SESSION)
        return 1;

    for(i=0;i<A->CurrentSize;i++){
        double *fit = Leader[i].Fitness;
        printf("%lf %lf| %lf %lf %lf %lf %lf %lf\n",fit[0],fit[1]);
    }

    getchar();

    return 0;
}

int main(){

    time_t t;
    srand((unsigned) time(&t));

    S = InicializeSwarm(SIZE_P,SIZE_D,SIZE_O,SIZE_F,2,2,N_HI,N_LO,SESSION);

    SwarmOptimization(S,EvaluateParticles,TerminationCriteria);

    return 0;
}
