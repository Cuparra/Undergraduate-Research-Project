#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "ParticleSwarmOptimization.h"
#include "Function.h"

/*Atenção: Colocar os objetivos primeiro, depois a restrição.*/
//modificar o INF para por aproximadamente o menor dos low ou high limit
//Cada particula pode usar gradient descend method

#define SIZE_P      100
#define SIZE_D      10
#define SIZE_F      2
#define SIZE_O      2
#define N_HI        5
#define N_LO        -5
#define SESSION     500

Swarm *S;

int session = 0;

void FON(Particle *P){

    int i;
    double f1 = 0, f2 = 0;
    double *x = P->CurrentPosition;

    double con = 1.0/sqrt(5);

    f1 = 1 - exp( - ((x[0]-con)*(x[0]-con)+(x[1]-con)*(x[1]-con)+(x[2]-con)*(x[2]-con)+(x[3]-con)*(x[3]-con)+(x[4]-con)*(x[4]-con)));
    f2 = 1 - exp( - ((x[0]+con)*(x[0]+con)+(x[1]+con)*(x[1]+con)+(x[2]+con)*(x[2]+con)+(x[3]+con)*(x[3]+con)+(x[4]+con)*(x[4]+con)));

    P->CurrentFitness[0] = f1;
    P->CurrentFitness[1] = f2;

}

void EvaluateParticles(Particle * particles){

    int i;

    for(i=0;i<SIZE_P;i++)
        ZDT4(&particles[i]);

}

int TerminationCriteria(Swarm *S){

    int i,j;

    if(++session < SESSION)
        return 1;

    /*for(i=0;i<SIZE_P;i++){

        Archive *A              = S->ParetoFront;
        ParetoParticle *Leader  = A[i].Leader;
        int CurrentSize         = A[i].CurrentSize;

        for(j=0;j<CurrentSize;j++){
            double *fit = Leader[j].Fitness;
            printf("%lf %lf |%lf\n",fit[0],fit[1],Leader[j].Density);
        }
        printf("\n");
        getchar();
     }*/


    return 0;
}

int main(){

    int i,j;
    time_t t;

    srand((unsigned) time(&t));

    FILE *file = fopen("L_ZDT4","w");

    for(i=0;i<100;i++){

        session = 0;
        printf("%d\n",i);

        S = InicializeSwarm(SIZE_P,5,2,2,2,2,1,0,SESSION);
        SwarmOptimization(S,EvaluateParticles,TerminationCriteria);

        for(j=0;j<50;j++){
            ParetoParticle *p = S->ParetoFront[j].Leader;
            fprintf(file,"%lf %lf\n", p[0].Fitness[0], p[0].Fitness[1]);
        }
    }

    close(file);

    return 0;
}
