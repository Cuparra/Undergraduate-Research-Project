#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "ParticleSwarmOptimization.h"
#include "Inicialization.h"
#include "Function.h"

/*Atenção: Colocar os objetivos primeiro, depois a restrição.*/

#define SIZE_P      100
#define SIZE_F      2
#define SIZE_O      2
#define N_HI        1
#define N_LO        0
#define SESSION     1000

Swarm *S;

int session = 0;

void EvaluateParticles(Particle * particles){

    int i;

    for(i=0;i<100;i++)
        ZDT1(&particles[i]);

}


int TerminationCriteria(Swarm *S){

    int i,j;

    if(++session < SESSION)
        return 1;


    for(i=0;i<SIZE_P;i++){

        Archive *A              = S->ParetoFront;
        ParetoParticle *Leader  = A[i].Leader;
        int CurrentSize         = A[i].CurrentSize;

        for(j=0;j<CurrentSize;j++){
            double *fit = Leader[j].Fitness;
            printf("%d %lf %lf |%lf\n",j,fit[0],fit[1]);

        }
        printf("\n");
        getchar();
     }

    return 0;
}


int main(){

    int i,k,j;
    time_t t;

    srand((unsigned) time(&t));

    FILE *file = fopen("ZDT2","w");

    for(i=0;i<2;i++){
        session = 0;

        S = InicializeSwarm(SIZE_P,30,2,2,2,2,1,0,SESSION);
        SwarmOptimization(S,EvaluateParticles,TerminationCriteria);

        for(j=0;j<5;j++){
            ParetoParticle *p = S->ParetoFront[j].Leader;

            for(k=0;k<10;k++)
                fprintf(file,"%lf %lf\n", p[k].Fitness[0], p[k].Fitness[1]);
        }
    }

    return 0;
}
