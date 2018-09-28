#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "ParticleSwarmOptimization.h"

#define SIZE_P      100
#define SIZE_D      10
#define SIZE_O      2
#define N_HI        5
#define N_LO        -5
#define SESSION     500

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

//  0<= x <= 1, 100 variables, 2 functions
void ZDT2(Particle *P){

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
    h = 1 - (f1/g)*(f1/g);

    // F2 function
    f2 = g*h;

    P->CurrentFitness[0] = f1;
    P->CurrentFitness[1] = f2;
}

//  0<= x <= 1, 30 variables, 2 functions
void ZDT3(Particle *P){

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
    h = 1 - sqrt(f1/g) - (f1/g)*sin(10*M_PI*f1);

    // F2 function
    f2 = g*h;

    P->CurrentFitness[0] = f1;
    P->CurrentFitness[1] = f2;
}

void ZDT4(Particle *P){

    int i;
    double *x = P->CurrentPosition;
    double f1 = 0, f2 = 0, h = 0, g = 91;

    // changing the variable's interval...
    x[0] = (x[0]+5)/10;

    // F1 function
    f1 = x[0];

    // G function
    for(i=1;i<10;i++)
        g += (x[i]*x[i] - 10.0*cos(4.0*M_PI*x[i]));

    // H function
    h = 1 - sqrt(f1/g);

    // F2 function
    f2 = g*h;

    P->CurrentFitness[0] = f1;
    P->CurrentFitness[1] = f2;
}

//  0 <= x <= 1, 10 variables, 4 functions
void ZDT6(Particle *P){

    int i;
    double *x = P->CurrentPosition;
    double f1 = 0, f2 = 0, h = 0, g = 0;

    // F1 function
    f1 = 1 - exp(-4*x[0]) * pow ( sin(6*M_PI*x[0]) , 6 );

    // G function
    for(i=1;i<10;i++)
        g += x[i];
    g = 1 + 9 * pow( g/9 , 0.25 );
    // H function
    h = 1 - (f1/g) * (f1/g);

    // F2 function
    f2 = g*h;

    P->CurrentFitness[0] = f1;
    P->CurrentFitness[1] = f2;
}

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

    if(session++ < SESSION)
        return 1;

    /*
    Archive *A                      = S->ParetoFront;
    ParetoParticle *Leader          = A->Leader;

    for(i=0;i<A->CurrentSize;i++){
        double *fit = Leader[i].Fitness;
        printf("%lf %lf\n",fit[0],fit[1]);
    }*/

    return 0;
}

int main(){

    int i,k,j;
    time_t t;

    srand((unsigned) time(&t));

    FILE *file = fopen("ZDT4","w");
    //FILE *file1 = fopen("FON1","w");

    for(i=0;i<100;i++){

        printf("%d\n",i);

        session = 0;

        S = InicializeSwarm(SIZE_P,SIZE_D,SIZE_O,2,2,N_HI,N_LO,SESSION,MINIMIZATION);
        SwarmOptimization(S,EvaluateParticles,TerminationCriteria);

        Archive *A                      = S->ParetoFront;
        ParetoParticle *Leader          = A->Leader;

        for(k=0;k<50;k++)
            fprintf(file,"%lf %lf\n", Leader[k].Fitness[0], Leader[k].Fitness[1]);

      }

    fclose(file);
    //fclose(file1);

    return 0;
}
