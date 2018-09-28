
#include "Function.h"


//  0<= x <= 1, 100 variables, 2 functions
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

//  -5 <= x <= 5, 10 variables, 4 functions
void ZDT4(Particle *P){

    int i;
    double *x = P->CurrentPosition;
    double f1 = 0, f2 = 0, h = 0, g = 0;

    // F1 function
    f1 = 1.0*(x[0]+5.0)/10.0;

    // G function
    for(i=1;i<SIZE_D;i++)
        g += (x[i]*x[i] - 10.0*cos(4.0*M_PI*x[i]));
    g += 1 + 10*(SIZE_D-1);

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


double ErrorZDT3(Swarm *S, FILE* file){

    int i,j;
    double error  = 0;
    int sizeOfSw  = S->SizeOfSwarm;
    int sizeOfAr  = S->SizeOfArchive;
    int sizeOfSol = sizeOfSw*sizeOfAr;

    for(i=0;i<sizeOfSw;i++){

        ParetoParticle *p = S->ParetoFront[i].Leader;

        for(j=0;j<sizeOfAr;j++){

            double *sol = p[j].Fitness;

            double result = (1 - sqrt(sol[0]) - sol[1])*(1 - sqrt(sol[0]) - sol[1]);

            if(result <= 0.000001)
                ++error;

        }
    }

    error /= sizeOfSol;
    fprintf(file,"%lf\n",error);
    return error;
}

double ErrorZDT6(Swarm *S, FILE* file){

    int i,j;
    double error  = 0;
    int sizeOfSw  = S->SizeOfSwarm;
    int sizeOfAr  = S->SizeOfArchive;
    int sizeOfSol = sizeOfSw*sizeOfAr;

    for(i=0;i<sizeOfSw;i++){

        ParetoParticle *p = S->ParetoFront[i].Leader;

        for(j=0;j<sizeOfAr;j++){

            double *sol = p[j].Fitness;

            double result = (1 - sol[0]*sol[0] - sol[1])*(1 - sol[0]*sol[0] - sol[1]);

            if(result <= 0.1)
                ++error;

        }
    }

    error /= sizeOfSol;
    fprintf(file,"%lf\n",error);
    return error;
}

