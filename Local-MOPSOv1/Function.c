

#include "Function.h"

//  0<= x <= 10, 6 variables,  2 functions, 6 constraints
void OsyczkaFunction(Particle *P){

    double *x = P->CurrentPosition;

    x[2] = x[2]/2.5 + 1;
    x[4] = x[4]/2.5 + 1;
    x[3] = x[3]/1.66666;

    double fit1 = -25*(x[0]-2)*(x[0]-2) - (x[1]-2)*(x[1]-2) - (x[2]-1)*(x[2]-1) - (x[3]-4)*(x[3]-4) - (x[4]-1)*(x[4]-1);

    double fit2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] + x[5]*x[5];

    double rest1 = x[0] + x[1] - 2;
    P->CurrentFitness[2] = rest1 = (rest1 >= 0) ? 0 : rest1*rest1;

    double rest2 = 6 - x[0] - x[1];
    P->CurrentFitness[3] = rest2 = (rest2 >= 0) ? 0 : rest2*rest2;

    double rest3 = 2 - x[1] + x[0];
    P->CurrentFitness[4] = rest3 = (rest3 >= 0) ? 0 : rest3*rest3;

    double rest4 = 2 - x[0] + 3*x[1];
    P->CurrentFitness[5] = rest4 =(rest4 >= 0) ? 0 : rest4*rest4;

    double rest5 = 4 - (x[2]-3)*(x[2]-3) - x[3];
    P->CurrentFitness[6] = rest5 = (rest5 >= 0) ? 0 : rest5*rest5;

    double rest6 = (x[4]-3)*(x[4]-3) + x[5] - 4;
    P->CurrentFitness[7] = rest6 = (rest6 >= 0) ? 0 : rest6*rest6;

    double restriction = rest1+rest2+rest3+rest4+rest5+rest6;

    P->CurrentFitness[0] = (restriction < LIMIT) ? fit1 : 600;
    P->CurrentFitness[1] = (restriction < LIMIT) ? fit2 : 600;

}

//  0<= x <= 1, 30 variables, 2 functions
void ZDT1(Particle *P){

    int i;
    double *x = P->CurrentPosition;
    double f1 = 0, f2 = 0, h = 0, g = 0;

    // F1 function
    f1 = x[0];

    // G function
    for(i=1;i<30;i++)
        g += x[i];
    g = 1 + (9.0/29)*g;

    // H function
    h = 1 - sqrt(f1/g);

    // F2 function
    f2 = g*h;

    P->CurrentFitness[0] = f1;
    P->CurrentFitness[1] = f2;
}

//  0<= x <= 1, 30 variables, 2 functions
void ZDT2(Particle *P){

    int i;
    double *x = P->CurrentPosition;
    double f1 = 0, f2 = 0, h = 0, g = 0;

    // F1 function
    f1 = x[0];

    // G function
    for(i=1;i<30;i++)
        g += x[i];
    g = 1 + (9.0/29)*g;

    // H function
    h = 1 - (x[0]/g)*(x[0]/g);

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
    for(i=1;i<30;i++)
        g += x[i];
    g = 1 + (9.0/29)*g;

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

// -4 <= x <= 4, 30 variables, 2 functions
void FonFleFunction(Particle *P){

    int i;
    double f1 = 0, f2 = 0;
    double *x = P->CurrentPosition;

    // F1 function
    for(i=0;i<30;i++)
        f1 += (x[i] - CONSTANT) * (x[i] - CONSTANT);

    f1 = 1 - exp(-f1);

    // F2 function
    for(i=0;i<30;i++)
        f1 += (x[i] + CONSTANT) * (x[i] + CONSTANT);

    f1 = 1 - exp(-f1);

    P->CurrentFitness[0] = f1;
    P->CurrentFitness[1] = f2;
}

// -7 <= x <= 4, 2 variables, 2 functions, 3 constraints
void TestFunction4(Particle *P){

    double f1 = 0, f2 = 0;
    double restriction = 0;
    double *x = P->CurrentPosition;
    double rest1 = 0, rest2 = 0, rest3 = 0;

    f1      = x[0]*x[0] - x[1];
    f2      = -0.5*x[0] - x[1] - 1;

    rest1   = 6.5 - (x[0]/6.0) - x[1];
    rest2   = 7.5 - 0.5*x[0] - x[1];
    rest3   = 30 - 5*x[0] - x[1];

    P->CurrentFitness[2] = ( rest1 = (rest1 >= 0 ? 0 : rest1*rest1) );
    P->CurrentFitness[3] = ( rest2 = (rest2 >= 0 ? 0 : rest2*rest2) );
    P->CurrentFitness[4] = ( rest3 = (rest3 >= 0 ? 0 : rest3*rest3) );

    restriction = rest1 + rest2 +rest3;

    P->CurrentFitness[0] = (restriction < 0.00001) ? f1 : 600;
    P->CurrentFitness[1] = (restriction < 0.00001) ? f2 : 600;
}

// 0 <= x <= 1, 2 variables, 2 functions, 2 constraints
void CTP1Function(Particle *P){

    double f1 = 0, f2 = 0;
    double restriction = 0;
    double rest1 = 0, rest2 = 0;
    double *x = P->CurrentPosition;

    f1      = x[0];
    f2      = (1 + x[1]) * exp( - ( x[0]/(1+x[1]) ) );

    rest1   =  f2 / ( 0.858 * exp(-0.541*f1) );
    rest1   =  f2 / ( 0.728 * exp(-0.295*f1) );

    P->CurrentFitness[2] = ( rest1 = (rest1 >= 1 ? 0 : rest1*rest1) );
    P->CurrentFitness[3] = ( rest2 = (rest2 >= 1 ? 0 : rest2*rest2) );

    restriction = rest1 + rest2;

    P->CurrentFitness[0] = (restriction < 0.00001) ? f1 : 600;
    P->CurrentFitness[1] = (restriction < 0.00001) ? f2 : 600;
}

// 0 <= x <= 5, 2 variables, 2 functions, 2 constraints
void BihnKornFunction(Particle *P){

    double f1 = 0, f2 = 0;
    double restriction = 0;
    double rest1 = 0, rest2 = 0;
    double *x = P->CurrentPosition;

    x[1]    =  3*(x[1]/5);

    f1      =  4*x[0]*x[0] + 4*x[1]*x[1];
    f2      =  (x[0] - 5)*(x[0] - 5) + (x[1] - 5)*(x[1] - 5);

    rest1   =  (x[0] - 5)*(x[0] - 5) - x[1]*x[1];
    rest2   =  (x[0] - 8)*(x[0] - 8) + (x[1] + 3)*(x[1] + 3);

    P->CurrentFitness[2] = ( rest1 = (rest1 <= 25  ? 0 : rest1*rest1) );
    P->CurrentFitness[3] = ( rest2 = (rest2 >= 7.7 ? 0 : rest2*rest2) );

    restriction = rest1 + rest2;

    P->CurrentFitness[0] = (restriction < 0.00001) ? f1 : 600;
    P->CurrentFitness[1] = (restriction < 0.00001) ? f2 : 600;
}

// -M_PI <= x <= M_PI, 2 variable, 2 functions
void PoloniFunction(Particle *P){

    double f1 = 0, f2 = 0;
    double *x = P->CurrentPosition;

    double B1   = 0.5*sin(x[0]) - 2*cos(x[0]) + sin(x[1]) - 1.5*cos(x[1]);
    double B2   = 1.5*sin(x[0]) - cos(x[0]) + 2*sin(x[1]) - 0.5*cos(x[1]);

    f1          = 1 + (A1 - B1)*(A1 - B1) + (A2 - B2)*(A2 - B2);
    f2          = (x[0] + 3)*(x[0] + 3) + (x[1] + 1)*(x[1] + 1);

    P->CurrentFitness[0] = f1;
    P->CurrentFitness[1] = f2;
}

double ErrorRation(Swarm *S, FILE* file){

    int i,j;
    double error  = 0;
    int sizeOfSw  = S->SizeOfSwarm;
    int sizeOfAr  = S->SizeOfArchive;
    int sizeOfSol = sizeOfSw*sizeOfAr;

    for(i=0;i<sizeOfSw;i++){

        ParetoParticle *p = S->ParetoFront[i].Leader;

        for(j=0;j<sizeOfAr;j++){

            double *sol = p[j].Fitness;

            double result = (1 - sqrt(sol[0]) - sol[0]*sin(10*M_PI*sol[0]) -sol[1])*(1 - sqrt(sol[0]) - sol[0]*sin(10*M_PI*sol[0]) -sol[1]);

            if(result <= 0.000001){

                if(sol[0] >= 0 && sol[0] <= 0.0830015349)
                    ++error;
                else if(sol[0] >= 0.4093136748 && sol[0]<=0.4538821041)
                    ++error;
                else if(sol[0] >= 0.6183967944 && sol[0]<=0.6525117038)
                    ++error;
                else if(sol[0] >= 0.8233317983 && sol[0]<=0.8518328654)
                    ++error;
            }

            //else
            //printf("%d-%lf %lf | %lf\n",++f,sol[0],sol[1], 1 - sqrt(sol[0]) );
        }
    }

    error /= sizeOfSol;
    fprintf(file,"%lf\n",error);
    //getchar();
    return error;
}

