#include "Swarm.h"

#define TRUE        1
#define NEUTRAL     0
#define FALSE      -1

#define INF         0.1

#define CalculateTotalDensity(Total)    int tempSize = *CurrentSize;             \
                                        Total = 0;                              \
                                                                                \
                                        while(tempSize--)                       \
                                            Total += Leader[tempSize].Density;


#define Insert(x)       int LeaderIdx;                                       \
                        int include;                                        \
                                                                          \
                        if(*CurrentSize < SizeOfArchive){                      \
                            LeaderIdx   = ++*CurrentSize-1;                     \
                            include     = TRUE;                             \
                        }else{                                              \
                            LeaderIdx   = rand() % (SizeOfArchive-1);               \
                            int density = CalculateDensity(A, &particle[x], SizeOfFitness);       \
                            include     = density > Leader[LeaderIdx].Density;      \
                        }                                                           \
                                                                                        \
                        if(include){                                                  \
                                                                                        \
                            for(k=0;k<SizeOfFitness;k++)                                    \
                                Leader[LeaderIdx].Fitness[k] = particle[x].CurrentFitness[k];   \
                            for(k=0;k<SizeOfDimention;k++)                                 \
                                Leader[LeaderIdx].Position[k] = particle[x].CurrentPosition[k]; \
                            InsertCrowdingSort(A,SizeOfFitness,LeaderIdx);                                          \
                            CrowdingDistance(A, SizeOfFitness);                         \
                            flag = TRUE;                                                  \
                        }                                                                \
                        ++i;


#define Delete(x)       --*CurrentSize;                                     \
                        double *tempPos = Leader[x].Position;                \
                        double *tempFit = Leader[x].Fitness;              \
                                                                        \
                        for(k=x;k<SizeOfArchive-1;k++){                 \
                            Leader[k].Position = Leader[k+1].Position;      \
                            Leader[k].Fitness  = Leader[k+1].Fitness;       \
                        }                                               \
                                                                        \
                        Leader[SizeOfArchive-1].Position = tempPos;              \
                        Leader[SizeOfArchive-1].Fitness  = tempFit;         \
                        DeleteCrowdingSort(A,SizeOfFitness,x);                               \
                        CrowdingDistance(A, SizeOfFitness);                 \
                        flag = TRUE;


#define Dominate(x,y)   int XbestY       = 0;                                    \
                        int YbestX       = 0;                                    \
                        double *Yfitness = Leader[y].Fitness;                     \
                        double *Xfitness = particle[x].CurrentFitness;                     \
                                                                                \
                        for(k=0;k<SizeOfFitness;k++){                           \
                                                                                \
                            XbestY += ( Compare(Xfitness[k], Yfitness[k]) ? 1 : 0 );  \
                            YbestX += ( Compare(Yfitness[k], Xfitness[k]) ? 1 : 0 );  \
                        }                                                       \
                                                                                \
                        if(XbestY > 0 && YbestX == 0)                           \
                            result = TRUE;                                      \
                        else if (YbestX > 0 && XbestY == 0)                     \
                            result = FALSE;                                     \
                        else                                                    \
                            result = NEUTRAL;


#define SwapVec(x,y)        int temp = x; x = y; y = temp;


#define BubbleSort(vec,k)   int change = TRUE;                                                  \
                                                                                                \
                            while(change == TRUE){                                              \
                                                                                                \
                                change = FALSE;                                                 \
                                                                                                \
                                for (i = 0; i < CurrentSize - 1; i++){                          \
                                                                                                \
                                    int idx1 = vec[i];                                          \
                                    int idx2 = vec[i+1];                                        \
                                                                                                \
                                    if (Leader[idx1].Fitness[k] > Leader[idx2].Fitness[k]) {    \
                                        SwapVec(vec[i],vec[i+1]);                               \
                                        change = TRUE;                                          \
                                    }                                                           \
                                }                                                               \
                            }

static inline double Absolute(double num){

    return num >= 0 ? num : -num;

}

int RouletteWheel(Archive *A){

    int num                 = -1;
    double stopPoint        = 0.0;
    double TotalDensity     = A->TotalDensity;
    int CurrentSize         = A->CurrentSize;
    ParetoParticle *Leader  = A->Leader;
    double randomNumber     = ( (double) rand() / RAND_MAX);

    do{
        stopPoint += ( Leader[++num].Density / TotalDensity );
    }while(num < CurrentSize - 1 && stopPoint < randomNumber);

    return num;
}

int TournamentSelection(Archive *A){

    int i;
    ParetoParticle *Leader  = A->Leader;
    int CurrentSize         = A->CurrentSize;
    int best                = rand()%CurrentSize;

    int SizeOfSelection     = 0.4*CurrentSize + 1;

    for(i=0; i<SizeOfSelection; i++){

        int temp = rand()%CurrentSize;

        if(Leader[best].Density < Leader[temp].Density)
            best = temp;
    }

    return best;
}

static double CalculateDensity( Archive *A, Particle *p, int SizeOfFitness ){

    int i;
    int CurrentSize         = A->CurrentSize;
    ParetoParticle *Leader  = A->Leader;
    double density          = 0;

    for(i=0;i<SizeOfFitness;i++){
        int j       = 0;
        int *vec    = &A->VectorLeader[j][0];

        while(j < CurrentSize && p->CurrentFitness[i] > Leader[vec[j]].Fitness[i])
            ++j;

        if(j == 0 || j == CurrentSize)
            density += INF;
        else
            density += Absolute( (Leader[vec[j]].Fitness[i] - Leader[vec[j-1]].Fitness[i])/(Leader[vec[CurrentSize-1]].Fitness[i] - Leader[vec[0]].Fitness[i] + 1) );
    }

    return density;
}

static void CrowdingDistance(Archive *A, int SizeOfFitness){

    int i,j;
    ParetoParticle *Leader  = A->Leader;
    int CurrentSize         = A->CurrentSize;

    for(i=0 ;i<CurrentSize; i++)
        Leader[i].Density = 0;

    for(i=0;i<SizeOfFitness;i++){
        int *vec = &A->VectorLeader[i][0];

        for(j=1;j<CurrentSize-1;j++){

            Leader[vec[j]].Density += Absolute( (Leader[vec[j+1]].Fitness[i] - Leader[vec[j-1]].Fitness[i])/(Leader[vec[CurrentSize-1]].Fitness[i] - Leader[vec[0]].Fitness[i]) );
        }

        Leader[vec[0]].Density              += INF;
        Leader[vec[CurrentSize-1]].Density  += INF;
    }

}

static void InsertCrowdingSort(Archive *A, int SizeOfFitness, int k){

    int i;
    ParetoParticle *Leader  = A->Leader;
    int CurrentSize         = A->CurrentSize;
    int **VectorLeader      = A->VectorLeader;

    if(k == CurrentSize-1){
        for(i=0;i<SizeOfFitness;i++)
            VectorLeader[i][CurrentSize-1] = k;
    }

    for(i=0;i<SizeOfFitness;i++){
        int j=k;
        int *vec = &VectorLeader[i][0];

        while(j>0 && Leader[vec[j]].Fitness[i] < Leader[vec[j-1]].Fitness[i]){
            SwapVec(vec[j],vec[j-1]);
            --j;
        }

        j = k;
        while(j<CurrentSize -1 && Leader[vec[j]].Fitness[i] > Leader[vec[j+1]].Fitness[i]){
            SwapVec(vec[j],vec[j+1]);
            ++j;
        }
    }

}

static void DeleteCrowdingSort(Archive *A, int SizeOfFitness, int k){

    int i;
    int CurrentSize     = A->CurrentSize;
    int **VectorLeader  = A->VectorLeader;

    for(i=0;i<SizeOfFitness;i++){
        int j = 0;
        int *vec = &VectorLeader[i][0];

        while( k!=vec[j] )
            ++j;

        for(j=j;j<CurrentSize;j++)
            vec[j] = vec[j+1];

        for(j=0;j<CurrentSize;j++){
            if(vec[j] > k)
                vec[j] += -1;
        }

    }

}

static void CrowdingSort(Archive *A, int SizeOfFitness ){

    int i,k;
    ParetoParticle *Leader          = A->Leader;
    int CurrentSize                 = A->CurrentSize;
    int **VectorLeader              = A->VectorLeader;

    for(k=0;k<SizeOfFitness;k++){

        int *vec = &VectorLeader[k][0];
        BubbleSort(vec,k);

    }

}

/*Ver isto em baixo*/
void UpdateArchive(Swarm *S){

    int i = 0, j = 0, k;    //i index do S, j index do A
    Archive *A                      = S->ParetoFront;
    ParetoParticle *Leader          = A->Leader;
    int (*Compare)(double,double)   = S->Compare;
    Particle *particle              = S->Particles;
    int *CurrentSize                = &A->CurrentSize;
    int SizeOfArchive               = A->SizeOfArchive;
    int SizeOfFitness               = S->SizeOfFitness;
    int SizeOfDimention             = S->SizeOfDimention;
    int SizeOfSwarm                 = S->SizeOfSwarm;

    int flag                        = FALSE;

    CrowdingSort(A,SizeOfFitness);
    CrowdingDistance(A,SizeOfFitness);

    do{
        int result;
        Dominate(i,j);

        if(result == TRUE){

            Delete(j); //Inclui --CurrSize;

            if(j >= *CurrentSize){
                Insert(i); //Inclui ++i;
            }

        }else if(result == FALSE){
            j = 0;
            ++i;

        }else{
            if (j < *CurrentSize - 1)
                ++j;
            else{
                Insert(i); //Inclui ++i;
                j = 0;
            }
        }

    }while(i < SizeOfSwarm);

    if(flag == TRUE){
        CalculateTotalDensity(A->TotalDensity);
    }

}
