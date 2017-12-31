#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NOT_ACTIVED -1
#define PROCESSING   0
#define ACTIVED      1

#define DISPLACEMENT 2

#define MIN(x,z)    (x > z ? z : x)
#define MAX(x,y)    (x > y ? x : y)

#define ADJUST(x,y,z)   x = MIN(x,z); \
                        x = MAX(x,y); \

#define ADJUST_INPUT(X) X % (MaxNodeInput - 1) + 2

#define OR  0
#define XOR 1
#define AND 2

#define MOD(x,y) (x%y < 0 ? x%y + y : x%y);

int *ExptOuput;
int usedNodes;
int NumFunc;
int NumInput;
int NumOutput;
int NumCol;
int NumRow;
int LevelBack;
int len_string;
int MaxNodeInput;
int TotalNodes; //it includs NumInput and NumOuput

int     *stack;
double  *string;    //size = NumInput + NumOutput + NumCol*NumRow*(MaxNodeInput+2)
int     *actNodes;
int     *Output;
int     *input;

int M, N;

// M = NumInput + NumCol*NumRow;
// N = NumInput + NumCol*NumRow*(MaxNodeInput+2)

// Procedimento: Encontrar os vértices ativos, e depois a cada coluna ,de forma crescente, calcular a função.
// Function | Num of Inputs | Input1 | ... | InputN

int Convert(int NodeId, int Node){

    int NodeId_Row  = (NodeId - NumInput) % NumRow; // get the vectice's row position
    int NodeId_Col  = (NodeId - NumInput) / NumRow; // get the vectice's column position

    int active      = LevelBack <= NodeId_Col;

    int min         = (NodeId - LevelBack*NumRow - NodeId_Row)*active;  // get the minimum vetice allowed
    int max         = NodeId - NodeId_Row - 1 ;                         // get the maximum

    return Node % (max - min + 1) + min*active; // turn the vertice into a allowed one
}

void Decoding(){

    int i=0,k=-1,p=0,q=0;

    for( i=NumInput; i<TotalNodes ;++i )
        actNodes[i] = NOT_ACTIVED;

    for(i=0;i<NumInput;i++)
        actNodes[i] = ACTIVED;

    while( ++k < NumOutput ){ // 5

        q = q + 1;

        stack[p++]      = M + k;
        actNodes[M + k] = ACTIVED;
        int NodeId      = (int) string[N+k];

        ADJUST(NodeId, NumInput, M-1);

        actNodes[NodeId]= PROCESSING;
        stack[p++]      = NodeId;
        string[N+k]     = NodeId;

        while( q != p){ // 4

            int l       = -1;
            int nodeId  = stack[q++];

            if(actNodes[nodeId] == PROCESSING){ // 3

                int strNode         = (nodeId-NumInput)*(MaxNodeInput+DISPLACEMENT) + NumInput + DISPLACEMENT;
                int nodeInput       = ADJUST_INPUT((int) string[strNode - 1]);
                actNodes[nodeId]    = ACTIVED;
                string[strNode - 1] = nodeInput;

                while( ++l < nodeInput ){ // 2

                    int tempNodeId      = Convert(nodeId,(int) string[strNode + l]);
                    string[strNode + l] = tempNodeId;

                    if(actNodes[tempNodeId] == NOT_ACTIVED){ // 1

                        stack[p++]              = tempNodeId;
                        actNodes[tempNodeId]    = PROCESSING;

                    } //1
                } // 2
            } // 3
        } // 4
    } // 5

}

int funcOR(int n){

    int i = 0;

    for(i=0;i<n;i++){
        if(input[i] == 1)
            return 1;
    }

    return 0;
}

int funcAND(int n){

    int i = 0;

    for(i=0;i<n;i++){
        if(input[i] == 0)
            return 0;
    }

    return 1;
}

int funcXOR(int n){

    int i;
    int result = 0;

    for(i=0;i<n;i++)
        result+= input[i];

    if(result == i || result == 0)
        return 0;
    else
        return 1;
}

void Processing(){

    int node;
    int Size_loop = TotalNodes - NumOutput;

    usedNodes = NumInput;

    for(node=NumInput; node<Size_loop ;node++){

        if(actNodes[node] == ACTIVED){

            int k           = -1;
            int j           = -1;
            int strNode     = (node-NumInput)*(MaxNodeInput+DISPLACEMENT) + NumInput +  DISPLACEMENT;
            int nodeInput   = string[strNode - 1];
            int func        = MOD( (int) string[strNode], NumFunc );

            while(++j < nodeInput){ // verifica

                int varId   = string[strNode + j]; // get i's inputs
                input[++k]  = Output[varId];
            }


            if(func == OR)
                Output[node] = funcOR(nodeInput);
            else if(func == XOR)
                Output[node] = funcXOR(nodeInput);
            else
                Output[node] = funcAND(nodeInput);

            ++usedNodes;
        }
    }

    //It puts the result on the outputs
    for(node=Size_loop;node<TotalNodes;node++){

        int nodeId     = string[N+node-Size_loop];
        Output[node]   = Output[nodeId];

    }

}

// Minimize the fitness
double EvaluateFitness(){

    int i;
    double result = 0;

    for(i=0;i<NumOutput;i++){

        int node    = M + i;
        result      += (ExptOuput[i] - Output[node])*(ExptOuput[i] - Output[node]);
        printf("%d %d\n",ExptOuput[i],Output[node]);
    }

    return result;
}

void SetInput(){

    int i;

    for(i=0;i<NumInput;i++){
        string[i] = rand()%2;
        Output[i] = string[i];
        printf("I%d: %d ",i,Output[i]);
   }

   printf("\n");
}

void SetParameters(int nInput, int nOutput, int column, int row, int lb, int maxInput, int nFunc, int *ValueOut){

    int i;
    NumInput        = nInput;
    NumOutput       = nOutput;
    NumCol          = column;
    NumRow          = row;
    LevelBack       = lb;
    MaxNodeInput    = maxInput;
    NumFunc         = nFunc;

    M               = NumInput + NumCol*NumRow;
    TotalNodes      = NumCol*NumRow + NumInput + NumOutput;
    N               = NumInput + NumCol*NumRow*(MaxNodeInput+DISPLACEMENT);
    len_string      = NumCol*NumRow*(MaxNodeInput+DISPLACEMENT) + NumInput + NumOutput;

    string          = (double*)malloc(len_string*sizeof(double));
    actNodes        = (int*)malloc(TotalNodes*sizeof(int));
    stack           = (int*)malloc(TotalNodes*sizeof(int));
    Output          = (int*)malloc(TotalNodes*sizeof(int));
    input           = (int*)malloc((NumInput+DISPLACEMENT)*sizeof(int));
    ExptOuput       = (int*)malloc((nOutput)*sizeof(int));

    for(i=0;i<nOutput;i++)
        ExptOuput[i] = ValueOut[i];
    printf("%d ",ExptOuput[i]);

}

int main(void){

    int i;
    time_t t;
    int ValueOut[3] = {1,0,1};

    srand((unsigned) time(&t));

    SetParameters(6,3,7,3,1,2,3,ValueOut);
    SetInput();

    srand((unsigned) time(&t));
    for(i=NumInput; i<len_string ;i++)
        string[i] = rand()%(TotalNodes+1);

    Decoding();

    Processing();

    printf("%lf\n",EvaluateFitness());

    return 0;
}

/* Debug Functions */
void PrintGraph(){

    int i,j;

    for(i=NumInput; i<TotalNodes - NumOutput ;i++){

        int idx          = (i - NumInput)*(MaxNodeInput+DISPLACEMENT) + NumInput + DISPLACEMENT;
        int NodeInput    = (int) string[idx-1];

        if(actNodes[i] == ACTIVED){

            printf("\n%d <- ", i);

            for(j=0;j<NodeInput;j++){

                int num = (int) string[idx+j];
                printf("%d ", num);
            }

        }
     }

    for(i=0;i<NumOutput;i++){

        int idx = NumInput + NumCol*NumRow*(MaxNodeInput+DISPLACEMENT) + i;
        printf("\nSaida <- %d",(int) string[idx]);

    }

    printf("\n\n");
}

void PrintOutput(int node, int nodeInput, int func){

    int j;

    printf("Node: %d | Func: %d | Result: %d | ",node,func,Output[node]);
    for(j=0;j<nodeInput;j++)
        printf("%d ",input[j]);

    printf("\n");

}

void PrintResults(int node, int Size_loop, int TotalNodes){

    for(node=Size_loop;node<TotalNodes;node++)
        printf("Node: %d | Result: %d\n",node,Output[node]);


}
