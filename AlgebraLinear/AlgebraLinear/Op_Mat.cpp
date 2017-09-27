#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <time.h>

void preenche(int l, int c, float **mat)
{
    for (int i = 0; i < l; i++) {
        for (int j = 0; j < c; j++) {
            printf("\n Digite a posicao [%d][%d] da matriz: ", i+1, j+1);
            scanf("%f", &mat[i][j]);
        }
    }
}
void preencheVet(int tam, float *vet)
{
    printf("\n");
    for (int i = 0; i < tam; i++) {
            printf("\n Digite a posicao [%d] do vetor: ", i+1);
            scanf("%f", &vet[i]);
        }
  
}
void impVet(int tam, float *vet)
{
    for (int i = 0; i < tam; i++) {
            printf("-[%f]-\n", vet[i]);
        }
  
}
float numAl(int Maior, int Menor){
      //srand((unsigned)time(0)); //para gerar números aleatórios reais.
      int aleatorio = rand()%(Maior-Menor+1) + Menor;
      return aleatorio;
}
float** copia(int M, int N, float **matA){
    float** matB;
    matB = (float **)malloc(sizeof(float)*M);
    for (int i = 0; i < M; i++)
        matB[i] = (float *)malloc(sizeof(float)*(N));

    for(int i = 0; i<M;i++)
        for(int j=0;j<N;j++)
            matB[i][j]=matA[i][j];
   
    return matB;

}
void preencheAleatorio(int l, int c, float **mat, int Max, int Min)
{
    for (int i = 0; i < l; i++)
        for (int j = 0; j < c; j++)
            mat[i][j]= numAl(Max,Min);
      
}
void imp(int l, int c, float **mat){
  
    for (int i = 0; i < l; i++){
        printf("\n");
        for (int j = 0; j < c; j++)
            printf("-[%f]-", mat[i][j]);
    }
    printf("\n");
}
float** mult(int M, int N, int P, float **matA, float **matB)
{
    float **matC;

    matC = (float **)malloc(sizeof(float)*M);
    for (int i = 0; i < M; i++)
        matC[i] = (float *)malloc(sizeof(float)*P);


  

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < P; j++) {
          
            matC[i][j] = 0;
            for (int k = 0; k < N; k++) {
                matC[i][j] += matA[i][k] * matB[k][j];
              
            }
            //printf("\ni = %d; j=%d; C[i][j] = %f", i,j, matC[i][j]);
        }
      
    }

    return matC;
}
float** MatIdentidade(int N)
{
    float **mat;
    mat = (float **)malloc(sizeof(float)*N);
    for (int i = 0; i < N; i++)
        mat[i] = (float *)malloc(sizeof(float)*N);

    for(int i =0; i< N;i++)
        for(int j =0; j < N; j++){
            if(i==j){
                mat[i][j]=1;
            }else{
                mat[i][j]=0;
            }

        }

  

    return mat;
}
float** Transposta (int l, int c, float** A) {
  int i, j;
 
  float **matA;

    matA = (float **)malloc(sizeof(float)*l);
    for (int i = 0; i < l; i++)
        matA[i] = (float *)malloc(sizeof(float)*(c));
  
  
  for(int i=0; i<l;i++){
      for(int j=0; j<c;j++){
          matA[i][j]=A[j][i];
      }
  }
  return matA;
}
float** TrocaLinha(int N, int L1, int L2, float **matA){
    float** Id, **Result;
    Id = MatIdentidade(N);
    Id[L1][L1] = 0;
    Id[L1][L2] = 1;
    Id[L2][L2] = 0;
    Id[L2][L1] = 1;

    imp(N,N,Id);

    Result = mult(N,N,N,Id,matA);
    return Result;
}
float** TrocaColuna (int M,int N, int P, int C1, int C2, float **matA){
    float** Id, **Result;
    Id = MatIdentidade(N);
    Id[C1][C1] = 0;
    Id[C1][C2] = 1;
    Id[C2][C2] = 0;
    Id[C2][C1] = 1;
    Result = mult(M,N,P,matA,Id);
    return Result;
}
float* Gauss (int N, float **matA, float *vet2){
  
    float* vet1 = (float *)malloc(sizeof(float)*N);
   
    //Vetor auxiliar para alterar índices de incognitas em caso de pivotação
    int* vetAux = (int *)malloc(sizeof(int)*N);
    for(int k=0;k<N;k++){
        vetAux[k] = k;
    }
    int i,j;
    float alfa;
    //Primeiro Passo
    for(j=0;j<=(N-2);j++){
        for(i=(j+1);i < N;i++){
            //printf("\nValor j = %d, i = %d", j,i);
            if(matA[i][j]!=0){
                if(matA[j][j]==0){ //Pivotação
                    float MJ = 0;
                    int aux=0;
                   
                    for(int k=j+1; k<N;k++){
                        if(matA[j][k]>MJ){
                            MJ = matA[j][k];
                            aux = k;
                            }
                    }
                    int NC = aux; //Nova Coluna
                    matA = TrocaColuna(N,N,N,j,NC,matA);
                    //imp(N,N,matA);
                    vetAux[j]=NC;
                    vetAux[NC]=j;
                }
                alfa = (-matA[i][j])/matA[j][j];
                matA[i][j] = 0;
                for(int k=(j+1); k < N; k++){
                    matA[i][k] = matA[i][k] + (matA[j][k]*alfa);
                }
               
                vet2[i] = vet2[i] + (vet2[j]*alfa);
                //printf("\nAlfa = %f; Valor A[j][j] = %f; valor B[i] = %f",alfa, matA[j][j],vet2[i]);
            }
        }
    }
    //Resubstituição
   //imp(N,N,matA);
   //impVet(N,vet2);
    for(i = (N-1); i >=0; i--){
        float soma = 0;
        for(int k = i+1; k < N; k++)
            soma += matA[i][k] * vet2[k];
        vet2[i] = (vet2[i]-soma)/matA[i][i];
      
    }

    //impVet(N,vet2);

    for(int k=0;k<N;k++){
        vet1[k]=vet2[vetAux[k]];
    }
  
    return vet1;
}
float* GaussJordan(int N, float **matA, float *vet2){
  
    //Vetor Referência p/ Pivotacao
  
    int* vetAux = (int *)malloc(sizeof(int)*N);
    for(int k=0;k<N;k++){
        vetAux[k] = k;
    }
    //Fim Vetor Referencia

  
    for(int j = 0; j < N;j++){
        matA[j][N] = vet2[j];
    }//Loop para concatenar matriz e vetor

    //imp(N,N+1,matA);
    //Primeira Parte (Etapa Gauss)

    for(int j = 0; j < N; j ++){

        //Fazer pivotacao;
        if(matA[j][j]==0 && j<=N-2){ //Pivotação
                float MJ = 0;
                int aux=0;
                   
                for(int k=j+1; k<N;k++){
                    if(matA[j][k]>MJ){
                        MJ = matA[j][k];
                        aux = k;
                    }
                }
                int NC = aux;  
                //Nova Coluna
                matA = TrocaColuna(N,(N+1),(N+1),j,NC,matA);
                vetAux[j]=NC;
                vetAux[NC]=j;
              
        }

        for(int k = (j+1); k<=N; k++){
            matA[j][k]=matA[j][k]/matA[j][j];
        }
        matA[j][j]=1;

        for(int i=(j+1);i<N;i++){
            for(int k=j+1;k<=N;k++){
                matA[i][k]=matA[i][k]-(matA[i][j]*matA[j][k]);
              
            }
            matA[i][j]=0;
        }

    }
    // FIM Primeira Parte
    //imp(N,N+1,matA);

    //Segunda parte:

    for(int j=(N-1);j>=1;j--){
        for(int i=(j-1);i>=0;i--){
            for(int k = N; k>i;k--){
                matA[i][k] = matA[i][k] - (matA[i][j]*matA[j][k]);
            }
        }
    }
  
    //FIM segunda parte;

    //Atribuição de resultados ao vetor;
    for(int j=0; j<N;j++){
        vet2[j]=matA[j][N];
    }

    //Reordenação do Vetor de Resultados no caso de Pivotação
    float* Result;
    Result = (float *)malloc(sizeof(float)*N);

    for(int k=0;k<N;k++){
        Result[k]=vet2[vetAux[k]];
    }
  

    return Result;
}
float** InversaGJ(int N, float **matA){
   
    float** Inv, **id;
    id=MatIdentidade(N);

    Inv = (float **)malloc(sizeof(float)*N);
    for (int i = 0; i < (N*2); i++)
        Inv[i] = (float *)malloc(sizeof(float)*(N));
    
    for(int i = 0; i < N;i++){
        for(int j=0;j<(N*2);j++){
            if(j<N)
                Inv[i][j]=matA[i][j];
            else
                Inv[i][j]=id[i][j-N];
        }
    }//Loop para concatenar matrizes

    imp(N,N*2,Inv);
    //Primeira Parte (Etapa Gauss)

    for(int j = 0; j < N; j ++){

        for(int k = (j+1); k<(N*2); k++){
            Inv[j][k]=Inv[j][k]/Inv[j][j];
        }
        Inv[j][j]=1;

        for(int i=(j+1);i<N;i++){
            for(int k=j+1;k<(N*2);k++){
                Inv[i][k]=Inv[i][k]-(Inv[i][j]*Inv[j][k]);
              
            }
            Inv[i][j]=0;
        }

    }
    // FIM Primeira Parte
    //imp(N,N+1,matA);
   
    //Segunda parte:

    for(int j=(N-1);j>=1;j--){
        for(int i=(j-1);i>=0;i--){
            for(int k = N; k>i;k--){
                Inv[i][k] = Inv[i][k] - (Inv[i][j]*Inv[j][k]);
            }
        }
    }
    imp(N,N*2,Inv);
    return Inv;
}
float** DecLU(int N, float **A){
   
    for(int j=0; j<N;j++){
   
        for(int i=0;i<=j;i++){
            float soma = 0;
            for(int k=0;k<=(i-1);k++)
                soma += A[i][k] * A[k][j];
            A[i][j] = A[i][j]-soma;
        }
        for(int i=(j+1);i<N;i++){
            float soma = 0;
            for(int k=0;k<=(j-1);k++)
                soma += A[i][k]*A[k][j];
            A[i][j]=(A[i][j]-soma)/A[j][j];
        }

    }

    return A;
}
float** DecLU2(int N, float **A){
    for(int j=0;j<=(N-2);j++){
        float alfa;
        for(int i = (j+1);i < N; i++){
            alfa = (A[i][j])/A[j][j];
            A[i][j] = alfa;
            for(int k=(j+1); k < N; k++)
                A[i][k] = A[i][k] - (A[j][k]*alfa);
        }
    }
    return A;
}
float* SELDecLU(int N, float **matA, float *vet1){

    int i,j;
    float alfa;
    float** L;
    float** U;
    U = DecLU(N,matA);
    L = MatIdentidade(N);

    for(i=1;i<N;i++)
        for(j=0;j<i;j++){
            L[i][j]=U[i][j];
            U[i][j]=0;
        }

    //RetroSubstituição com L para encontrar vetor intermediário/
   
    for(i = 1; i < N; i++){
        float soma = 0;
        for(int k = 0; k < i; k++)
            soma += L[i][k] * vet1[k];
        vet1[i] = (vet1[i]-soma);   
               
    }
    //RetroSubstituição com U para encontrar vetor de Incognitas
   
    for(i = (N-1); i >=0; i--){
        float soma = 0;
        for(int k = i+1; k < N; k++)
            soma += U[i][k] * vet1[k];
        vet1[i] = (vet1[i]-soma)/U[i][i];
      
    }

    return vet1;
}
float** Cholesky(int N, float** A){
    float** S = MatIdentidade(N);
    for(int j = 0; j<N; j++){
        float soma;
        for(int i=0;i<j;i++){
            soma=0;
            for(int k=0;k<i;k++)
                soma+=S[i][k]*S[j][k];
            S[j][i]=(A[i][j]-soma)/S[i][i];
        }
        soma=0;
        for(int k=0;k<j;k++)
            soma+=S[j][k]*S[j][k];
        float temp = A[j][j]-soma;
        if(temp>0)
            S[j][j]= sqrt(temp);   
        else
            printf("\nErro em S[%d][%d] -> (Ajj - soma) < 0. A matriz nao eh definida positiva.\n",j+1,j+1);
    }
    return S;
}
float* SELCholensky(int N, float** A, float *B){
   
    A = Cholesky(N,A);
    float** S2 = Transposta(N,N,A);

    //RetroSubstituição com S para encontrar vetor intermediário/

    for(int i = 0; i < N; i++){
        float soma = 0;
        for(int k = 0; k < i; k++)
            soma += A[i][k] * B[k];
        B[i] = (B[i]-soma)/A[i][i];   
               
    }

    //RetroSubstituição com U para encontrar vetor de Incognitas
   
    for(int i = (N-1); i >=0; i--){
        float soma = 0;
        for(int k = i+1; k < N; k++)
            soma += S2[i][k] * B[k];
        B[i] = (B[i]-soma)/S2[i][i];
      
    }
   
    return B;
}
bool igual(int M,int N, float** LU, float** A){
    bool ret = true;
    for(int i=0; i<M;i++)
        for(int j=0;j<N;j++)
            if(LU[i][j]!=A[i][j])
                ret = false;
    return ret;
}
bool ConfereLU(int N, float** LU, float** A){
    float** L = MatIdentidade(N);
    float** U = copia(N,N,LU);
    for(int i=1;i<N;i++)
        for(int j=0;j<i;j++){
            L[i][j]=U[i][j];
            U[i][j]=0;
        }
   
    float** C = mult(N,N,N,L,U);
    bool Ret = igual(N,N,C,A);

    return Ret;
}
int main()
{
   
    int N = 6;

    float **matA;

    matA = (float **)malloc(sizeof(float)*N);
    for (int i = 0; i < N; i++)
        matA[i] = (float *)malloc(sizeof(float)*(N));
  
   
    //Sistemas de equações lineares:
    /*
    preenche(N, N, matA);
    float *vetX,*vetB;
    vetB = (float *)malloc(sizeof(float)*N);
    vetX = (float *)malloc(sizeof(float)*N);

    preencheVet(N,vetB);
   
   
    //vetX = Gauss(N,matA,vetB);
    //vetX = GaussJordan(N, matA, vetB);
    //vetX = SELDecLU(N,matA,vetB);
    //vetX = SELCholensky(N,matA,vetB);

    impVet(N,vetX);
    */
   
    //preenche(N, N, matA);
    //InversaGJ(N,matA);
   
    //Teste LU
    /*
    preenche(N, N, matA);
    float **matB = copia(N,N,matA);
    matA = DecLU(N,matA);
    if(ConfereLU(N,matA,matB))
        imp(N,N,matA);
    else
        printf("\n Erro durante a decomposição, L*U != A\n");
    */
	

   
    //Cholesky
    
    preenche(N,N,matA);
    float ** S = Cholesky(N,matA);
    float ** ST = Transposta(N,N,S);
    float** multS= mult(N,N,N,S,ST);
    
	printf("\nCholensky concluudo:\n");
    imp(N,N,matA);
    printf("\nS=\n");
    imp(N,N,S);
    printf("\nST=\n");
    imp(N,N,ST);
	printf("\nMult SxSt\n");
	imp(N,N,multS);
	
    

    /*
    //Decomposição L e U com temporizacao:
    clock_t Ticks[2];
   
    preencheAleatorio(N,N,matA,100,1);
    //imp(N,N,matA);   
   
    Ticks[0] = clock();
   
    matA = DecLU(N,matA);
   
    //Clock;
    Ticks[1] = clock();
    double Tempo = (Ticks[1] - Ticks[0]) * 1000.0 / CLOCKS_PER_SEC;
    printf("Tempo gasto: %g ms. \n", Tempo);
       
    */
    system("pause");
    return 0;
}