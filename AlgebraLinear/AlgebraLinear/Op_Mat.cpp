#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <time.h>
#include <math.h>

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
void copiaVal(int M, int N, float **matA,float** matB){
    for(int i = 0; i<M;i++)
        for(int j=0;j<N;j++)
            matB[i][j]=matA[i][j];
}
float* copiaVet(int M, float* A){
	float* v = (float *)malloc(sizeof(float)*M);
	for(int i=0;i<M;i++)
		v[i]=A[i];
	return v;
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
  float **matA;
  matA = (float **)malloc(sizeof(float)*c);
  for (int i = 0; i < c; i++){
	matA[i] = (float *)malloc(sizeof(float)*(l));
	for(int j=0;j<l;j++)
		matA[i][j]=A[j][i];
  }
	return matA;
}

void Transposta2 (int l, int c, float** A) {
	float** B = copia(l,c,A);
  for(int i = 0; i < l; i++){
	for(int j=0;j<c;j++)
		A[i][j] = B[j][i];
  }
}
float** VetorColuna(int N, float* V){
	float** v;
	v = (float**)malloc(sizeof(float)*(N));
	for (int i = 0; i < N; i++){
        v[i] = (float *)malloc(sizeof(float));
		v[i][0] = V[i];
		
	}
	return v;
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
float NormaVetor(int N, float* V){
	float norma = 0;
	for(int i=0;i<N;i++)
		norma+=(V[i]*V[i]);
	norma = sqrt(norma);
	return norma;
}
float** sub(int M, int N, float** A, float** B){
	float** C = (float **)malloc(sizeof(float)*M);
    for (int i = 0; i < M; i++){
        C[i] = (float *)malloc(sizeof(float)*(N));
		for(int j=0; j<N;j++)
			C[i][j]=A[i][j]-B[i][j];
	}
	return C;
}
float NormaMatriz(int M, int N, float** V){
	float norma = 0;
	for(int i=0;i<M;i++)
		for(int j=0;j<N;j++)
			norma+=(V[i][j]*V[i][j]);
	return sqrt(norma);
}

float NormaMatriz_inf(int M, int N, float**V){
	float max = 0;
	for(int i=0;i<M;i++)
		for(int j=0;j<N;j++){
			float ab = abs(V[i][j]);
			if(max<ab)
				max = ab;
		}
	return max;		
}
float* subVetor(int N, float* A, float* B){
	float* v = (float *)malloc(sizeof(float)*N);
	for(int i=0;i<N;i++)
		v[i]=A[i]-B[i];
	return v;
}
float* somaVetor(int N, float* A, float* B){
	float* v = (float *)malloc(sizeof(float)*N);
	for(int i=0;i<N;i++)
		v[i]=A[i]+B[i];
	return v;
}
float ProdutoEscalar(int N, float* A, float* B){
	float soma = 0;
	for(int i =0;i<N;i++)
		soma+= A[i]*B[i];
	return soma;
}
void Escalar(int N, float* V, float Esc){
	for(int i=0;i<N;i++)
		V[i] *= Esc;
}
void Escalar_Mat(float E, int M, int N, float**A){
	for(int i=0;i<M;i++)
		for(int j=0;j<N;j++)
			A[i][j]*=E;
}
float** householder(int N, float** n){
	float** E = MatIdentidade(N); //E = I -2n*nt
	float** nt = Transposta(N,1,n);
	float** H = mult(N,1,N,n,nt);
	Escalar_Mat(2,N,N,H);
	E=sub(N,N,E,H);
	return E;
}


float* Gauss (int N, float **A, float *b){

	float** matA = copia(N,N,A);
	float* vet2 = copiaVet(N,b);
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
   ///imp(N,N,matA);
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
float* GaussJordan(int N, float **A, float *b){
	float** matA = copia(N,N,A);
	float* vet2 = copiaVet(N,b);
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
float** DecLU(int N, float **Mat){
    float** A = copia(N,N,Mat);
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

float* SELDecLU(int N, float **matA, float *b){
	float*vet1 = copiaVet(N,b);
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

float* Solve_LU(int N, float **LU, float *vet2){
	
    int i,j;
    float alfa;
    float** L;
    float** U;
	U = copia(N,N,LU);
    L = MatIdentidade(N);
	float *vet1 = (float*)malloc(sizeof(float)*N);
	
    for(i=1;i<N;i++){
		vet1[i-1] = vet2[i-1];
		for(j=0;j<i;j++){
            L[i][j]=U[i][j];
            U[i][j]=0;
        }
	}vet1[N-1]=vet2[N-2];
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
float* SELCholensky(int N, float** A, float *b){
    float*B = copiaVet(N,b);
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

float** MMQ(int N, int M, float* x, float*y ){
	float** matA;
	matA = (float **)malloc(sizeof(float)*(M+1));
    for (int i = 0; i < (M+1); i++)
        matA[i] = (float *)malloc(sizeof(float)*(M+1));
	
	float* v = (float *)malloc(sizeof(float)*((M*2)+1));
	float* v2 = (float *)malloc(sizeof(float)*(M+1));
	v[0]=N;

	matA[0][0]=N;
	
	for(int i=0;i<((M*2)+1);i++){
		float X=0;
		float Y=0;
		if(i==0){
			for(int j = 0;j<N;j++)
				Y += y[j]*pow(x[j],i);
			v2[i]=Y;
		}else{
			for(int j = 0;j<N;j++)
				X+=pow(x[j],i);
			v[i]=X;
			if(i<=M){
				for(int j = 0;j<N;j++)
				Y += y[j]*pow(x[j],i);
			
				v2[i]=Y;
				matA[0][i]= v[i];
				matA[i][0]= v[i];
			}
		}
	}
	for(int i=1;i<=M;i++){
		for(int j=1;j<=M;j++)
			matA[i][j]=v[j+i];
		
	}
	//impVet((2*M)+1,v);
	imp(M+1,M+1,matA);
	impVet(M+1,v2);

	float *solGauss = Gauss(M+1,matA,v2);
	impVet(M+1,v2);
	return matA;
}

float** GaussSeidel_Kint(int N, float** A, float *b, int K){
	
	float** x = (float **)malloc(sizeof(float)*N);
		for (int i = 0; i < N; i++){
			x[i]=(float *)malloc(sizeof(float)*(K));
			x[i][0]=0;
		}

		for(int k=1;k<N;k++){
			
			for(int i=0;i<N;i++){
				float S = 0,S2=0;
				for(int j=0;j<(i-1);j++)
					S+= (A[i][j]*x[j][k-1]);
				
				for(int j=(i+1);j<N;j++)
					S2+= (A[i][j]*x[j][k]);
				
				x[i][k]= (-S-S2+b[i])/A[i][i];
			}
		}
		return x;

}

float* Jacobi(int N, float** A, float *b, float e, int kMax){
	float*x,*x_novo;

	x=(float *)malloc(sizeof(float)*N);
	x_novo=(float *)malloc(sizeof(float)*N);

	bool tool=false;
	float n1=0,n2=0;
	int k=0;
	while(k < kMax && !tool)
	{
		for(int i=0;i<N;i++){
		x_novo[i]=0;
		}

		for(int i=0; i<N;i++){
			float soma = 0;
			for(int j=0;j<N;j++){
				//printf("\n soma = %f", soma);
				if(i!=j)
					soma += (A[i][j]*x[j]);
			}
			x_novo[i] = (b[i]-soma)/A[i][i];
		}

		n1 = NormaVetor(N,x);
		n2 = NormaVetor(N,x_novo);

		if(abs(n2-n1)<e)
			tool=true;
		else
			for(int i=0;i<N;i++){
				x[i] = x_novo[i];
			}
		k++;
	}
	printf("\n Num Interacoes = %d",k);
	return x_novo;

}
float* Gauss_Seidel(int N, float** A, float *b, float e, int kMax){
	float*x,*x_novo;

	x=(float *)malloc(sizeof(float)*N);
	x_novo=(float *)malloc(sizeof(float)*N);
	
	for(int i=0;i<N;i++){
		x[i]=0;
	}

	bool tool=false;
	float n1=0,n2=0;
	int k=0;
	while(k < kMax && !tool)
	{

		for(int i=0; i<N;i++){
			float soma = 0;
			for(int j=0;j<i;j++)
				soma += (A[i][j]*x[j]);
			for(int j=i+1;j<N;j++)
				soma += (A[i][j]*x_novo[j]);

			x_novo[i] = (b[i]-soma)/A[i][i];
		}

		n1 = NormaVetor(N,x);
		n2 = NormaVetor(N,x_novo);

		if(abs(n2-n1)<e)
			tool=true;
		else
			for(int i=0;i<N;i++){
				x[i] = x_novo[i];
			}
		k++;
	}
	
	printf("\n Num Interacoes = %d",k);
	return x_novo;

}

float* SOR(int N, float** A, float *b, float e, int kMax, float W){
	float*x,*x_novo;

	x=(float *)malloc(sizeof(float)*N);
	x_novo=(float *)malloc(sizeof(float)*N);
	
	for(int i=0;i<N;i++){
		x[i]=0;
	}

	bool tool=false;
	float n1=0,n2=0;
	int k=0;

	while(k < kMax && !tool)
	{

		for(int i=0; i<N;i++){
			float soma1 = 0;
			float soma2 =0;
			for(int j=0;j<i;j++)
				soma1 += (A[i][j]*x[j]);
			for(int j=i+1;j<N;j++)
				soma2 += (A[i][j]*x_novo[j]);

			float GS = (b[i]-soma1-soma2);
			x_novo[i] = ((1-W)*x[i]) + ((W/A[i][i])*GS);
		}

		n1 = NormaVetor(N,x);
		n2 = NormaVetor(N,x_novo);

		if(abs(n2-n1)<e)
			tool=true;
		else
			for(int i=0;i<N;i++){
				x[i] = x_novo[i];
			}
		k++;
	}
	
	printf("\n Num Interacoes = %d",k);
	return x_novo;

}


float** QR_givens(int M, int N, float** Mat, float **R){
	float** q = MatIdentidade(M);
	float **r = copia(M,N,Mat);
	for(int j=0;j<N;j++){
		for(int i=j+1;i<M;i++){
			float raio,c,s;
			raio = sqrt(r[i][j]*r[i][j] + r[j][j]*r[j][j]);
			s = r[i][j]/raio;
			c = -r[j][j]/raio;
			float** J = MatIdentidade(M);
			J[j][j]=c; J[i][i]=c;
			J[i][j]=-s;J[j][i]=s;
			J = Transposta(M,M,J);
			r = mult(M,M,N,J,r);
			q = mult(M,M,M,J,q);
		}

	}
	q = Transposta(M,M,q);
	copiaVal(M,N,r,R);
	return q;
}
float** QR_HouseHolder(int M, int N, float** Mat, float **R){

	float** Q = MatIdentidade(M);
	float **r = copia(M,N,Mat);
	float** H = MatIdentidade(M);
	float*v = (float*)malloc(sizeof(float)*N);
	float norm=0;
	float*n;

	for(int j=0;j<N-1;j++){

		for(int i=0;i<N;i++)
			if(i<j)
				v[i]=0;
			else{
				v[i]=r[i][j];
			}
		norm = NormaVetor(N,v);
		//printf("V: \n");
		//impVet(N,v);
		n = copiaVet(N,v);
		
		n[j] -= norm;
		norm = NormaVetor(N,n);
		
		Escalar(N,n,1/norm);
		
		//printf("N: \n");
		//impVet(N,n);
		H = householder(N,VetorColuna(N,n));

		//imp(M,N,H);

		r = mult(N,N,N,H,r);
		//imp(M,N,r);

		Q = mult(M,N,N,Q,H);
		
	}
	//Q = Transposta(M,M,Q);
	copiaVal(N,N,r,R);
	return Q;
}


float** Gram_Schmidt(int M, int N, float** A){
	float** Q = Transposta(M,N,A);
	Escalar(M,Q[0],1/NormaVetor(M,Q[0]));
	
	for(int i=1;i<N;i++){
		float *w = copiaVet(M,Q[i]);
		w = subVetor(M,w,w);
		for(int j=0;j<i;j++){
			float* v =  copiaVet(M,Q[j]);
			float x = ProdutoEscalar(M,v,Q[i]);
			Escalar(M,v,x);
			w = somaVetor(M,v,w);
		}
		Q[i]=subVetor(M,Q[i],w);
		Escalar(M,Q[i],1/NormaVetor(M,Q[i]));

	}
	//Q = Transposta(N,M,Q);
	return Q;
}
float* SEL_QR(int M, int N, float** A, float* b){
	float** Q = Gram_Schmidt(M,N,A); // Q tranposta 
	//imp(N,M,Q);
	float** r = mult(N,M,N,Q,A);
	float** B = mult(N,M,1,Q,VetorColuna(M,b));
	float* v = (float *)malloc(sizeof(float)*N);
	//imp(N,N,r);
	//imp(N,1,B);


	for(int i=(N-1);i>0;i--){
		float soma = 0;
		for(int k = i+1; k < N; k++)
			soma += r[i][k]*B[k][0];
		B[i][0] = (B[i][0]-soma)/r[i][i]; 
	}
	//imp(N,1,B);

	for(int i=0;i<N;i++)
		v[i]=B[i][0];

	//float** T = mult(M,N,N,Transposta(N,M,Q),r);
	//imp(M,N,T);

	return v;
}

void Maior_Auto_Valor_Vetor(int Int, int M, float** A, float e, float*v, float &l){
	float **q,**qt,**y,**Lamb;

	float l1=0,l2=0;
	qt = (float **)malloc(sizeof(float));
	
	qt[0] = copiaVet(M,v);
	Escalar(M,qt[0],1/NormaVetor(M,qt[0]));
	q=Transposta(1,M,qt);
	y = mult(M,M,1,A,q);
	Lamb = mult(1,M,1,qt,y);

	int k=0;
	float error=5;
	while(k!= Int && error>e){
	
		l1 = Lamb[0][0];
		qt= Transposta(M,1,y);
		Escalar(M,qt[0],1/NormaVetor(M,qt[0]));
		q = Transposta(1,M,qt);
		y = mult(M,M,1,A,q);
		Lamb = mult(1,M,1,qt,y);
		l2 = Lamb[0][0];
		error = abs((l2-l1)/l2);
		k++;

		/*
		printf("Q");
		imp(M,1,q);
		printf("Y");
		imp(M,1,y);
		printf("Qt");
		imp(1,M,qt);

		printf("Valor k =%d Teste autoValor = %f, error = %f\n",k,l2,error);
		*/
	}
	l=l2;
	//printf("valor k =%d Teste autoValor = %f",k,l2);
	for(int i=0;i<M;i++)
		v[i]=y[i][0];
	
}
void PR(int Int, int N, float** A, float e, float* v, float &l){
	float **qt, **Y, **lamb;

	bool Tol = false;
	int k =0;
	float L_anterior =0;
	float L_corrente=0;
	qt = (float **)malloc(sizeof(float));
	qt[0] = copiaVet(N,v);
	
	Escalar(N,qt[0],1/NormaVetor(N,qt[0]));
	Y = mult(N,N,1,A,VetorColuna(N,qt[0]));
	lamb = mult(1,N,1,qt,Y); 

	L_corrente=lamb[0][0];


	while(k!=Int && !Tol){
		L_anterior=L_corrente;
		qt = Transposta(N,1,Y);
		//imp(1,N,qt);
		Escalar(N,qt[0],1/NormaVetor(N,qt[0]));
		//imp(1,N,qt);
		Y = mult(N,N,1,A,VetorColuna(N,qt[0]));
		//imp(N,1,Y);
		lamb = mult(1,N,1,qt,Y); 
		//imp(1,1,lamb);
		L_corrente=lamb[0][0];
		

		if(abs((L_anterior - L_corrente)/L_anterior) < e)
			Tol=true;

		k++;
	}
	l=L_corrente;
	printf("\n %d iteracoes ",k); 
	for(int i=0;i<N;i++)
		v[i]=Y[i][0];
}

void invPR(int Int, int N, float** A, float e, float* v, float &l){
 	float** LU = DecLU(N,A);

	float **qt, *Y;
	float norm=0;
	bool Tol = false;
	int k =0;
	float L_anterior =0;
	float L_corrente=0;
	
	qt = (float **)malloc(sizeof(float)*N);

    qt[0]=copiaVet(N,v);
    Escalar(N,qt[0],1/NormaVetor(N,v));

	Y = Solve_LU(N,LU,qt[0]); 

	L_corrente=ProdutoEscalar(N,Y,qt[0]);
	

	while(k!=Int && !Tol){
		

		norm =NormaVetor(N,Y);
		for(int i=0;i<N;i++)
			qt[0][i]= (Y[i]/norm);

		Y = Solve_LU(N,LU,qt[0]);
		
		L_corrente=ProdutoEscalar(N,Y,qt[0]);

		if(L_corrente==L_anterior)
			Tol=true;
		k++;
		L_anterior=L_corrente;
	}

	l=1/L_corrente;

	norm =NormaVetor(N,Y);
	printf("\n %d iteracoes ",k); 
	for(int i=0;i<N;i++)
		v[i]= Y[i]/norm;


}

void Menor_Auto_Valor_Vetor(int Int, int M, float** A, float e, float*v, float &l){

	float **qt,**y,**Lamb;
	float* temp;
	float l1=0,l2=0;
	float** LU = DecLU(M,A);

	qt = (float **)malloc(sizeof(float));
	
	qt[0] = copiaVet(M,v);
	float normaQ = NormaVetor(M,qt[0]);
	Escalar(M,qt[0],(1/normaQ));
	
	temp= Solve_LU(M,LU,qt[0]); //Gauss(M,A,qt[0]);
	y = VetorColuna(M,temp);
	Lamb = mult(1,M,1,qt,y);

	int k=0;
	float error=1;
	while(k!=Int && error>e){
	
		l1 = Lamb[0][0];
		qt = Transposta(M,1,y);
		normaQ = NormaVetor(M,qt[0]);
		Escalar(M,qt[0],(1/normaQ));
		temp = Solve_LU(M,LU,qt[0]);//•Gauss(M,A,qt[0]);
		y = VetorColuna(M,temp);
		Lamb = mult(1,M,1,qt,y);
		l2 = Lamb[0][0];
		error = abs((l2-l1)/l2);
		k++;
		//printf("Valor k =%d Teste autoValor = %f, error = %f\n",k,1/l2,error);
		
	}
	l=1/l2;
	
	for(int i=0;i<M;i++)
		v[i]=y[i][0];
	
}
void Desloc_Auto_Valor_Vetor(int Int, int M, float** A, float e, float*v, float &l, float x){
	float** MI = MatIdentidade(M);
	for(int i=0;i<M;i++){
		MI[i][i] = x;
	}
	float Lb =0;
	float** Ab = sub(M,M,A,MI);
	//imp(M,M,Ab);
	invPR(Int,M,Ab,e,v,Lb);
	l=Lb+x;
}

float** S_householder(int N, float **A, float** Vetores){
	float** H_fim = MatIdentidade(N);
	float **H, *v, *n, norm;
	float **AT = copia(N,N,A);
	v = (float*)malloc(sizeof(float)*N);
	norm=0;
	//imp(N,N,AT);
	for(int j=0; j<(N-2);j++){
		for(int i=0;i<N;i++)
			if(i<=j){
				v[i]=0;
			}else{
				v[i]= A[i][j];
				
			}
	   // impVet(N,v);
		norm = NormaVetor(N,v);
		//printf("\nTeste Norm %f",norm);
		n = copiaVet(N,v);
		n[j+1]-=norm;
		norm = NormaVetor(N,n);
		if(norm!=0)
			Escalar(N,n,1/norm);
		//impVet(4,v);
		H = householder(N,VetorColuna(N,n));
		//imp(N,N,H);
		//imp(N,N,AT);
		AT = mult(N,N,N,AT,H);
		AT = mult(N,N,N,H,AT);
		
		//imp(N,N,AT);
		H_fim = mult(N,N,N,H_fim,H);
		
		//imp(N,N,H_fim);
	}
	//printf("\n AT e H(H1*H2*H3...*Hn-2) \n");
	//imp(N,N,AT);
	//imp(N,N,H_fim);
	copiaVal(N,N,H_fim,Vetores);
	return AT;
}

float** Diag_QR(int N, int Int, float** Mat, float** Auto_Vet){
	int k=0;
	float** r = MatIdentidade(N);
	float** A = copia(N,N,Mat);
	float** Q = MatIdentidade(N);
	float** q;
	float e = 0.000001;
	bool Inter = true;
	float n2 = 0;
	for(int i=0;i<N;i++)
		n2+=NormaVetor(N,A[i]);
	int cont = 0;
	while(k!=Int && Inter){
		float n1 = n2;
		q = QR_givens(N,N,A,r);	
		Q = mult(N,N,N,Q,q);
		A = mult(N,N,N,r,q);
		n2 = 0;
		
		for(int i=0;i<N;i++)
			n2+=NormaVetor(N,A[i]);

		//printf("\nNorma de A[k-1] = %f\n", n1);
		//printf("\nNorma de A[k] = %f\n", n2);
		if(abs((n2-n1))<e){
			cont++;
			if(cont==29)
				Inter = false;
		}else{
			cont=0;
		}
		k++;

	}
	//printf("\n Valor K %d; \n", k);
	printf("\n %d iteracoes para Diagonalizacao QR ",k); 
	for(int i=0;i<(N-1);i++){
		if(abs(A[i+1][i])>e)
			printf("\n Auto Valores complexos \n");
	}

	copiaVal(N,N,Q,Auto_Vet);
	return A;
}

void AutoValores_Vetores(int N, float** Mat, float** Vetores, float* valores){
	float** H = MatIdentidade(N);
	float** Q = MatIdentidade(N);

	float** T = S_householder(N,Mat,H);
	if(T[N-1][0]==0){
		float** D = Diag_QR(N,100000,T,Q);;
		float** Vet = mult(N,N,N,H,Q);
		copiaVal(N,N,Vet,Vetores);
		for(int i=0;i<N;i++)
			valores[i]=D[i][i]; 
	}else{
		printf("\n\nErro no processo de Tridiagonalizacao");
		float** D = Diag_QR(N,100000,Mat,Q);
		copiaVal(N,N,Q,Vetores);
		for(int i=0;i<N;i++)
			valores[i]=D[i][i];
	
	}
	
}

void SVD(int M, int N, float** Mat, float** U, float** S, float** V){

 
	float** AT = Transposta(M,N,Mat);
	//imp(M,N,AT);
	float** AtA = mult(N,M,N,AT,Mat);
	float* autoVal = (float*)malloc(sizeof(float)*N);
	AutoValores_Vetores(N,AtA,V,autoVal);
	float** u = mult(M,N,N,Mat,V);
	Transposta2(N,N,V);
	for(int i=0; i<N; i++){
		S[i][i]= sqrt(autoVal[i]);
	}
	for(int i=0; i<M; i++){
		for(int j=0;j<N;j++){
			U[i][j] = u[i][j]/S[j][j];
		}
	}
	
}

float* Conj_Gradiente(int N, float **A, float*b, float*x0, float e){
	float **g, **gt, **x, **S, **St, **temp, Lamb=0, Beta=0;
	
	x = (float**)malloc(sizeof(float)*N);
	g = (float**)malloc(sizeof(float)*N);

	for(int i=0; i<N;i++){
	x[i] = (float *)malloc(sizeof(float));
		x[i][0] = x0[i];
	}

	temp = mult(N,N,1,A,x);
	

	for(int i=0; i<N;i++){
	g[i] = (float *)malloc(sizeof(float));
		g[i][0] = b[i]-temp[i][0];
	}

	S = copia(N,1,g);

	int k = 0;
	float norm = NormaMatriz(N,1,g);

	while(k<5000 && norm>e){
		St = Transposta(N,1,S);
		
		float** AS = mult(N,N,1,A,S);

		temp = mult(1,N,1,St,g);
		float St_g = temp[0][0];
		

		temp = mult(1,N,1,St,AS);
		float St_A_S = temp[0][0];

		

		Lamb = St_g/St_A_S;
	

		for(int i=0;i<N;i++)
			x[i][0] = x[i][0] + (Lamb * S[i][0]);

		//printf("Lanb = %f",Lamb);
	

		temp = mult(N,N,1,A,x);

		for(int i=0;i<N;i++)
			g[i][0] = b[i]-temp[i][0];

		gt = Transposta(N,1,g);
		temp = mult(1,N,1,gt,AS);
		float gT_A_S = temp[0][0];

		Beta = -(gT_A_S/St_A_S);

		for(int i=0;i<N;i++)
			S[i][0]=g[i][0]+ (Beta*S[i][0]);

		norm = NormaMatriz(N,1,g);
		k++;
	}

	printf("\n %d iteracoes ",k); 
	float *Ret = (float*)malloc(sizeof(float)*N);

	for(int i=0; i<N;i++)
		Ret[i] = x[i][0];

	return Ret;


}

int main()
{
int M =4;
int N =4;

// float matA[6][6] = {{20,10,0,0,0,0},{10,20,10,0,0,0},{0,10,20,10,0,0},{0,0,10,20,10,0},{0,0,0,10,20,10},{0,0,0,0,10,20}};
 //float B[6] = {-10,10,20,20,10,-10};	

 float matA[4][4] = {{3,2,1,-1},{0,-3,-5,7},{0,1,0,3},{0,2,4,0}};  // Não é positiva definida
float B[4] = {5,7,6,15};

 
 //float matA[3][3] = {{12.0,-51.0,4}, {6.0,167.0,-68.0}, {-4.0,24.0,-41.0}}; // Não é Simetrica
//float matA[4][4] = {{7, 3, -1, 2},{3, 8, 1, -4},{-1, 1, 4, -1},{2, -4, -1, 6}};

//float matA[5][3] = {{4,3,2},{3,4,3},{2,3,4},{1,2,3},{0,1,2}};




 float *b = (float*)malloc(sizeof(float)*M);
 float *x = (float*)malloc(sizeof(float)*M);


 float **A = (float **)malloc(sizeof(float)*M);
    for (int i = 0; i < M; i++){
        A[i] = (float *)malloc(sizeof(float)*(N));
		A[i] = matA[i];
		b[i] = B[i];
	}
	

	int K = 99; 


	while(K!=0){


		printf("\n Algebra Linear, digite: \nDigite 1 - Resolucao de Sistemas Lineares  \nDigite 2 - Decomposicoes\nDigite 3 - Autovalores e Autovetores \nDigite 0 - Sair\n\n\n");
		scanf("%d", &K);

		if(K==1){
			float *x0;
			

			int SEL = 0;
			while(SEL!=9){
			printf("\nDigite : \n1-Gauss\n2-Gauss-Jordan\n3-LU\n4-Cholenky\n5-Jacobi\n6-Gauss_Seidel\n7-SOR\n8-Gradiente\n9-Voltar\n");

			scanf("%d", &SEL);

			float e=0; 
		     int kMax =0;
			float w =0;

			printf("\n Para o Sistema A x = b");



			switch (SEL)
			{
			   case 1:
				 printf("\nEliminacao de Gauss : ");
				 printf("\n Matriz A: ");
				 imp(N,N,A);
				 printf("\n Vetor b: \n");
				 impVet(N,b);
				 x = Gauss(N,A,b);
				 printf("\nVetor Solucao:\n");
				 impVet(N,x);

			   break;
					
			   case 2:
				 printf("\nGauss-Jordan: ");
				 x = GaussJordan(N,A,b);
				 printf("\n Matriz A: ");
				 imp(N,N,A);
				 printf("\n Vetor b: \n");
				 impVet(N,b);
				 printf("\nVetor Solucao:\n");
				 impVet(N,x);
				break;
			   case 3: 
				    printf("\nLU : ");
				 printf("\n Matriz A: ");
				 imp(N,N,A);
				 printf("\n Vetor b: \n");
				 impVet(N,b);
				 x = SELDecLU(N,A,b);
				 printf("\nVetor Solucao:\n");
				 impVet(N,x);
				 break;

				 case 4: 
				   printf("\nCholensky : ");
				 printf("\n Matriz A: ");
				 imp(N,N,A);
				 printf("\n Vetor b: \n");
				 impVet(N,b);
				 x = SELCholensky(N,A,b);
				 printf("\nVetor Solucao:\n");
				 impVet(N,x);
				 break;

			   case 5: 
				 e=0;kMax=0;	
				 printf("\n Digite: ");
				 printf("\n Tolerancia: ");
				 scanf("%f", &e);
				 printf("\n Max Interacoes: ");
				 scanf("%d", &kMax);

				 printf("\n Jacobi: ");
				 printf("\n Matriz A: ");
				 imp(N,N,A);
				 printf("\n Vetor b: \n");
				 impVet(N,b);
				 x = Jacobi(N,A,b,e,kMax);
				 printf("\nVetor Solucao:\n");
				 impVet(N,x);

			   break;

			   case 6: 
				 e=0;kMax=0;	
				 printf("\n Digite: ");
				 printf("\n Tolerancia: ");
				 scanf("%f", &e);
				 printf("\n Max Interacoes: ");
				 scanf("%d", &kMax);

				 printf("\n Gauss_seidel: ");
				 printf("\n Matriz A: ");
				 imp(N,N,A);
				 printf("\n Vetor b: \n");
				 impVet(N,b);
				 x = Gauss_Seidel(N,A,b,e,kMax);
				 printf("\nVetor Solucao:\n");
				 impVet(N,x);

			   break;

			   case 7: 
				 e=0;kMax=0;	
				 printf("\n Digite: ");
				 printf("\n Tolerancia: ");
				 scanf("%f", &e);
				 printf("\n Max Interacoes: ");
				 scanf("%d", &kMax);
				 printf("\n W: ");
				 scanf("%f", &w);

				 printf("\n SOR: ");
				 printf("\n Matriz A: ");
				 imp(N,N,A);
				 printf("\n Vetor b: \n");
				 impVet(N,b);
				 x = SOR(N,A,b,e,kMax,w);
				 printf("\nVetor Solucao:\n");
				 impVet(N,x);

			   break;
			     case 8: 
				 x0 = (float*)malloc(sizeof(float)*N);
				 for(int i=0;i<N;i++)
					x0[i]=1;
				 printf("\n Tolerancia: ");
				 scanf("%f", &e);
				 printf("\n Gradiente: ");
				 printf("\n Matriz A: ");
				 imp(N,N,A);
				 printf("\n Vetor b: \n");
				 impVet(N,b);
				 x = Conj_Gradiente(N,A,b,x0,e);
				 printf("\nVetor Solucao:\n");
				 impVet(N,x);

			   break;

			}
			}
		}
		
		if(K == 2){

			int Dec = 0;
			float** L, **U, **S, **St, **Q, **R,**Mat, **Temp,**V;
			while(Dec!=7){
			printf("\nDigite : \n1-LU\n2-Cholensky\n3-QR_Gram Schmidt\n4-QR_HouseHolder\n5-QR_Givens\n6-SVD\n 7- Voltar\n");
			scanf("%d", &Dec);

				switch (Dec)
				{
				   case 1:
					   if(M==N){
						   printf("\n Matriz A: ");
						   imp(M,N,A);
						   L = MatIdentidade(N); 
						   U = DecLU(N,A);
							for(int i=1;i<N;i++)
								for(int j=0;j<i;j++){
									L[i][j]=U[i][j];
									U[i][j]=0;
								}

							printf("\nL:");
							imp(M,N,L);
							printf("\nU:");
							imp(M,N,U);

							Temp = mult(M,N,M,L,U);
							printf("\nTeste L*U:");
							imp(M,N,Temp);
					   }else
						   printf("\n M!=N ");
					   

				   break;
				   case 2: 
					   if(M==N){
						   printf("\n Matriz A: ");
						   imp(M,N,A);
						   S = Cholesky(N,A);
						   St = Transposta(N,N,S);

						   printf("\nS:");
						   imp(N,N,S);
						   printf("\nSt:");
						   imp(N,N,St);

						    Temp = mult(M,N,M,S,St);
							printf("\nTeste S*St:");
							imp(M,N,Temp);
					   }else
						   printf("\n M!=N ");
					break;

					case 3: 
	                   Mat = copia(M,N,A);
					   printf("\n Matriz A: ");
					   imp(M,N,Mat);
					   Q = Gram_Schmidt(M,N,Mat); 
					   R = mult(N,M,N,Q,Mat);
					   Q = Transposta(N,M,Q);
					   printf("\n Matriz Q: ");
					   imp(M,N,Q);
					   printf("\n Matriz R: ");
					   imp(N,N,R);
					   Temp = mult(M,N,M,Q,R);
					    printf("\nTeste Q*R:");
						imp(M,N,Temp);

					break;
					case 4: 
					   printf("\n Matriz A: ");
					   imp(M,N,A);
					   R = MatIdentidade(M);
					   Q = QR_HouseHolder(M,N,A,R);
					   printf("\n Matriz Q: ");
					   imp(M,N,Q);
					   printf("\n Matriz R: ");
					   imp(N,N,R);
					   Temp = mult(M,N,N,Q,R);

					   printf("\nTeste Q*R:");
					   imp(M,N,Temp);
					break;
					case 5: 
					   printf("\n Matriz A: ");
					   imp(M,N,A);
					   R=MatIdentidade(M);
					   Q = QR_givens(M,N,A,R);
					   printf("\n Matriz Q: ");
					   imp(M,N,Q);
					   printf("\n Matriz R: ");
					   imp(N,N,R);
					   Temp = mult(M,N,N,Q,R);
					    printf("\nTeste Q*R:");
						imp(M,N,Temp);
					break;

					case 6:
						
						U = MatIdentidade(M);
						S = MatIdentidade(M);
						V = MatIdentidade(M);

						SVD(M,N,A,U,S,V);

						printf("\n Matriz Original");
						imp(M,N,A);
						printf("\nU:");
						imp(M,N,U);
						printf("\nS:");
						imp(N,N,S);
						printf("\nVt:");
						imp(N,N,V);

						printf("\n Teste U*S*Vt :");
						Temp = mult(M,N,N,U,S);
						Temp = mult(M,N,N,Temp,V);
						imp(M,N,Temp);

				}
			
			}
		}
		if(K==3){
		
			float** Vetores, **T, **D, *Valores;
			int x = 0;
			float* v;
			float l=0;
			float error=0;
			int Max =0;
			float Desl =0;
			
			while(x!=7){
				printf("\nDigite : \n1-Potencia Regular \n2-Potencia Inversa \n3-Potencia com Deslocamento \n4-Tridiagonalizacao\n5-Diagonalizacao QR \n6-Processos 4 e 5 para Todos os Autovetores e autovalores \n 7- Voltar\n");
				scanf("%d", &x);
				switch (x)
				{
					case 1:
					    error=0;Max=0;l = 0;
						printf("\nDigite a Tolerancia\n");
						scanf("%f", &error);
						printf("\nDigite o numero maximo de iteracoes\n");
						scanf("%d", &Max);
						v = (float*)malloc(sizeof(float)*N);
						for(int i=0;i<N;i++)
							v[i]=1;
						
						PR(Max,N, A, error, v,l);
						printf("\nMatriz A\n");
						imp(N,N,A);
						printf("\nAuto Vetor\n");
						impVet(N,v);
						printf("\nAutoValor %f\n",l);
						//free(v);
					break;
					case 2:
						error=0;Max=0;l = 0;
						printf("\nDigite a Tolerancia\n");
						scanf("%f", &error);
						printf("\nDigite o numero maximo de interacoes\n");
						scanf("%d", &Max);
						v = (float*)malloc(sizeof(float)*N);
						for(int i=0;i<N;i++)
							v[i]=1;
						l = 0;	
						invPR(Max,N, A, error, v,l);
						printf("\nMatriz A\n");
						imp(N,N,A);
						printf("\nAuto Vetor\n");
						impVet(N,v);
						printf("\nAutoValor %f\n",l);

					break;
					case 3:
					    error=0;Max=0;l = 0;Desl=0;
						printf("\nDigite a Tolerancia\n");
						scanf("%f", &error);
						printf("\nDigite o numero maximo de iteracoes\n");
						scanf("%d", &Max);
						printf("\nDigite o Deslocamento\n");
						scanf("%f", &Desl);

						v = (float*)malloc(sizeof(float)*N);
						for(int i=0;i<N;i++)
							v[i]=1;

						Desloc_Auto_Valor_Vetor(Max,N, A, error, v,l, Desl);
						printf("\nMatriz A\n");
						imp(N,N,A);
						printf("\nAuto Vetor\n");
						impVet(N,v);
						printf("\nAutoValor %f\n",l);

					break;
					case 4:
						Vetores = MatIdentidade(N);
						T = S_householder(N,A,Vetores);
						printf("\n Matriz Original \n");
						imp(N,N,A);
						printf("\nT");
						imp(N,N,T);
						printf("\nAutovetores de T");
						imp(N,N,Vetores);
						//free(Vetores);
						//free(T);

					break;
					case 5:
					    Max=0;
					    Vetores = MatIdentidade(N);
						printf("\nDigite o numero maximo de iteracoes\n");
						scanf("%d", &Max);
						Vetores = MatIdentidade(N);
						D = Diag_QR(N,Max,A,Vetores);
						printf("\n Matriz Original \n");
						imp(N,N,A);
						printf("\nD");
						imp(N,N,D);
						printf("\nAutovetores de D");
						imp(N,N,Vetores);
						//free(Vetores);
						//ree(D);

					break;
					case 6:

					Vetores = MatIdentidade(N);
					Valores = (float*)malloc(sizeof(float)*N);
					AutoValores_Vetores(N,A,Vetores,Valores);

					printf("\n Matriz Original \n");
					imp(N,N,A);
					printf("\nAutovetores");
					imp(N,N,Vetores);
					printf("\nAutovalores: \n");
					impVet(N,Valores);
					//free(Vetores);
					//free(Valores);

					break;
				}

			}
		}

	}


//


	/*
	//float B[4][4]= {{18.324146,0.001371,0.000357,0.000001},{0.001370,8.959310,-0.831104,-0.000003},{0.000358,-0.831104,-7.037904,-0.000019},{0.000000,-0.000003,-0.000019,-2.245561}};
	float A[4][4] = {{4,0,0,0},{0,2,4,0},{0,0,2,0},{0,0,0,2}};
	//float A[6][6] = {{20,10,0,0,0,0},{10,20,10,0,0,0},{0,10,20,10,0,0},{0,0,10,20,10,0},{0,0,0,10,20,10},{0,0,0,0,10,20}};
	/*20 10  0  0 0   0
    10 20 10  0 0   0
     0 10 20 10 0   0
     0 0  10 20 10  0
     0 0   0 10 20 10
     0 0   0  0 10 20*/
	
	//float A[3][5] = {{4,2,2,1,0},{3,1,2,2,1},{2,3,2,3,2}};
/*
	float A[4][4] = {{4,0,0,0},{0,2,4,0},{0,0,2,0},{0,0,0,2}};
    int M = 4;
	int N = 4;

    float **matA;

    matA = (float **)malloc(sizeof(float)*M);
    for (int i = 0; i < M; i++){
        matA[i] = (float *)malloc(sizeof(float)*(N));
		matA[i] = A[i];
	}
	//float* b = (float*)malloc(sizeof(float)*N);
	//preencheVet(N,b);

	float* v = (float*)malloc(sizeof(float)*M);
	for(int i=0;i<M;i++)
		v[i]=1;
	float l=0;
	
	float e = 0.000001;


	Maior_Auto_Valor_Vetor(1000,M, matA, e, v,l);
	//Menor_Auto_Valor_Vetor(10000,M,matA,e,v,l);
	//Desloc_Auto_Valor_Vetor(10000,M,matA,e,v,l,6);

	printf("\nAuto Vetor\n");
	impVet(M,v);
	printf("\nAutoValor %f\n",l);


	/*

	float *x = SOR(N, matA, b, 0.00001, 15000,0.5);
	impVet(N,x);

	imp(M,N,matA);

	float** Result = mult(M,N,1,matA,VetorColuna(N,x));
	imp(N,1,Result);

	/*

	float**R = (float **)malloc(sizeof(float)*M);
    for (int i = 0; i < M; i++){
        R[i] = (float *)malloc(sizeof(float)*(N));
	}
	float** Q = QR_HouseHolder(M, N,matA, R); //QR_givens(M,N,matA,R);//

	imp(M,N,matA);

	imp(M,N,Q);
	imp(M,N,R);
	float** Mat = mult(M,N,M, Q, R);
	imp(M,N,Mat);

	//preenche(M,N,matA);
	/*
	if(N>M){
		matA = Transposta(M,N,matA);
		float x = M;
		M=N;
		N=x;
	}*/

	//imp(4,4,matA);
	
	
	/*
	float* b = (float*)malloc(sizeof(float)*M);
	preencheVet(4,b);
	//impVet(4,b);
	float *x = Gauss(4, matA,b);
	//impVet(4,x);

	/*
	float** U = MatIdentidade(M);
	float** S = MatIdentidade(N);
	float** V = MatIdentidade(N);

	SVD(M,N,matA,U,S,V);

	printf("\n Matriz Original \n");
	imp(M,N,matA);
	printf("\nU:\n");
	imp(M,N,U);
	printf("\nS:\n");
	imp(N,N,S);
	printf("\nV:\n");
	imp(N,N,V);
	
	printf("\nU*S*Vt:\n");
	float** Temp = mult(M,N,N,U,S);
	Temp = mult(M,N,N,Temp,V);
	imp(M,N,Temp);


	/*
	float** Vetores = MatIdentidade(N);
	float* Valores = (float*)malloc(sizeof(float)*M);
	AutoValores_Vetores(N,matA,Vetores,Valores);
	printf("\n Matriz Original \n");
	imp(M,N,matA);
	printf("\nAuto Vetores");
	imp(N,N,Vetores);
	printf("\nAuto valores: \n");
	impVet(N,Valores);
	*/

	/*
	float* vetor = (float*)malloc(sizeof(float)*M);
	float x = 0;
	Maior_Auto_Valor_Vetor(500,M,matA,0.000001, vetor, x);
	printf("\n Auto vetor (Potencia regular) : \n");
	impVet(M,vetor);
	printf("\n Auto valor associado : %f ", x);
	*/
	//float** Teste = mult(4,4,1,matA,Vetores);
	//imp(N,1,Teste);

	/*
	float** T = S_householder(M,matA,Vetores);
	printf("Teste tridiagonal\n");
	
	
	imp(N,N,T);
	float** D = Diag_QR(M,500,T, Vetores);

	printf("Teste Diagonalizacao\n");
	imp(N,N,D);
	imp(N,N,Vetores);
	*/
	//printf("Teste Tri-diagonal\n");
	
	//float** VetT = MatIdentidade(N);
	//float** T = S_householder(N,matA, VetT);
	//imp(N,N,T);
	//imp(N,N,VetT);
	//printf("Teste Diagonalizacao\n");


	//float** D2 = Diag_QR(N,500,matA,Vetores);
	//imp(N,N,D2);
	//imp(N,N,Vetores);
	

	//S_householder(M, matA);
	
	

	/*
	float* v = (float*)malloc(sizeof(float)*M);
	for(int i=0;i<M;i++)
		v[i]=1;
	float l=0;
	
	float e = 0.000001;


	//Maior_Auto_Valor_Vetor(1000,M, matA, e, v,l);
	//Menor_Auto_Valor_Vetor(10000,M,matA,e,v,l);
	Desloc_Auto_Valor_Vetor(10000,M,matA,e,v,l,6);

	printf("\nAuto Vetor\n");
	impVet(M,v);
	printf("\nAutoValor %f\n",l);

	float** Valida = mult(M,M,1,matA,VetorColuna(M,v));;

	printf("\n Valida: \n ");
	imp(M,1,Valida);
	*/

    return 0;
}
