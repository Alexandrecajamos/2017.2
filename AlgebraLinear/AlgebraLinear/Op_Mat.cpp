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

float** Jacobi_Kint(int N, float** A, float *b, int K){
	
	float** x = (float **)malloc(sizeof(float)*N);
		for (int i = 0; i < N; i++){
			x[i]=(float *)malloc(sizeof(float)*(K));
			x[i][0]=0;
		}

		for(int k=1;k<N;k++){
			
			for(int i=0;i<N;i++){
				float som = 0;
				for(int j=0;j<N;j++)
					som+= (A[i][j]*x[j][k-1]);
			
				x[i][k]= (-som+b[i])/A[i][i];
			}
		}
		return x;

}
float* Jacobi_error(int N, float** A, float *b, float e){

	float *Xk = (float *)malloc(sizeof(float)*(N));
	float *Xkb = (float *)malloc(sizeof(float)*(N));
	for(int i=0;i<N;i++)
		Xk[i]=0;
	bool Error = true;
	while(Error){
		for(int i=0;i<N;i++){
			float som = 0;
			for(int j=0;j<N;j++)
				som+= (A[i][j]*Xk[j]);
			Xkb[i]= (-som+b[i])/A[i][i];
		}
		int cont=0;
		float Er = NormaVetor(N,Xkb);

		if(Er<e)
			Error=false;	
		//Se o Erro ainda não foi satisfeito, enstão copiar vetor Kb p/ K;
		if(Error)
			for(int i=0;i<N;i++)
				Xk[i]=Xkb[i];
	}
	return Xkb;
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
	Menor_Auto_Valor_Vetor(Int,M,Ab,e,v,Lb);
	l=Lb+x;
}

float** S_householder(int N, float **A, float** Vetores){
	float** H_fim = MatIdentidade(N);
	float **H, *v, *n, norm;
	float **AT = copia(N,N,A);
	v = (float*)malloc(sizeof(float)*N);
	norm=0;
	
	for(int j=0; j<(N-2);j++){
		for(int i=0;i<N;i++)
			if(i<=j)
				v[i]=0;
			else
				v[i]=AT[i][j];
		//impVet(N,v);
		norm = NormaVetor(N,v);
		n = copiaVet(N,v);
		n[j+1]-=norm;
		norm = NormaVetor(N,n);
		Escalar(N,n,1/norm);
		//impVet(N,n);

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
	printf("\n Valor K %d; \n", k);
	
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
	float** D = Diag_QR(N,500,T,Q);;
	float** Vet = mult(N,N,N,H,Q);
	copiaVal(N,N,Vet,Vetores);
	for(int i=0;i<N;i++)
		valores[i]=D[i][i]; 
	
}

void SVD(int M, int N, float** Mat, float** U, float** S, float** V){
	float** AT = Transposta(M,N,Mat);
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

int main()
{

	//float B[4][4]= {{18.324146,0.001371,0.000357,0.000001},{0.001370,8.959310,-0.831104,-0.000003},{0.000358,-0.831104,-7.037904,-0.000019},{0.000000,-0.000003,-0.000019,-2.245561}};
	//float A[4][4] = {{1,2,6,8},{2,5,7,-3},{6,7,9,4},{8,-3,4,3}};

	float A[4][3] = {{1,2,1},{1,1,1},{1,1,1},{1,1,1}};


    int M = 4;
	int N = 3;

    float **matA;

    matA = (float **)malloc(sizeof(float)*M);
    for (int i = 0; i < M; i++){
        matA[i] = (float *)malloc(sizeof(float)*(N));
		matA[i] = A[i];
	}
	//preenche(M,N,matA);
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
	
	printf("\nU*S*V:\n");
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

    system("pause");
    return 0;
}