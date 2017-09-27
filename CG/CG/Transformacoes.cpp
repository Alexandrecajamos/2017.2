#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <time.h>

#define PI 3.14159265

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
float** mult(int M, int N, int P, float **A, float **B)
{
    float **C;

    C = (float **)malloc(sizeof(float)*M);
    for (int i = 0; i < M; i++)
        C[i] = (float *)malloc(sizeof(float)*P);
	for (int i = 0; i < M; i++)
        for (int j = 0; j < P; j++) {
            C[i][j] = 0;
            for (int k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
               //printf("\ni = %d; j=%d; C[i][j] = %f", i,j, matC[i][j]);
        }
    return C;
}
void pEscalar(float E, int M, int N, float**A){
	for(int i=0;i<M;i++)
		for(int j=0;j<N;j++)
			A[i][j]*=E;
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
  for (int i = 0; i < c; i++)
	matA[i] = (float *)malloc(sizeof(float)*(l));
  for(int i=0; i<l;i++)
      for(int j=0; j<c;j++)
          matA[j][i]=A[i][j];
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
float** VetorColuna(int N, float* V){
	float** v;
	v = (float**)malloc(sizeof(float)*(N));
	for (int i = 0; i < N; i++){
        v[i] = (float *)malloc(sizeof(float));
		v[i][0] = V[i];
		
	}
	return v;
}

float NormaVetor(int N, float* V){
	float norma = 0;
	for(int i=0;i<N;i++)
		norma+=(V[i]*V[i]);
	norma = sqrt(norma);
	return norma;
}
void Escalar(int N, float* V, float Esc){
	for(int i=0;i<N;i++)
		V[i] *= Esc;
}
float* Normal(int N, float* A, float* B){
	float* n = (float *)malloc(sizeof(float)* N);
	n[0] = (A[1]*B[2])-(A[2]*B[1]); 
	n[1] = (A[2]*B[0])-(A[0]*B[2]);
	n[2] = (A[0]*B[1])-(A[1]*B[0]);
	float norma = NormaVetor(N,n);
	Escalar(N,n,(1/norma));
	return n;
}
//(Float Normal) Obtem vetor unitário perpendicular à 2 vetores;
float** Mat_Left_Quat(float* Q){
	float** q = MatIdentidade(4);
	for(int i=0;i<4;i++){
		q[i][i]=Q[3];
		q[3][i]=-Q[i];
		q[i][3]=Q[i];
	}
	q[0][1]= -Q[2];
	q[0][2]=  Q[1];
	q[1][0]=  Q[2]; 
	q[1][2]= -Q[0];
	q[2][0]= -Q[1]; 
	q[2][1]=  Q[0];
	
	return q;
}
float** Mat_Right_Quat(float* Q){
	float** q = MatIdentidade(4);
	for(int i=0;i<4;i++){
		q[i][i]=Q[3];
		q[3][i]=-Q[i];
		q[i][3]=Q[i];
	}
	q[0][1]=  Q[2];
	q[0][2]= -Q[1];
	q[1][0]= -Q[2]; 
	q[1][2]=  Q[0];
	q[2][0]=  Q[1]; 
	q[2][1]= -Q[0];
	
	return q;
}
//Matrizes de transformação
float** Escala(int N, float*V){
	float** E = MatIdentidade(N);
	for(int i=0; i<N;i++)
		E[i][i]=V[i];
	return E;
}
//Cisalhamento;
float** Cisalhamento(float A, int I1, int I2){
	float T = tan(A*PI/180);
	float** C = MatIdentidade(4);
	C[I2][I1] = T;
	return C;
}

float** Translacao(int N, float* V){
	float** T = MatIdentidade(N+1);
	for(int i=0;i<N;i++)
		T[i][N]=V[i];
	return T;
}
float** EspelhoArb(int N, float** V){
	float** E = MatIdentidade(N+1); //E = I -2n*nt
	float** nt = Transposta(N,1,V);
	float** H = mult(N,1,N,V,nt);
	pEscalar(-2,N,N,H); 
	H[N]=(float*)malloc(sizeof(float)*N);
	for(int i=0;i<=N;i++){
		H[i][N]=0;
		H[N][i]=0;
	}
	E=sub(N+1,N+1,E,H);
	return E;
}
float** Rotacao(int N, int E,float A){
	float** M = MatIdentidade(N);
	float S = sin(A*PI/180);
	float C = cos(A*PI/180);
	for(int i=0;i<N;i++){
		if(i!=E)
			for(int j=0;j<N;j++)
				if(j!=E)
					if(i!=j)
						if((i>E && i!=0 && j!=1) || (j<E && j!=0))
							M[i][j]=-S;
						else
							M[i][j]=S;	
						
					else
						M[i][j]=C;
	}
	return M;
}
float** Rotacao2(int N, int E,float A, float Sen, float Cos){
	float** M = MatIdentidade(N);
	float S = Sen;
	float C = Cos;
	for(int i=0;i<N;i++){
		if(i!=E)
			for(int j=0;j<N;j++)
				if(j!=E)
					if(i!=j)
						if((i>E && i!=0 && j!=1) || (j<E && j!=0))
							M[i][j]=-S;
						else
							M[i][j]=S;	
						
					else
						M[i][j]=C;
	}
	return M;
}
float** RotacaoArb(int N, float A, float* V){
	float** Rot;
	float **Rx, **Ry, **Rz;
	float S, C, L, H;
	float aux = ((V[1]*V[1])+(V[2]*V[2]));
	L=sqrt(aux);
	S = V[1]/L;
	C = V[2]/L;
	Rx = Rotacao2(N,0,A,S,C);
	printf("\nRotacao em x\n");imp(N,N,Rx);
	H = NormaVetor(N,V);
	S = L/H;
	C = V[0]/H;
	Ry = Rotacao2(N,1,A,S,C);
	Ry = Transposta(N,N,Ry);
	printf("\nRotacao em y-\n");imp(N,N,Ry);
	Rz = Rotacao(N,2,A);
	printf("\nRotacao em z\n");imp(N,N,Rz);
	Rot = mult(N,N,N,Ry,Rx);
	Rot = mult(N,N,N,Rz,Rot);
	Ry = Transposta(N,N,Ry);
	printf("\nRotacao em y\n");imp(N,N,Ry);
	Rx = Transposta(N,N,Rx);
	printf("\nRotacao em x-\n");imp(N,N,Rx);
	Rot = mult(N,N,N,Ry,Rot);
	Rot = mult(N,N,N,Rx,Rot);

	return Rot;
}
float** RQ(float A, float* V){
	float *q1, *q2, **Lq1, **Rq2, **rot;
	q1 = (float *)malloc(sizeof(float)*4);
	q2 = (float *)malloc(sizeof(float)*4);
	float mV=NormaVetor(3,V);
	//printf("Norma Vetor = %f\n",mV);
	for(int i=0;i<3;i++){
		float aux = (sin(A*PI/180)*(V[i]/mV));
		q1[i]= aux;
		q2[i]=-aux;
	}
	q1[3] = cos(A*PI/180);
	q2[3] = cos(A*PI/180);
	//impVet(4,q1);
	//impVet(4,q2);
	Lq1 = Mat_Left_Quat(q1);
	//imp(4,4,Lq1);
	Rq2 = Mat_Right_Quat(q2);
	//imp(4,4,Rq2);
	rot = mult(4,4,4,Lq1,Rq2);
	return rot;
}

int main(){

	int N = 4;
	float* p = (float *)malloc(sizeof(float)*N);
	preencheVet(N,p);
	float** C = Cisalhamento(30,1,0);
	
	imp(N,N,C);

	
	/*
	float p[4] = {2,2,2,1};
	float **P = VetorColuna(4,p);
	
	//Translação 1;
	float VT1[3] = {0,-10,0};
	float** T1 = Translacao(N,VT1);
	imp(N+1,N+1,T1);
	//Espelho
	float A[3] = {6,0,-4};
	float B[3] = {0,10,-4};
	float* n = Normal(N,A,B);
	n[N]=1;
	float **VN = VetorColuna(4,n);
	float** E = EspelhoArb(N,VN);
	imp(N+1,N+1,E);
	//Translação 2
	float VT2[3] = {0,10,0};
	float** T2 = Translacao(N,VT2);
	imp(N+1,N+1,T2);


	float** trans = mult(N+1,N+1,N+1,T1,E);
	trans = mult(N+1,N+1,N+1,trans,T2);
	imp(4,1,P);
	imp(N+1,N+1,trans);
	P = mult(4,4,1,trans,P);
	imp(4,1,P);
	*/

	system("pause");
    return 0;
}