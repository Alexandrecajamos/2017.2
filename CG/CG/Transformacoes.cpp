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
	n[3] = 0;
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
	pEscalar(2,N,N,H);
	H[N]=(float*)malloc(sizeof(float)*N);
	for(int i=0;i<=N;i++){
		H[i][N]=0;
		H[N][i]=0;
	}
	E=sub(N+1,N+1,E,H);
	return E;
}
float** Rotacao(int N, int E,float A){
	float** M = MatIdentidade(N+1);
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

float** Rotacao2z(int N, int E,float CosA){
	float** M = MatIdentidade(N+1);
	float S = -sqrt(1-(CosA*CosA));
	float C = CosA;
	printf("\nCa = %f",C);
	printf("\nSenA = %f",S);
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
	S = V[0]/H;
	C = L/H;
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

	float **Id = MatIdentidade(N+1);
	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++)
			Id[i][j]=Rot[i][j];
	return Id;
}
float** RQ(float A, float* V){
	A = A/2;
	float *q1, *q2, **Lq1, **Rq2, **rot;
	q1 = (float *)malloc(sizeof(float)*4);
	q2 = (float *)malloc(sizeof(float)*4);
	float mV=NormaVetor(4,V);
	//printf("Norma Vetor = %f\n",mV);
	for(int i=0;i<3;i++){
		float aux = (sin(A*PI/180)*(V[i]/mV));
		q1[i]= aux;
		q2[i]=-aux;
	}
	q1[3] = cos(A*PI/180);
	q2[3] = cos(A*PI/180);
	impVet(4,q1);
	impVet(4,q2);
	Lq1 = Mat_Left_Quat(q1);
	//imp(4,4,Lq1);
	Rq2 = Mat_Right_Quat(q2);
	//imp(4,4,Rq2);
	rot = mult(4,4,4,Lq1,Rq2);
	return rot;
}
float** RQ2(float CosA, float* V){

	float *q1, *q2, **Lq1, **Rq2, **rot;
	float CA = 0.5*(CosA+1);
	CA = sqrt(CA);
	//printf("\nCa = %f",CA);
	float SA = 0.5*(1-CosA);
	SA = sqrt(SA);
	//printf("Sa = %f\n",SA);

	q1 = (float *)malloc(sizeof(float)*4);
	q2 = (float *)malloc(sizeof(float)*4);
	float mV=NormaVetor(4,V);
	//printf("Norma Vetor = %f\n",mV);
	for(int i=0;i<3;i++){
		float aux = (SA*(V[i]/mV));
		q1[i]= aux;
		q2[i]=-aux;
	}
	q1[3] = CA;
	q2[3] = CA;
	printf("Quaternio e Conjugado\n\n\n");
	impVet(4,q1);
	printf("\n");
	impVet(4,q2);
	printf("\n Matrizes L e R \n");
	Lq1 = Mat_Left_Quat(q1);
	imp(4,4,Lq1);
	Rq2 = Mat_Right_Quat(q2);
	imp(4,4,Rq2);
	rot = mult(4,4,4,Lq1,Rq2);
	return rot;
}

struct Material{
	float Ar,Ag,Ab,Dr,Dg,Db,Er,Eg,Eb;
	int m;
};

struct Cor{
	float IR,IG,IB;
};
struct Ponto{
	float Coord[4];
};
struct Aresta{
	Ponto P1;
	Ponto P2;
};
struct Face{
	Ponto P1;
	Ponto P2;
	Ponto P3;
	Material M;
};
struct Obj {
	Ponto* Pontos;
	Aresta* A;
	Face* F;
	int QtdPontos, QtdArestas,QtdFaces;
	Ponto CentroCirc;
	float R;

};
struct Luz{
	Ponto Pos;
	Cor I;
};



Ponto PreencheP(float* V){
	Ponto P;
	for(int i=0; i< 4;i++)
		P.Coord[i]= V[i];	
	return P;
}
void CalcCirc(Obj* O){

	float MX=O->Pontos[0].Coord[0],MY=O->Pontos[0].Coord[1],MZ= O->Pontos[0].Coord[2];
	float mX=MX,mY=MY,mZ=MZ;

	for(int i=0;i<O->QtdPontos;i++){

		float x = O->Pontos[i].Coord[0];
		float y = O->Pontos[i].Coord[1];
		float z = O->Pontos[i].Coord[2];

		if(x<mX)
			mX = x;
		if(x>MX)
			MX = x;
		if(y<mY)
			mY = y;
		if(y>MY)
			MY = y;
		if(z<mZ)
			mZ = z;
		if(z>MZ)
			MZ = z;

	}
	//printf("Teste Círculo de obj, Maior x = %f, y =%f e z =%f, Menor x = %f, y = %f, z = %f", MX,MY,MZ,mX,mY,mZ); 
	float Centro[4] = {(MX+mX)/2,(MY+mY)/2,(MZ+mZ)/2,1};
	//impVet(4,Centro);
	O->CentroCirc = PreencheP(Centro);
	float dx = abs(MX)+abs(mX);
	float dy = abs(MY)+abs(mY);
	float dz = abs(MZ)+abs(mZ);
	float d = 0;
	if(dx>dy && dx>dz)
		d = dx;
	else if(dy>dz)
		d = dy;
	else 
		d = dz;
	O->R=d/2;
}

Ponto* VetorPontos(int N){
	Ponto* P = (Ponto *)malloc(sizeof(Ponto)*N);
	for(int i =0; i<N;i++){
		preencheVet(4,P[i].Coord);
	}
	return P;
}
Obj PreencheObj(int QtdP){
	Obj obj;
	Ponto* P = VetorPontos(QtdP);
	obj.Pontos=P;
	obj.QtdPontos = QtdP;
	return obj;
}
Obj CopiaObj(Obj A){
	Obj O;
	O.Pontos=(Ponto *)malloc(sizeof(Ponto)*A.QtdPontos);
	for(int i=0; i< A.QtdPontos;i++)
		O.Pontos[i] = PreencheP(A.Pontos[i].Coord);

	O.QtdPontos=A.QtdPontos;
	return O;
}
Obj Transforma(Obj O, float** M, int N){
	Obj R;
	R = CopiaObj(O);
	for(int i=0;i<=O.QtdPontos;i++){
		float** V = VetorColuna(N, R.Pontos[i].Coord);

		V = mult(N,N,1,M,V);

		for(int j=0;j<N;j++)
			R.Pontos[i].Coord[j]=V[j][0];
	}
	return R;
}
void ImpObj(Obj O){
	printf("\n Imprimindo %d pontos do objeto",O.QtdPontos);
	for(int i=0;i<O.QtdPontos;i++){
		printf("\n Ponto %d: ", i);
		for(int j=0;j<4;j++){
			printf("|%f|",O.Pontos[i].Coord[j]);
		}
	}printf("\n");
}


float intersecciona(Obj O, Ponto Pij, float* normal){
	float T = -1;
	float tmin = 999;
	///IMplementar
	//1 - Conferir se intersecciona Esfera;
	//2 - Percorrer Faces 
				//Descartar Faces em que o produto escalar do vetor unitário do Raio com o vetor unitário normal à Face for positivo
	
	CalcCirc(&O);
	float a = ProdutoEscalar(4,Pij.Coord,Pij.Coord);
	float b = ProdutoEscalar(4,Pij.Coord,O.CentroCirc.Coord);
	b = -2*b;
	float c = ProdutoEscalar(4,O.CentroCirc.Coord,O.CentroCirc.Coord);
	c = c-(O.R*O.R);
	float Delta = ((b*b)-4*a*c);
	

	if(Delta>=0){
	
		for(int i=0;i<O.QtdFaces;i++){
			Face F = O.F[i];

			float* v1 = subVetor(4,F.P2.Coord,F.P1.Coord);
			float* v2 = subVetor(4,F.P3.Coord,F.P1.Coord);
			float *nF = Normal(4,v1,v2);
			
			
			float o[4] = {0,0,0,1};
			float *nR = subVetor(4,Pij.Coord,o);
			float nRaio = NormaVetor(4,nR);
			Escalar(4,nR,1/nRaio);
			
			if(ProdutoEscalar(4,nR,nF)<0){

				float* V2 = subVetor(4,F.P3.Coord,F.P2.Coord);
				float** A = MatIdentidade(4);

				for(int k=0;k<3;k++){
					A[k][0]=v2[k];
					A[k][1]=V2[k];
					A[k][2]=Pij.Coord[k];
				}
				float* lamb = Gauss(4,A,F.P3.Coord);
				if(lamb[2]>=0  && lamb[2]<tmin){ 
					
					float l3 = 1-(lamb[0]+lamb[1]);
					if(lamb[0]>=0 && lamb[0]<=1 && lamb[1]>=0 && lamb[1]<=1 && l3>=0 && l3<=1){
						T = lamb[2];
						normal[0]=nF[0];
						normal[1]=nF[1];
						normal[2]=nF[2];
						normal[3]=nF[3];
					}
				}
			}

		}
	}

	
	return T;
}
struct Observador{
	float coord[4];
	float i[4];
	float j[4];
	float k[4];
};
Observador ObsCalc(float *PO, float* LooK_At,float* Avup){
	Observador Obs;
	float* K = subVetor(4,PO,LooK_At);
	float nK = NormaVetor(3,K);
	Escalar(3,K,1/nK);
	float *X = Normal(3,subVetor(3,Avup,PO),K);
	float *J = Normal(3,K,X);
	for(int i=0; i< 4;i++){
		Obs.coord[i]=PO[i];
		Obs.k[i]=K[i];
		Obs.i[i]=X[i];
		Obs.j[i]=J[i];
	}
	return Obs;
}
float** CAMtoWord(Observador O){
	float** M = MatIdentidade(4);
	for(int i=0;i<4;i++){
		M[i][0]=O.i[i];
		M[i][1]=O.j[i];
		M[i][2]=O.k[i];
		M[i][3]=O.coord[i];
	}
	return M;
}
float** WtoCam(Observador O){
	float** M = MatIdentidade(4);
	for(int i=0;i<3;i++){
		M[0][i]=O.i[i];
		M[1][i]=O.j[i];
		M[2][i]=O.k[i];
	}
	M[0][3]= -ProdutoEscalar(4, O.i,O.coord);
	M[1][3]= -ProdutoEscalar(4, O.j,O.coord);
	M[2][3]= -ProdutoEscalar(4, O.k,O.coord);

	return M;
}
struct JanelaVis{
	float d,W,H;
	int N,M;
};
struct Cenario{
	Observador O;
	Obj *Objetos;
	int QtdObjetos;
};
Ponto** PixelsCoord(JanelaVis J){
	float DX,DY;
	DX=J.W/J.M;
	DY=J.H/J.N;

	Ponto** Pix = (Ponto **)malloc(sizeof(Ponto)*J.N);
    for (int i = 0; i < J.N; i++)
        Pix[i] = (Ponto *)malloc(sizeof(Ponto)*J.M);
		
	for(int i=0;i<J.N;i++){
		float Yi= (J.H/2)-(DY/2)-(i*DY);
		for(int j=0;j<J.M;j++){
			float Xj = (-J.W/2)+(DX/2)+(j*DX);
			Pix[i][j].Coord[0]=Xj;
			Pix[i][j].Coord[1]=Yi;
			Pix[i][j].Coord[2]=-J.d;
			Pix[i][j].Coord[3]=1;
		}
		
	}
	return Pix;
}
void Ray(JanelaVis J, Obj O, Observador Ob){
	float DX,DY;
	DX=J.W/J.M;
	DY=J.H/J.N;

	
	Ponto Point;
	float C[4] = {3,3,3,1};
	Point = PreencheP(C);
	Cor b;
	b.IB=1;b.IG=1;b.IR=1;
	Luz L;
	L.I=b;
	L.Pos=Point;

	
	Cor P; P.IB=1;P.IG=1;P.IR=1;

	Cor A; A.IB=0.2;A.IG=0.2;A.IR=0.2;

	Ponto** Pix = (Ponto **)malloc(sizeof(Ponto)*J.N);
    for (int i = 0; i < J.N; i++)
        Pix[i] = (Ponto *)malloc(sizeof(Ponto)*J.M);
		
	for(int i=0;i<J.N;i++){
		float Yi= (J.H/2)-(DY/2)-(i*DY);
		for(int j=0;j<J.M;j++){
			float Xj = (-J.W/2)+(DX/2)+(j*DX);
			Pix[i][j].Coord[0]=Xj;
			Pix[i][j].Coord[1]=Yi;
			Pix[i][j].Coord[2]=-J.d;
			Pix[i][j].Coord[3]=1;
			float* n = (float *)malloc(sizeof(float)*4);
			float t = intersecciona(O, Pix[i][j],n);
	
			if(t != -1 && t>0.9){
				//printf("|t = %f |\n",t);
				//printf("Teste Luz : Pos x=%f,y=%f;z=%f \n", L.Pos.Coord[0],  L.Pos.Coord[1],  L.Pos.Coord[2]);
				//printf("Teste Luz : Cor IR=%f,IG=%f;IB=%f \n", L.I.IR,L.I.IG,L.I.IB);
				//printf("Teste Material: Ra R = %f, Rd g = %f, Re b = %f", O.F[0].M.Ar, O.F[0].M.Dg, O.F[0].M.Eb);
				//impVet(4,n);
				
				Escalar(4, Pix[i][j].Coord,t);

				Cor Ia, Id, Ie;
				
				Ia.IR = O.F[0].M.Ar*A.IR;
				Ia.IG = O.F[0].M.Ag*A.IG;
				Ia.IB = O.F[0].M.Ab*A.IB;
				
				float *l = subVetor(4, L.Pos.Coord, Pix[i][j].Coord);
				Escalar(4,l,1/NormaVetor(4,l));
				float e = ProdutoEscalar(4,n,l);

				Id.IR = (O.F[0].M.Dr*L.I.IR)*e; 
				Id.IG = (O.F[0].M.Dg*L.I.IR)*e;
				Id.IB = (O.F[0].M.Db*L.I.IR)*e;
				
				
				float* v = subVetor(4, Ob.coord, Pix[i][j].Coord);
				float* r = copiaVet(4,n);
				Escalar(4,r,2*e);
				r = subVetor(4,r,l);
				float e2 = ProdutoEscalar(4,v,r);
				e2 = pow(e2, O.F[0].M.m);

				Ie.IR = (O.F[0].M.Er*L.I.IR)*e2;
				Ie.IG = (O.F[0].M.Eg*L.I.IG)*e2;
				Ie.IB = (O.F[0].M.Eb*L.I.IB)*e2;
				


				P.IR = Ia.IR + Id.IR + Ie.IR;
				P.IG = Ia.IG + Id.IG + Ie.IG;
				P.IB = Ia.IB + Id.IB + Ie.IB;

				if(P.IR>1)
					P.IR=1;
				if(P.IG>1)
					P.IG=1;
				if(P.IB>1)
					P.IB=1;
			
			}else{
				P.IR = 0.1;
				P.IG = 0.1;
				P.IB = 0.1;
			}

		

			printf("Pixel[%d][%d] = R: %f, G: %f, B:%f \n",i,j, P.IR,P.IG,P.IB);


		}
		//printf("\n");
	}
	
}

int main(){

	//Define Observador: 
	
	float PO[4] = {10,11,14,1}; // Posição
	float LA[4] = {1,0,4,1}; // Look_At
	float Avup[4] = {7,8,5,1}; // Ponto auxiliar View Up; 
	
	Observador Obs = ObsCalc(PO,LA,Avup); //Calcula vetores i, j, k;

	float** Twc = WtoCam(Obs); // Matriz Word to Cam
	float** Tcw = CAMtoWord(Obs); // Cam to Word;
	imp(4,4,Twc);
	//float** conf = mult(4,4,4,Twc,Tcw); imp(4,4,conf);
	
	//Definir Janela de visualização

	JanelaVis J; 
	J.W=10; J.H=10;
	J.d=4;  
	
	J.M=1000; J.N=1000;

	Ponto** pixs = PixelsCoord(J);
	Escalar(3,pixs[20][5].Coord,20);
	//impVet(4,pixs[20][5].Coord);


	//Definição do Objeto
	

		//Definição do Vetor de pontos p/ objeto;

	Ponto*P;  P = (Ponto *)malloc(sizeof(Ponto)*4); //Quantidade de pontos p/ objeto;

	float P1[4] = {-5,-5,-4,1};
	P[0] = PreencheP(P1);
	float P2[4] = {5,5,-4,1};
	P[1] = PreencheP(P2);
	float P3[4] = {0,5,-20,1};
	P[2] = PreencheP(P3);
	P[3] = PreencheP(pixs[20][5].Coord);

	Obj obj; 
	obj.Pontos=P;
	obj.QtdPontos=4;
	CalcCirc(&obj);

	ImpObj(obj);

    float* e1 = subVetor(4, obj.Pontos[1].Coord,obj.Pontos[0].Coord); // Vetor P1P2
	float* e2 = subVetor(4, obj.Pontos[2].Coord,obj.Pontos[0].Coord); // Vetor P1P3
	//impVet(4,e1);
	//impVet(4,e2);
	float* n = Normal(4,e1,e2);
	impVet(4,n);
	float** Es = EspelhoArb(3,VetorColuna(4,n));
	float o[4] = {0,0,0,1};
	float* t1 = subVetor(4,o ,obj.Pontos[0].Coord);// Vetor OP1
	float **T1 = Translacao(3,t1);
	float **T2 = Translacao(3,obj.Pontos[0].Coord);
	float ** T = mult(4,4,4,T2,Es);
	T = mult(4,4,4,T,T1);
	obj = Transforma(obj,T,4);
	
	ImpObj(obj);

	imp(4,4,T1);
	imp(4,4,Es);
	imp(4,4,T2);



	//Transforma Objeto em coordenadas de câmera;
	/*
	Obj O_Cam = Transforma(obj,Twc,4);  

	//Definindo material padrão
	Material M;
	M.Ar = 0.2; M.Ag = 0; M.Ab = 0; M.Dr = 0; M.Dg = 1; M.Db = 0; M.Er = 0; M.Eg = 1; M.Eb = 0; M.m = 2;

	//Construção do Vetor de Faces do Objeto

	Face *faces = (Face *)malloc(sizeof(Face));

	Face F1;

	F1.P1=O_Cam.Pontos[0]; F1.P2=O_Cam.Pontos[3]; F1.P3=O_Cam.Pontos[1]; F1.M = M; 	faces[0]=F1;



	// Atribui vetor de faces ao Objeto;
	O_Cam.F=faces;
	O_Cam.QtdFaces=12;


	
	//Calcula coordenadas dos pontos na janela de visualização
		//Ponto** Pixs = PixelsCoord(J,Oc);

	float i,j;
	i=250; j=250;

	float DX,DY;
	DX=J.W/J.M;
	DY=J.H/J.N;

	Ponto Pixel;
	float Yi= (J.H/2)-(DY/2)-(i*DY);
	float Xj = (-J.W/2)+(DX/2)+(j*DX);
			
	Pixel.Coord[0]=Xj;
	Pixel.Coord[1]=Yi;
	Pixel.Coord[2]=-J.d;
	Pixel.Coord[3]=1;
	
	
	//impVet(4,Pixel.Coord);

	float * n = (float*)malloc(sizeof(float*)*4);

	float t = intersecciona(O_Cam,Pixel,n);
	
	if(t != -1){
		printf("\n\nDeu certo\n\n");
		printf("\nTeste t = %f " ,t);
		Escalar(4, Pixel.Coord,t);
	
	}
	
	
	 //impVet(4,n);




	/*
	Obj obj, objT; 
	
	Ponto*P; 
	P = (Ponto *)malloc(sizeof(Ponto)*3);
	float P1[4] = {0,0,0,1};// 349049 A = 9, B = 4; C = 0; D = 9; E = 4; F = 3; 
	P[0] = PreencheP(P1);
	float P2[4] = {0,0,11,1}; //9
	P[1] = PreencheP(P2);
	float P3[4] = {7,0,0,1}; //11
	P[2] = PreencheP(P3);
	float P4[4] = {0,1,0,1};//6
	P[3] = PreencheP(P4);
	
	obj.Pontos=P;
	obj.QtdPontos=4;
	ImpObj(obj);
	float r50 = sqrt(50.0);
	float ex,ey,ez;
	ex = r50/P3[0];ey = r50/P4[1];ez = r50/P2[2];
	float ve[4] = {ex,ey,ez,1};
	float **E = Escala(4,ve);
	imp(4,4,E);
	objT = Transforma(obj,E,4);
	ImpObj(objT);


	float* t1 = subVetor(4,P1,objT.Pontos[2].Coord);
	float **T1 = Translacao(3,t1);
	float P3B[4] = {300,50,0,1}; // A = 9, B = 4; C = 0; D = 9; E = 4; F = 3; 
	float* t2 = subVetor(4,P3B,P1);
	float **T2 = Translacao(3,t2);
	
	objT = Transforma(objT,T1,4);
	ImpObj(objT);

	float *v = subVetor(4,objT.Pontos[3].Coord,objT.Pontos[2].Coord);
	//impVet(4,v);
	//float *M = somaVetor(4,objT.Pontos[3].Coord,objT.Pontos[2].Coord);
	//Escalar(4,M,0.5);
	//float c45 = sqrt(2.0)/2;
	float v1[4] = {0,0,-1,0};//{c45,c45,0,0};
	//float *v2 = subVetor(4,objT.Pontos[1].Coord,M);
	//float nv2 = NormaVetor(4,v2);
	//Escalar(4,v2,1/nv2);
	float *F1 = subVetor(4,objT.Pontos[2].Coord,objT.Pontos[1].Coord);
	float *F2 = subVetor(4,objT.Pontos[3].Coord,objT.Pontos[1].Coord);
	float* v2 = Normal(4,F1,F2);
	float O = ProdutoEscalar(4, v1,v2);
	//printf("Cos O = %f",O);
	float** Q = RQ2(O,v);
	
	objT = Transforma(objT,Q,4);
	ImpObj(objT);
	
	
	float* a = subVetor(4,objT.Pontos[1].Coord,objT.Pontos[2].Coord);
	float* b = subVetor(4,P3B,objT.Pontos[2].Coord);
	float na = NormaVetor(4,a);
	float nb = NormaVetor(4,b);
	Escalar(4,a,1/na);
	Escalar(4,b,1/nb);
	float teta = ProdutoEscalar(4,a,b);
	printf("O = %f", teta);
	float **Rz = Rotacao2z(3,2, teta);
	
	objT = Transforma(objT,Rz,4);
	ImpObj(objT);
	

	//float** T = mult(4,4,4,T2,Rz);
	//T= mult(4,4,4,T,Q);
	//T = mult(4,4,4,T,T1);
	
	
	objT = Transforma(objT,T2,4);
	
	ImpObj(objT);

	printf("\n Matrizes em ordem T1, RQ, Rz, T2 \n");
	imp(4,4,T1);
	imp(4,4,Q);
	imp(4,4,Rz);
	imp(4,4,T2);
	

	//Questão 3 
	
	float* e1 = subVetor(4, objT.Pontos[1].Coord,objT.Pontos[0].Coord); // Vetor P1P2
	float* e2 = subVetor(4, objT.Pontos[3].Coord,objT.Pontos[0].Coord); // Vetor P1P4
	impVet(4,e1);
	impVet(4,e2);
	float* n = Normal(4,e1,e2);
	impVet(4,n);
	float** Es = EspelhoArb(3,VetorColuna(4,n));
	t1 = subVetor(4, P1,objT.Pontos[0].Coord);// Vetor OP1
	T1 = Translacao(3,t1);
	T2 = Translacao(3,objT.Pontos[0].Coord);
	float ** T = mult(4,4,4,T2,Es);
	T = mult(4,4,4,T,T1);
	objT = Transforma(objT,T,4);
	
	ImpObj(objT);

	imp(4,4,T1);
	imp(4,4,Es);
	imp(4,4,T2);
	imp(4,4,T);
	//
	*/
	
	
	system("pause");
    return 0;
}