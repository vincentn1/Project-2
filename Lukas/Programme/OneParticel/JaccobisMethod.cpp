#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

double findmax(double **A , int * i,int * j, int n ){
	double maximum=0.0;
	for(int l=0;l<n;l++){
		for(int k=0;k<n;k++){
			if(k!=l){
				if(fabs(A[l][k])>maximum){
					maximum=fabs(A[l][k]);
					*i=l;
					*j=k;
				}
			}
		}
	}

	//cout << "Maximum OFF  "<<maximum<<"  "<<*i<<"  "<<*j<<endl;
	return(maximum);
}

void rotate(double **A,double **R, int i, int j, int n){
	//setting up different Variables sin/cos/tan for the matrix
	double s,c,t;
	double tau;
	if(A[i][j]!=0){
		tau=(A[j][j]-A[i][i])/(2*A[j][i]);
		if(tau>=0){
			t=1.0/(tau+sqrt(1+tau*tau));
		}else{
			t=-1.0/(-tau+sqrt(1+tau*tau));
		}
		c=1/(sqrt(1+t*t));
		s=t*c;
	}else{
		c=1.0;
		s=0;
		//This means no Rotation .. Theta=0;
	}
	
	double aii,ajj,alj,ali,rlj,rli;
	aii=A[i][i];
	ajj=A[j][j];
	A[i][i]=aii*c*c-2.0*A[i][j]*c*s+ajj*s*s;
	A[j][j]=ajj*c*c+2.0*A[i][j]*c*s+aii*s*s;
	A[i][j]=0;
	A[j][i]=0;
	for(int l=0;l<n;l++){
		if(l!=i&&l!=j){
			alj=A[l][j];
			ali=A[l][i];
			A[l][i]=ali*c-alj*s;
			A[i][l]=A[l][i];
			A[l][j]=alj*c+ali*s;
			A[j][l]=A[l][j];
		}
		//Eigenvektoren manipulieren
		rli=R[l][i];
		rlj=R[l][j];
		R[l][i]=c*rli-s*rlj;
		R[l][j]=c*rlj+s*rli;
		
	}
}

void JacobiMethod(double **A, double **R,int n){
	//fill up the matrix for the Eigenvectors 
	//choose orthonmormalsystem (Identity)
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(i==j){
				R[i][j]=1;
			}else{
				R[i][j]=0;
			}
		}
	}
	//choose epsilon
	double epsilon=pow(10,-8);
	int l,k;
	int iterationnummer=0,i_max=(double) n * (double) n * (double) n;
	double maximumoffdiag=findmax(A,&l,&k,n);
	//cout <<endl<<findmax(A,&l,&k,n)<<endl;
	while(fabs(maximumoffdiag)>epsilon&&iterationnummer<i_max){
		maximumoffdiag=findmax(A,&l,&k,n);
		rotate(A,R,l,k,n);
		iterationnummer++;
	}
	cout << "iterationsschritte"<<iterationnummer;
}

inline double V(int i,double rmin,double h){
	return((rmin+i*h)*(rmin+i*h));
}

int main(){
	int n=200;
	double **A= new double *[n];
	for(int i=0;i<n;i++){
		A[i]= new double [n];
	}
	
	double **R= new double *[n];
	for(int i=0;i<n;i++){
		R[i]= new double [n];
	}
	double r_min=0,r_max=6,h;
	h=(r_max-r_min)/n;
	//matrix fill up
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			A[i][j]=0;
		}
	}
	
	for(int i=0;i<n;i++){
		A[i][i]=(2/(h*h))+V(i+1,r_min,h);
	}
	
		for(int i=0;i<n-1;i++){
		A[i+1][i]=-1/(h*h);
		A[i][i+1]=A[i+1][i];
	}
	//ending of matrix fill up
	
	/*for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout << A[i][j]<<"   ";
		}
		cout<< endl;
	}*/

	cout << endl;
	
	JacobiMethod(A,R,n);
	
	cout << endl;
	//testing
	/*
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout << A[i][j]<<"   ";
		}
		cout<< endl;
	}
	cout << endl;
		for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout << R[i][j]<<"   ";
		}
		cout<< endl;
	}
	*/
	ofstream Zieldatei("Daten.txt"),Zieldatei1("Eigenvector.txt"),dimension("dimension.txt");

	for(int i=0;i<n;i++){
		Zieldatei << A[i][i]<<endl;
	}
	Zieldatei.close();
	
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			Zieldatei1 << R[i][j]<<"   ";
		}
		Zieldatei1<< endl;
	}	
	Zieldatei1.close();
	dimension<<n;
	dimension.close();
	
	
	
	//testing ende
	
	for(int i=0;i<n;i++){
		delete [] A[i];
	}
	delete [] A;
	
	for(int i=0;i<n;i++){
		delete [] R[i];
	}
	delete [] R;
}



