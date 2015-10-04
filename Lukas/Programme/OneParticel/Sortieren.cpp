#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

int main(){
	fstream f,dimension,bla;
	int n;
	dimension.open("dimension.txt");
	dimension>>n;
    dimension.close();
    double *Eigenvalue= new double[n];
    double **R= new double *[n];
	for(int i=0;i<n;i++){
		R[i]= new double [n];
	}
    
    f.open("daten.txt");
    for(int i=0;i<n;i++){
    	f>>Eigenvalue[i];
	}
    f.close();
    
    bla.open("Eigenvector.txt");
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			bla >> R[i][j];
		}
	}
    bla.close();
    
    //finished Daten reading
    double k,Eigenwertnumerical;
    cout << "which Eigenvalue are you looking for?";
    cin >> k;
    int spalte=0;
    for(int i=0;i<n;i++){
    	if(Eigenvalue[i]>k-1&&Eigenvalue[i]<k+1){
    		Eigenwertnumerical=Eigenvalue[i];
    		spalte=i;
		}
	}

	ofstream psi("plot1.txt");
	for(int i=0;i<n;i++){
		psi << i+1<< "  "<<R[i][spalte]<<endl;
	}
	psi.close();
	
	for(int i=0;i<n;i++){
		cout << i+1<< "  "<<R[i][spalte]<<endl;
	}

	string teil1="plot [0:200][-0.2:0.2] \"plot1.txt\" using 1:2 title \"Radialfunktion with Energy= ";
	string teil2=" \" with lines \n	set xlabel \"XLABEL\" \n 	set ylabel \"u(r) -Radial part of the wave function \" \n	set term png \n	set output \"Radialfunktion.png\" \n	replot \n	set term wxt";
	
	ofstream file("plotskript.txt");
	
	file <<teil1<<Eigenwertnumerical<<teil2;
	
	file.close();
	
    system("start gnuplot plotskript.txt");
}
