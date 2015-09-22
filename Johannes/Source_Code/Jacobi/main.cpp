#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;

//function that returns maximum (not on the diagonal!) and gives its adress
double offmaximum(double ** matrixA, int size, int * row, int * col)
{
    double max=0;
    for(int i=0; i<size-1; i++){
        for(int j=0; j<i; j++){
            if (fabs(matrixA[i][j])>max){
                max=fabs(matrixA[i][j]);
                *row=i;
                *col=j;
            }
        }
    }
    return max;
}


//actually rotates a given matrix around axis row and col, given cos and sin of the angle
void rotate(double ** toRotate, double ** eigenvalueMatrix, int row, int col, int size, double cosi, double sinu){
    double a_kk=toRotate[row][row];
    double a_ll=toRotate[col][col];
    //first, change those element with index k and/or l; zeros hardcoded
    toRotate[row][row]=cosi*cosi*a_kk - 2.0*cosi*sinu*toRotate[row][col] + sinu*sinu*a_ll;//a_kk
    toRotate[col][col]=sinu*sinu*a_kk + 2.0*cosi*sinu*toRotate[row][col] + cosi*cosi*a_ll;;//a_ll
    toRotate[row][col]=0.0;//a_kl
    toRotate[col][row]=0.0;//a_lk
    //change other elements
    for(int i=0; i<size; i++){
        if(i!=row && i!=col){
            double a_irow=toRotate[i][row];
            double a_icol=toRotate[i][col];
            toRotate[i][row]=cosi*a_irow-sinu*a_icol;
            toRotate[row][i]=toRotate[i][row];
            toRotate[i][col]=cosi*a_icol+sinu*a_irow;
            toRotate[col][i]=toRotate[i][col];
            //others unchanged
        }

        //Manipulate eigenvectors
        double ev_irow=eigenvalueMatrix[i][row];
        double ev_icol=eigenvalueMatrix[i][col];
        eigenvalueMatrix[i][row]=(cosi*ev_irow)-(sinu*ev_icol);
        eigenvalueMatrix[i][col]=(cosi*ev_icol)+(sinu*ev_irow);
    }
    return;
}

int main()
{
    //decline variables
    const int rhomin=0;
    int n, rhomax;
    double h;
    double tolerance;
    //ask for number of gridpoints n and rhomax
    cout << "How many gridpoints?" << endl;
    cin >> n;
    cout << "Maximum value of rho?" << endl;
    cin >> rhomax;
    //calculate and set steplength h
    h=(rhomax-rhomin)/(double)n;
    //decline vectors
    double * rho = new double [n+1];
    double * V = new double[n+1];
    //initialize vectors rho(=coordinate) and V (=potential)
    for(int i=0; i<n+1; i++){
        rho[i]=rhomin+(i*h);
        V[i]=(rho[i]*rho[i]);
    }
    //decline and initialize matrix A according to project instructions, as well as matrix R (used for the eigenvectors)
    double ** A = new double*[n-1];
    double ** R = new double*[n-1];
    //temps to calculate as much as possible outside of loop!
    double temp1=1/(double)(h*h);
    double temp2=2*temp1;
    for(int i=0;i<n-1;i++){
        A[i]=new double[n-1];
        R[i]=new double[n-1];
        for(int j=0;j<n-1;j++){
            if(j==i){
                A[i][j]=temp2+V[i+1];
                R[i][j]=1;
            }else{
                if(abs(i-j)==1){
                    A[i][j]=temp1;
                    R[i][j]=0;
                }else{
                    A[i][j]=0;
                    R[i][j]=0;
                }
            }
        }
    }
    //pointers that are used for the locations of maximum, and a double that contains the max
    int rowmax;
    int colmax;
    double off_max=999;
    double t, tau, c, s;
    //tolerance: The largest value, that an off-diag element is supposed to have at the end. To be sure, set a maximum number of iterations (~n^3)
    tolerance=0.00000001;
    int maxIterations=n*n*n;
    int isIterations=0;
    /*
     *
     * Main part starts here!
     *
     */
    while(off_max>tolerance && isIterations<maxIterations){
        //get new maximum and its position
        off_max=offmaximum(A, n-1, &rowmax, &colmax);
        //calculate tau, t, c and s
        tau = (A[colmax][colmax] - A[rowmax][rowmax])/(2*A[rowmax][colmax]);
        double root=sqrt(1.0 + tau*tau);
        if ( tau > 0 ) {
            t = 1.0/(tau + root);
        } else {
            t = -1.0/(-tau + root);
        }
        c = 1/sqrt(1+t*t);
        s = c*t;
        rotate(A, R, rowmax, colmax, n-1, c, s);
        isIterations++;
    }
    /*
     *
     * End of main part!
     *
     */
    //Find 5 lowest EVs(=ev*), their location (=l*) and write on .txt-file
    ofstream valfile, vecfile;
    string valfilename = to_string(n)+"_"+to_string(rhomax)+"_"+to_string(tolerance)+"vals"+".txt";
    string vecfilename = to_string(n)+"_"+to_string(rhomax)+"_"+to_string(tolerance)+"vecs"+".txt";
    valfile.open (valfilename);
    vecfile.open (vecfilename);
    int l5, l4, l3, l2, l1;
    double ev1=999;
    double ev2=999;
    double ev3=999;
    double ev4=999;
    double ev5=999;
    for(int i=0; i<n-1; i++){
        if(A[i][i]<ev5){
            if(A[i][i]<ev4){
                ev5=ev4;
                l5=l4;
                if(A[i][i]<ev3){
                    ev4=ev3;
                    l4=l3;
                    if(A[i][i]<ev2){
                        ev3=ev2;
                        l3=l2;
                        if(A[i][i]<ev1){
                            ev2=ev1;
                            l2=l1;
                            ev1=A[i][i];
                            l1=i;
                        }else{
                            ev2=A[i][i];
                            l2=i;
                        }
                    }else{
                        ev3=A[i][i];
                        l3=i;
                    }
                }else{
                    ev4=A[i][i];
                    l4=i;
                }
            }else{
                ev5=A[i][i];
                l5=i;
            }
        }
    }
    //print out vectors to seperate file
    for(int i=0; i<n-1; i++){
        vecfile << rho[i+1] << "\t" << R[i][l1] << "\t" << R[i][l2] << "\t" << R[i][l3] << "\t" << R[i][l4] << "\t" << R[i][l5] << endl;
    }
    valfile << l1 << ": " << ev1 << "\t" << l2 << ": " << ev2 << "\t" << l3 << ": " << ev3 << "\t" << l4 << ": " << ev4 << "\t" << l5 << ": " << ev5 << endl;
    valfile.close();
    vecfile.close();

    //free allocated memory!
    for(int i=0; i<n-1; i++){
        delete [] A[i];
        delete [] R[i];
    }
    delete [] A;
    delete [] R;
    delete [] rho;
    delete [] V;
    return 0;
}
