/*
Jacobi's method for finding eigenvalues
eigenvectors of the symetric matrix A.
The eigenvalues of A will be on the diagonal
of A, with eigenvalue i being A[i][i].
The j-th component of the i-th eigenvector
is stored in R[i][j].
A: input matrix (n x n)
R: empty matrix for eigenvectors (n x n)
n: dimention of matrices
*/
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

//void jacobi_method ( double ** A, double ** R, int n );

double maxoffdiag ( double ** A, int * k, int * l, int n );
void rotate ( double ** A, int k, int l, int n );

// Setting up the eigenvector matrix
int main()
{
    ofstream table("matrixsize_against_iterations.txt");

    double e;
    double rho_max;
    int max;

    rho_max = 100;
    max = 300;  //Deciding the end of the for-loop

    for(int n = 2; n <= max+1; n++)
    {
    double h;
    double *rho_i = new double [n-1];
    double *d = new double [n-1];


    h = rho_max/n;
    e = -1/(h*h);

    //assigning roh_i`s
    for(int i = 0;i < n-1; i++)
    {
        rho_i[i] = (i+1) * h;
    }

    //assigning d[i]`s
    for(int i = 0; i < n-1; i++)
    {
        d[i] = (2/(h*h))+(rho_i[i]*rho_i[i]);
    }

    //assigning the matrix A
    double **A = new double *[n-1];
    for(int i = 0; i < n-1; i++)
    {
        A[i] = new double [n-1];
        for(int j = 0; j < n-1; j++)
        {
            A[i][j] = 0;
            if(i == j){ A[i][j] = d[i];}
            if(i - j == 1){ A[i][j] = e;}
            if(i - j == -1){ A[i][j] = e;}
        }
    }


    //From here on the we start solving the equation

    int k, l;
    double epsilon = 1.0e-8;
    int max_number_iterations =  n * n * n;
    int iterations = n-2;                           //it has this value because of the next for-loop
    double max_offdiag = maxoffdiag ( A, &k, &l, n );

    //makes the rotation for every first non-diagonal Matrixelement
    for(int i = 0; i < n-2 ; i++)
    {
        k = i;
        l = i +1;
        rotate ( A, k, l, n );
    }

    //The Jacobi algorithm gets executed
    while ( fabs(max_offdiag) > epsilon && iterations < max_number_iterations )
    {

        rotate ( A, k, l, n );
        max_offdiag = maxoffdiag ( A, &k, &l, n );
        iterations++;

    }
    cout << "Iterations: " << iterations << endl << endl;


    /* //Test: Gives out every calculated eigenvalue
   cout << "Number of iterations: " << iterations << endl;

    cout << "Eigenvalues: \n";
    for(int i = 0; i < n-1; i++)
    {
        cout << A[i][i] << endl;
    }
    cout << endl;*/

    table << n-1 << "     " << iterations << endl;


    delete [] rho_i;
    for (int i = 0; i < n-1; i++)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] d;
    }

    table.close();

    return 0;
}


// Function to find the maximum matrix element. Can you figure out a more
// elegant algorithm?
double maxoffdiag ( double ** A, int * k, int * l, int n )
{
    double max = 0.0;
    for ( int i = 0; i < n-1; i++ )
    {
        for ( int j = i + 1; j < n-1; j++ )
        {
            if ( fabs(A[i][j]) > max )
            {
                max = fabs(A[i][j]);
                *l = i;
                *k = j;
            }
        }
    }
    return max;
}


// Function to find the values of cos and sin and then to rotate
void rotate ( double ** A, int k, int l, int n )
{
    double s, c;
    if ( A[k][l] != 0.0 )
    {
        double t, tau;
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);
        if ( tau > 0 )
        {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else
        {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1+t*t);
        s = c*t;
    }
    else
    {
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];
    // changing the matrix elements with indices k and l
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
    A[k][l] = 0.0; // hard-coding of the zeros
    A[l][k] = 0.0;
    // and then we change the remaining elements
    for ( int i = 0; i < n-1; i++ )
    {
        if ( i != k && i != l )
        {
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }
    }
    return;
}

/*  // Cout of A for test

     for(int i = 0; i < n-1; i++)
    for(int j = 0; j < n-1; j++)
    {
       cout << A[i][j] << "     " ;

    }
    cout << endl;
} cout << endl << endl << endl; */
