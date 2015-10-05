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
#include <time.h>

using namespace std;

//void jacobi_method ( double ** A, double ** R, int n );

double maxoffdiag ( double ** A, int * k, int * l, int n );
void rotate ( double ** A, double ** R, int k, int l, int n );
void findeigenvalue(int n,double ** A, int *position,int option, double eigenvalue);
void createtextfile(int n, double h, double **R, int position);


int main()
{
    double eigenvalue = 0;
    int answer, position, option, loop = 1;
    int n;
    cout << "Select the total number of steps" << endl;
    cin >> n;

    double h;
    double rho_max;
    double *rho_i = new double [n-1];
    double *d = new double [n-1];
    double e;

    rho_max = 6;

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

    double **R = new double*[n-1];

    //assigning the matrix R
    for ( int i = 0; i < n-1; i++ )
    {
        R[i] = new double [n-1];
        for ( int j = 0; j < n-1; j++ )
        {
            if ( i == j )
            {
                R[i][j] = 1.0;
            }
            else
            {
                R[i][j] = 0.0;
            }
        }
    }


    //From here on the we start solving the equation

    int k, l;
    double epsilon = 1.0e-8;
    int max_number_iterations =  n * n * n;
    int iterations = n-2;                           //it has this value because of the next for-loop
    double max_offdiag = maxoffdiag ( A, &k, &l, n );

    clock_t start, finish;
    start = clock();


    //makes the rotation for every first non-diagonal Matrixelement
    for(int i = 0; i < n-2 ; i++)
    {
        k = i;
        l = i +1;
        rotate ( A, R, k, l, n );
    }

    //The Jacobi algorithm gets executed
    while ( fabs(max_offdiag) > epsilon && iterations < max_number_iterations )
    {

        rotate ( A, R, k, l, n );
        max_offdiag = maxoffdiag ( A, &k, &l, n );
        iterations++;

    }

    finish = clock();
    cout << "Needed time for the algorithm in seconds: " << ( ( finish - start ) /CLOCKS_PER_SEC ) << endl << endl;

    cout << "Iterations: " << iterations << endl << endl;


    /* //Test: Gives out every calculated eigenvalue
   cout << "Number of iterations: " << iterations << endl;

    cout << "Eigenvalues: \n";
    for(int i = 0; i < n-1; i++)
    {
        cout << A[i][i] << endl;
    }
    cout << endl;*/


    while(loop == 1)
    {

    cout << "Do yout want to\n 1: Find a specific Eigenvalue\n 2: Want to know the first 3 calculated eigenvalues\n 3: Give out all eigenvalues\n";
    cin >> answer;

    //Gives out your choosen eigenvalue and saves its eigenvector in a textfile
    if(answer == 1)
    {
        option = 1;
        findeigenvalue(n,A, &position,option,eigenvalue);
        createtextfile(n,h,R,position);
        cout << "The calculated value for this eigenvalue is: " << A[position][position] << endl << endl;
    }

    //Gives out the first 3 calculated eigenvalues
    if(answer == 2)
    {
        option = 2;
        for(int i = 3; i <= 11;i= i + 4)
        {
            eigenvalue = i;
            findeigenvalue(n,A, &position,option,eigenvalue);
            cout << "EV to " << i << " : " << A[position][position] << endl;

        }
        cout << endl;
    }

    //Gives out all calculated eigenvalues
    if(answer == 3)
    {
        for(int i = 0; i < n-1; i++)
        {
            cout << A[i][i] << endl;
        }
        cout << endl;

    }

    cout << "Do you want to do someshing else?\n 1: yes\n 2: no" << endl;
    cin >> loop;
    }


    delete [] rho_i;
    for (int i = 0; i < n-1; i++)
    {
        delete[] A[i];
    }
    delete[] A;
    for (int i = 0; i < n-1; i++)
    {
        delete[] R[i];
    }
    delete[] R;
    delete[] d;

    return 0;
}


//The function finds the position of that eigenvector in the matrix and saves it in the variable 'position'
void findeigenvalue(int n,double ** A, int *position, int option, double eigenvalue)
{

    double difference;

    difference = 2.0;

    if(option == 1)
    {
        cout << "Which eigenvalue would you like to find?" << endl << endl;
        cin >> eigenvalue;
        cout << endl;

        for(int i = 0; i < n-1; i++)
        {
            if(fabs(A[i][i] - eigenvalue) < difference)
            {
                difference = fabs(A[i][i] - eigenvalue);
                *position = i;
            }
        }
        if(difference >= 2)
        {
            cout << "The eigenvalue could not be found!" << endl;
        }
    }

    if(option == 2)
    {
        for(int i = 0; i < n-1; i++)
        {
            if(fabs(A[i][i] - eigenvalue) < difference)
            {
                difference = fabs(A[i][i] - eigenvalue);
                *position = i;
            }
        }
        if(difference >= 2)
        {
            cout << "The eigenvalue could not be found!" << endl;
        }
    }
}

//Creates a textfile which contains the eigenvector to the eigenvector you found in findeigenvalue
void createtextfile(int n,double h,double **R, int position)
{
    ofstream Textfile("Jacobi_eigenvector.txt");           //Saves eigenvectors in a Textfile

    for(int i = 0; i < n-1; i++)
    {
        Textfile << h*(i+1) << "    " << R[i][position] << endl;
    }
    Textfile.close();

    cout << "A file containing the Eigenvector to your eigenvalue was created!" << endl << endl;

    Textfile.close();

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
void rotate ( double ** A, double ** R, int k, int l, int n )
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
        // Finally, we compute the new eigenvectors
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
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
