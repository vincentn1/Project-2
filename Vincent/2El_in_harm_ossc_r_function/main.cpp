/*
This Programm is for evaluating the eigenfunctions of
two electrons in a harmonic osscillator with Coulomb interaction between each other.
The eigenfunctions are depending on r.
We are calculating in eV and nm.
*/

#include <iostream>
#include <cmath>
#include <fstream>
#include "lib.h"

using namespace std;

//void jacobi_method ( double ** A, double ** R, int n );

void findeigenvalue(int n,double *A, int *position, double eigenvalue);
void findeigenvalue2(int n,double *A, int *positioni);
void createtextfile(int n, double h, double **R, int position);
double pythag(double a, double b);
void tqli(double *d, double *e, int n, double **z);

int main()
{
    //initializing physical constants
    double k;                       //k reflects the strength of the oscillator potential
    double coulombinteraction;      //coulombinteraction = beta * e^2
    double phyfactor;               //phyfactor = (hquer*hquer)/m

    phyfactor = 7.61997*1e-6;
    coulombinteraction = 1.44;
    k = 1;

    //initializing non-physical values
    double eigenvalue = 0;
    int answer, position, loop = 1;
    int positioni[3];
    int n;
    cout << "Select the total number of steps" << endl;
    cin >> n;

    double h;
    double rho_max;
    double *rho_i = new double [n-1];
    double *d = new double [n-1];
    double *e = new double [n-1];

    rho_max = 2;

    h = rho_max/n;

    //assigning e[i]`s
    for(int i = 0;i < n-1; i++)
    {
        e[i] = -phyfactor/(h*h);
    }

    //assigning roh_i`s
    for(int i = 0;i < n-1; i++)
    {
        rho_i[i] = (i+1) * h;
    }

    //assigning d[i]`s
    for(int i = 0; i < n-1; i++)
    {
        d[i] = ((2*phyfactor)/(h*h))+0.25*k*(rho_i[i]*rho_i[i])+coulombinteraction*1/rho_i[i];
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
    tqli(d, e, n-1, R);

    while(loop == 1)
    {

    cout << "Do yout want to\n 1: Find a specific Eigenvalue\n 2: Want to know the first 3 calculated eigenvalues\n 3: Give out all eigenvalues\n";
    cin >> answer;

    //Gives out your choosen eigenvalue and saves its eigenvector in a textfile
    if(answer == 1)
    {
        findeigenvalue(n,d, &position, eigenvalue);
        createtextfile(n,h,R,position);
        cout << "The calculated value for this eigenvalue is: " << d[position] << endl << endl;
    }

    //Gives out the first 3 calculated eigenvalues
    if(answer == 2)
    {
        findeigenvalue2( n, d, positioni);
        for(int i = 0; i < 3; i++)
        {
            cout << i+1 << ". eigenvalue: " << d[positioni[i]] << endl;

        }
        cout << endl;
    }

    //Gives out all calculated eigenvalues
    if(answer == 3)
    {
        for(int i = 0; i < n-1; i++)
        {
            cout << d[i] << endl;
        }
        cout << endl;

    }

    cout << "Do you want to do someshing else?\n 1: yes\n 2: no" << endl;
    cin >> loop;
    }


    delete [] rho_i;
    for (int i = 0; i < n-1; i++)
    {
        delete[] R[i];
    }
    delete[] R;
    delete[] d;
    delete[] e;

    return 0;
}


//The function finds the position of that eigenvector in the matrix and saves it in the variable 'position'
void findeigenvalue(int n,double *A, int *position, double eigenvalue)
{

    double difference;

    difference = 2.0;

    cout << "Which eigenvalue would you like to find?" << endl << endl;
    cin >> eigenvalue;
    cout << endl;

    for(int i = 0; i < n-1; i++)
    {
        if(fabs(A[i] - eigenvalue) < difference)
        {
            difference = fabs(A[i] - eigenvalue);
            *position = i;
        }
    }
    if(difference >= 2)
    {
        cout << "The eigenvalue could not be found!" << endl;
    }
}

void findeigenvalue2(int n,double *A, int *positioni)
{
    for(int i = 0; i < 3; i++)
    {
        positioni[i] = 1;
    }
    for(int i = 0; i < 3; i++)
    {
        if(i == 0)
        {
            for(int j = 0; j < n-1; j++)
            {
                if(A[j] < A[positioni[i]])
                {
                    positioni[i] = j;
                }

            }
        }
        else
        {
            for(int j = 0; j < n-1; j++)
            {
                if((A[j] < A[positioni[i]]) && (A[j] > A[positioni[i-1]]))
                {
                    positioni[i] = j;
                }

            }
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

double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
// End: function pythag(), (C) Copr. 1986-92 Numerical Recipes Software )%.


       /*
       ** The function
       **              gauleg()
       ** takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

//This is the algoritm for this programm
void tqli(double *d, double *e, int n, double **z)
{
register int   m,l,iter,i,k;
double         s,r,p,g,f,dd,c,b;

for(i = 1; i < n; i++) e[i-1] = e[i];
e[n] = 0.0;
for(l = 0; l < n; l++) {
   iter = 0;
   do {
      for(m = l; m < n-1; m++) {
         dd = fabs(d[m]) + fabs(d[m+1]);
         if((double)(fabs(e[m])+dd) == dd) break;
      }
      if(m != l) {
         if(iter++ == 30) {
            printf("\n\nToo many iterations in tqli.\n");
            exit(1);
         }
         g = (d[l+1] - d[l])/(2.0 * e[l]);
         r = pythag(g,1.0);
         g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
         s = c = 1.0;
         p = 0.0;
         for(i = m-1; i >= l; i--) {
            f      = s * e[i];
            b      = c*e[i];
            e[i+1] = (r=pythag(f,g));
            if(r == 0.0) {
               d[i+1] -= p;
               e[m]    = 0.0;
               break;
            }
            s      = f/r;
            c      = g/r;
            g      = d[i+1] - p;
            r      = (d[i] - g) * s + 2.0 * c * b;
            d[i+1] = g + (p = s * r);
            g      = c * r - b;
            for(k = 0; k < n; k++) {
               f         = z[k][i+1];
               z[k][i+1] = s * z[k][i] + c * f;
               z[k][i]   = c * z[k][i] - s * f;
            } /* end k-loop */
         } /* end i-loop */
         if(r == 0.0 && i >= l) continue;
         d[l] -= p;
         e[l]  = g;
         e[m]  = 0.0;
      } /* end if-loop for m != 1 */
   } while(m != l);
} /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */

 /*
 ** The function
 **                tred2()
 ** perform a Housholder reduction of a real symmetric matrix
 ** a[][]. On output a[][] is replaced by the orthogonal matrix
 ** effecting the transformation. d[] returns the diagonal elements
 ** of the tri-diagonal matrix, and e[] the off-diagonal elements,
 ** with e[0] = 0.
 ** The function is modified from the version in Numerical recipe.
 */


