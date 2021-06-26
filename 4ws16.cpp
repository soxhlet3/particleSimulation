#include <iostream>
#include <cmath>
using namespace std;

double** Einlesen (int n)   {
double** a= new double*[n];
for (int i=0; i<n; i++)   {
a[i] = new double[n];
}
for (int i=0; i < n; i++)   {
   for (int j=0; j < n; j++)   {
   cout << "[" << i << "]" << "[" << j << "]=";
   cin >> a[i][j];
   }
}   

return a;

}


void Zerstoeren (double** matrix, int n)   {
for (int i=0; i < n; i++)   {
delete[] matrix[i];
}
delete[] matrix;
}

double** StreicheSpalte (double** matrix, int n, int k)   {

double** untermatrix = new double*[n-1];
for (int i=0; i < n-1; i++) {
untermatrix[i] = new double[n-1];
}
for (int i=0;i < n-1; i++)   {
   for (int j=0; j < k-1; j++)   {
   untermatrix[i][j]=matrix[i+1][j];
   }
   for (int j=k+1; j < n; j++)   {
   untermatrix[i][j-1]=matrix[i+1][j];
   }
}
   
return untermatrix;   

}



double Determinante (double** matrix, int n)   {
if (n==1)   {
return matrix[0][0];
}  
if (n>1)  {
double det=0.0;
for (int j=0; j<n; j++)   {
det += pow(-1,j)*matrix[0][j]*Determinante(StreicheSpalte(matrix, n, j),n);
//Zerstoeren(StreicheSpalte(matrix,n,j),n-1);
}

return det;

}


}


int main()   {
int n;
cin >> n;
double** a= new double*[n];
for (int i=0; i<n; i++)   {
a[i] = new double[n];
}
a = Einlesen (n);
cout << Determinante(a,n);
}
