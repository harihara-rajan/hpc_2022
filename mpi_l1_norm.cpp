#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
using namespace std; 

int main ( int argc, char *argv[] )
{
   int rank;
   int n = 10000;
   int np = 4; 
   vector<float> q(n, 1.432);
   int chunk = n/np;
   vector<float> sums(chunk);
   float l1_norm;
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   sums[rank] = 0;
   for (int j=0; j<chunk; j++)
   {
      sums[rank]+= q[j];
   }

   
   MPI_Reduce(&sums[rank], &l1_norm,1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

   if (rank == 0)
   {
      cout << l1_norm << "  ";
      cout << endl;
   }

   MPI_Finalize();
   return 0;
}