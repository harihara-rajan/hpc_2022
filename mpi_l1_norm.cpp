#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
using namespace std; 

int main ( int argc, char *argv[] )
{
   int rank;
   int n = 12;
   int np = 4; 
   // vector<long double> q(n, 1.432);
   long double q[12] = {0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4};
   long double to_sum[3]; 
   long double sums[4] = {0};

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);


   MPI_Scatter(q, 3, MPI_LONG_DOUBLE, to_sum, 3, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
   
   for (int j=0; j<3; j++)
   {
      sums[rank]+= to_sum[j];
   }

   MPI_Gather(&sums[rank],1,MPI_LONG_DOUBLE,sums,1,MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);
   
   if (rank == 0)
   {
      for(int j=1; j<4; j++)
      {
         sums[0] += sums[j];
      }
      cout << sums[0] << endl;
   }

   MPI_Finalize();
   return 0;
}