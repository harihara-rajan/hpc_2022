#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <time.h>
using namespace std; 

int main ( int argc, char *argv[] )
{
	time_t walltime = time(nullptr);
    srand(time(NULL));
    int rank;
    int n = 12;
    int np = 4; 
    int chunks = n/np;
    int L [n][n]; 
    int L_sub [chunks][n];
    MPI_Init(&argc,&argv);  
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Scatter(L, chunks*n, MPI_INT, L_sub, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i=0; i<chunks; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (rank == 0)
            {
                L_sub[i][j] = rand() %2;
            }
            else
            {
                L_sub[i][j] = rand()%1;

            }
        }
    }

    MPI_Gather(&L_sub[rank],chunks*n, MPI_INT, L_sub, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);
    // cout << "rank = " << rank << endl;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {   
            if (rank == 0)
            {
                L[i][j] = L_sub[i][j]; 
            }
            else
            {
                L[i][j] = L_sub[i][j]; 
            }
        }
    }
    if (rank == 0)
    {
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<n; j++)
            {   
                cout<< L[i][j] << " " ;
            }
            cout << endl ;
        }
    }
        MPI_Finalize();

}
