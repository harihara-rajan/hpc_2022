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


    int n_actual = 18;
    int i, j, rank, size, n;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rem = n_actual%size;

    if (rem!=0)
    {
        n = n_actual + (size-rem);
    } 

    else
    {
        n = n_actual;
    }   

    cout << "n =  " << n  << endl;
    int chunks = n/size;

    int L [n][n]= {0}; 
    int L_sub [chunks][n];

    MPI_Scatter(L, chunks*n, MPI_INT, L_sub, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i=0; i<chunks; i++)
    {
        for (int j=0; j<n_actual; j++)
        {
            if (rank == 0)
            {
                L_sub[i][j] = rand() %2;
            }
            else
            {
                L_sub[i][j] = rand();
                L_sub[i][j] = rand() % 2;

            }
        }
    }

    MPI_Gather(L_sub,chunks*n, MPI_INT, L, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        for (int i=0; i<n_actual; i++)
        {
            for (int j=0; j<n_actual; j++)
            {   
                cout<< L[i][j] << " " ;
            }
            cout << endl ;
        }
    }
        MPI_Finalize();

}
