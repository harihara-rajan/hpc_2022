#include<iostream>
#include <mpi.h>
using namespace std;

int main (int argc, char *argv[])
{
    int n_actual = 21;
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
    int e[n];
    int e_sub[chunks];

    /* Initialize e parallely*/
    MPI_Scatter(e, chunks, MPI_INT, e_sub, chunks, MPI_INT, 0, MPI_COMM_WORLD);

    for (i=0; i<chunks; i++)
    {
        e_sub[i] = 1;
    }

    MPI_Gather(&e_sub, chunks, MPI_INT, e, chunks, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank==0)
    {
        for(i=0; i<n_actual; i++)
        {
            cout << e[i] << " ";
        }
        cout << endl;
    }

    MPI_Finalize();

}