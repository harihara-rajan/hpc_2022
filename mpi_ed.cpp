#include<iostream>
#include <mpi.h>
using namespace std;

int main (int argc, char *argv[])
{
    int n = 12;
    int i, j, rank, size;
    int e[n];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int chunks = n/size;
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
        for(i=0; i<n; i++)
        {
            cout << e[i] << " ";
        }
        cout << endl;
    }

    MPI_Finalize();

}