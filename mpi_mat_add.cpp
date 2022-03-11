#include <iostream>
#include<mpi.h>
using namespace std; 

int main(int argc, char *argv [])
{

    int size, rank, i, j;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    float mat_2 [3][3] = {{0,1,1}, {1,1,1}, {1,0,1}};
    float Q [3][3] = {0};
    float n [3] = {3,0,5};
    float mat_2_sub[1][3];

    MPI_Scatter(mat_2, 1*3, MPI_FLOAT, mat_2_sub, 1*3, MPI_FLOAT, 0, MPI_COMM_WORLD);

    for (i=0; i<1; i++)
    {
        for (j=0; j<3; j++)
        {
            if (n[j] == 0)
            {
                mat_2_sub[i][j] = 0 * mat_2_sub[i][j];
            }
            else
            {
                mat_2_sub[i][j] = (1/float(n[j])) * mat_2_sub[i][j] ;
            }
        }
    }
    MPI_Gather(&mat_2_sub, 1*3, MPI_FLOAT, mat_2, 1*3, MPI_FLOAT, 0, MPI_COMM_WORLD);


    if (rank == 0)
    {
        for (i=0; i<3; i++)
        {
            for (j=0; j<3; j++)
            {
                cout << mat_2[i][j] << " ";
            }
            cout << endl;
        }
    }


    MPI_Finalize();
}