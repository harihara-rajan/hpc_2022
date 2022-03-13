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


    int n_actual = 5;
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

    int chunks = n/size;
    int L [n][n]= {0}; 
    int L_sub [chunks][n] = {0};
    int n_links[n] = {0};
    int n_links_total[n] = {0};
    int n_links_total_sub [n]= {0};
    int end_process = (n_actual / chunks) ;
    int rows_end_process = n_actual % chunks;
    int e[n];
    int d[n];
    int ed_matrix[n][n];
    int ed_matrix_sub[chunks][n];
    int d_sub[n];
    int e_sub[chunks];

    MPI_Scatter(L, chunks*n, MPI_INT, L_sub, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i=0; i<chunks; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (rank ==0)
            {
                L_sub[i][j] = rand() %2;
            }
            else
            {
                if (rank==end_process && n!= n_actual && i>= rows_end_process)
                {
                    L_sub[i][j] = 0;
                }
                else if (rank > end_process)
                {
                    L_sub[i][j] = 0;
                } 
                else
                {
                    L_sub[i][j] = rand();
                    L_sub[i][j] = rand() % 2;
                }

            }
        }
    }

    MPI_Gather(L_sub,chunks*n, MPI_INT, L, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        cout << "L Matrix " << endl;
        for (int i=0; i<n_actual; i++)
        {
            for (int j=0; j<n_actual; j++)
            {   
                cout<< L[i][j] << " " ;
            }
            cout << endl ;
        }
    }
    cout << endl;

    // computing n vectors parallely 
    MPI_Scatter(L, n*chunks, MPI_INT, L_sub, n*chunks, MPI_INT, 0, MPI_COMM_WORLD);
    
    for (i=0; i<chunks; i++)
    {
        for (j=0; j<n; j++)
        {
            n_links [j] += L_sub[i][j];
        }
    }

    MPI_Reduce(& n_links, n_links_total, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

    if (rank==0)
    {
        cout << "n_links " << endl;
        for (i=0; i<n_actual; i++)
        {
            cout << n_links_total[i] << " ";
        }
        cout << endl;
    }
    /* Initialize e parallely*/
    MPI_Scatter(e, chunks, MPI_INT, e_sub, chunks, MPI_INT, 0, MPI_COMM_WORLD);

    for (i=0; i<chunks; i++)
    {
        e_sub[i] = 1;
    }

    MPI_Gather(&e_sub, chunks, MPI_INT, e, chunks, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank==0)
    {
        cout << "e vector " << endl;
        for(i=0; i<n_actual; i++)
        {
            cout << e[i] << " ";
        }
        cout << endl;
    }

    /* initialise d vector parallely using gather and scatter*/
    MPI_Scatter(&n_links_total, chunks, MPI_INT, n_links_total_sub, chunks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(&d, chunks, MPI_INT, d_sub, chunks, MPI_INT, 0, MPI_COMM_WORLD);

    for (i=0; i<chunks; i++)
    {
        if (n_links_total_sub[i] == 0)
        {
            d_sub[i] = 1;
        }
        else
        {
            d_sub[i]=0;
        }
    }

    MPI_Gather(& d_sub, chunks, MPI_INT, d, chunks, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank==0)
    {
        cout << "d vector " << endl;
        for (i=0; i<n_actual; i++)
        {
            cout << d[i] << " ";
        }
        cout << endl;
    }

    /*outer prodct of e and d vector parallel implementation*/

    MPI_Scatter(e, chunks, MPI_INT, e_sub, chunks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(&d, chunks, MPI_INT, d_sub, chunks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(&ed_matrix, chunks*n, MPI_INT, ed_matrix_sub, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);

    for (i=0; i<chunks; i++)
    {
        for (j=0; j<n_actual; j++)
        {
            ed_matrix_sub[i][j] = e_sub[i] * d_sub[i];
        }
    }

    MPI_Gather(&ed_matrix_sub, chunks*n, MPI_INT, ed_matrix, n*chunks, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank ==0)
    {
        cout << "ed Matrix" << endl;
        for (i=0; i<n_actual; i++)
        {
            for (j=0; j<n_actual; j++)
            {
                cout << ed_matrix[j][i] << " " ;
            }
            cout << endl;
        }
    }

    MPI_Finalize();

}
