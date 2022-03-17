#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <time.h>
using namespace std; 

int main ( int argc, char *argv[] )
{
	time_t walltime = time(nullptr);
    // srand(time(NULL));

    int n_actual = 5;
    int i, j, rank, size, n;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // condition to check if number of webpages is exactly divisibly by size
    int rem = n_actual%size;
    if (rem!=0)
    {
        n = n_actual + (size-rem); // for padding if # webpages not exactly divisible by processes.
    } 

    else
    {
        n = n_actual;
    }   

    int chunks = n/size;           // 
    int num_power_iteration = 30; 
    
    int L [n][n]= {0};             // L matrix initialised with 0's
    int L_sub [chunks][n] = {0};   // L_sub matrix for scatter and gather operations 
    
    float Q [n][n]={0};            // Q matrix initialised with 0's
    float Q_sub[chunks][n]={0};    // Q_sub matrix for scatter and gather operations
    
    int e[n];                      // e vector 
    int e_sub[chunks]={0};         // e_sub vector for scatter and gather operations
    
    int d[n];                      // d vector
    int d_sub[n]={0};              // d_sub vector for scatter and gather operations

    float ed_matrix[n][n] = {0};   // ed Matrix initialised with 0's
    float ed_matrix_sub[chunks][n] = {0}; // ed_matrix for scatter and gather operations 
    
    float P [n][n]={0};            // P matrix initialised with 0's 
    float P_sub[chunks][n]={0};    // P_sub matrix for scatter and gather operations 
    
    int n_links[n] = {0};          // number of links from j'th webpage to i'th webpage
    int n_links_total[n] = {0};
    int n_links_total_sub [n]= {0};
    

    float res[n]={0};               // rank vector
    float res_sub[chunks] = {0};    // rank vector sub for scatter and gather operations 

    float vect[n] = {0};
    float vect_sub [chunks] = {0};

    float sum = 0;                 // l1 norm
    float sum_chunks=0;            // l1 norm of each chunk 

    int end_process = (n_actual / chunks) ; 
    int rows_end_process = n_actual % chunks;

    /*Initialising L matrix parallely*/
    MPI_Scatter(L, chunks*n, MPI_INT, L_sub, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);
    /*
    Scatter L matrix of shape [n][n] into L_sub of shape [chunks][n] to different processes 
    and initialise with random number (0 or 1) parallely and assemble the matrix together into
    L matrix using Gather operation
    */
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
                    L_sub[i][j] = rand() % 2;
                }

            }
        }
    }
    /*Gather operation collects L_Sub from every processes and assemble it to L matrix*/
    MPI_Gather(L_sub,chunks*n, MPI_INT, L, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) // Printing L matrix into the console.
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
    cout << endl;
    }

    // computing n vectors parallely 
    MPI_Scatter(L, n*chunks, MPI_INT, L_sub, n*chunks, MPI_INT, 0, MPI_COMM_WORLD);
    /*
    Scatter L matrix of shape [n][n] into L_sub of shape [chunks][n] to different processes. 
    within each process column sum is computed. Then by using Allreduce, column sums of each 
    chunks is added and then broadcasted to all other processes for future use. 
    */
    for (i=0; i<chunks; i++)
    {
        for (j=0; j<n; j++)
        {
            n_links [j] += L_sub[i][j];
        }
    }
    /*Allreduce --> perform sum over n_links from different process and send it to n_links_total and 
    finally broadcast it to all other processes.
    */
    MPI_Allreduce(& n_links, n_links_total, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

    if (rank==0) // printing n_links that represents connection from j'th webpage to i'th webpage
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
    /*
    e vector initially intialised with 0's and here chunks element of e vector of shape [n] 
    is scattered to e_sub. once scattered each processes will allocate 1's to the e_sub.
    */
    for (i=0; i<chunks; i++)
    {
        e_sub[i] = 1;
    }
    /*
    Gather operation collects e_sub from all processes and gather it to e vector.
    */
    MPI_Gather(&e_sub, chunks, MPI_INT, e, chunks, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank==0) // printing e vector to the console.
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
    /*
    d vector initially intialised with 0's and here "chunks" element of e vector of shape [n] 
    is scattered to d_sub. similary chunks element of n_link_total vector is also scattered to 
    n_links_total_sub. when particular element of n_links_total_sub is 1 then element corresponding to 
    d_sub will be assigned as 1
    */
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
    /*
    Gather operation gathers all d_sub and send it to d vector
    */
    MPI_Gather(& d_sub, chunks, MPI_INT, d, chunks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&d, n, MPI_INT, 0, MPI_COMM_WORLD); 


    if (rank==0) // printing d vector to console 
    {
        cout << "d vector " << endl;
        for (i=0; i<n_actual; i++)
        {
            cout << d[i] << " ";
        }
        cout << endl;
    }

    /*outer prodct of e and d vector parallel implementation

    MPI_Scatter(e, chunks, MPI_INT, e_sub, chunks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(&d, chunks, MPI_INT, d_sub, chunks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(&ed_matrix, chunks*n, MPI_INT, ed_matrix_sub, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);

    for (i=0; i<chunks; i++)
    {
        for (j=0; j<n; j++)
        {
            if (j < n_actual)
            {
                ed_matrix_sub[i][j] = (1/float(n_actual)) * e_sub[i] * d_sub[i];
            }
            else
            {
                ed_matrix_sub[i][j] = 0;
            }
        }
    }

    MPI_Gather(&ed_matrix_sub, chunks*n, MPI_INT, ed_matrix, n*chunks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ed_matrix, n*n, MPI_INT, 0, MPI_COMM_WORLD); */

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            ed_matrix[i][j] = (1/float(n_actual)) * (e[i] * d [j]); // performing outer product (serial programing)
        }

    }
    MPI_Bcast(&ed_matrix, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (rank ==0) // printing ed matrix in the console.
    {
        cout << "ed Matrix" << endl;
        for (i=0; i<n; i++)
        {
            for (j=0; j<n; j++)
            {
                cout << ed_matrix[i][j] << " " ;
            }
            cout << endl;
        }
    cout << endl;
    }
    /* parallelise Q matrix */
    MPI_Scatter(Q, chunks*n, MPI_FLOAT, Q_sub, chunks*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(L, chunks*n, MPI_INT, L_sub, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(&n_links_total, chunks, MPI_INT, n_links_total_sub, chunks, MPI_INT, 0, MPI_COMM_WORLD);


    for (i=0; i<chunks; i++)
    {
        for (j=0; j<n; j++)
        {

            if (j > n_actual)
            {
                Q_sub[i][j] =  0;
            }
            else if (n_links_total[j]==0)
            {
                Q_sub[i][j] =  0;
            }
            else
            {
                Q_sub[i][j] = (1/float(n_links_total[j])) * L_sub[i][j] ;
            }
        }
    }

    MPI_Gather(&Q_sub, chunks*n, MPI_FLOAT, Q, chunks*n, MPI_FLOAT, 0, MPI_COMM_WORLD);


    if (rank == 0)
    {
        cout << "Q Matrix" << endl;
        for (i=0; i<n_actual; i++)
        {
            for (j=0; j<n_actual; j++)
            {
                cout << Q[i][j] << " ";
            }
            cout << endl;
        }
    cout << endl;
    }

    /* P matrix implimentation parallely */
    MPI_Scatter(Q, chunks*n, MPI_FLOAT, Q_sub, chunks*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(ed_matrix, chunks*n, MPI_INT, ed_matrix_sub, chunks*n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(P, chunks*n, MPI_FLOAT, P_sub, chunks*n, MPI_FLOAT, 0, MPI_COMM_WORLD);


    for(i=0; i<chunks; i++)
    {
        for (j=0; j<n; j++)
        {
            P_sub[i][j] = Q_sub[i][j] + ed_matrix_sub[i][j];
        }
    }
    MPI_Gather(&P_sub, chunks*n, MPI_FLOAT, P, chunks*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (rank==0)
    {
        cout << "P Matrix"<< endl;
        for (i=0; i<n_actual; i++)
        {
            for (j=0; j<n_actual; j++)
            {
                cout << P[i][j] <<  "  ";
            }
            cout << endl;
        }
    cout << endl;
    }

    for(i=0; i<n; i++)
    {
        if (i<n_actual)
        {
            res[i] = 1/float(n_actual);
        }
        else
        {
            res[i] = 0;
        }
    }

    /* POWER ITERATION */
    for (int k=0; k< num_power_iteration; k++)
    {
        sum_chunks = 0;
        MPI_Scatter(P, chunks*n, MPI_FLOAT, P_sub, chunks*n, MPI_FLOAT, 0, MPI_COMM_WORLD );
        MPI_Scatter(vect, chunks, MPI_FLOAT, vect_sub, chunks, MPI_FLOAT, 0, MPI_COMM_WORLD );

        for (int i=0; i<chunks; i++) // 1 should be replaced with chunks in original programming
        {
            for (int j=0; j<n_actual;j++)
            {
                vect_sub[i] += P_sub[i][j] * res[j];
            }
        }
        
        MPI_Gather(&vect_sub, chunks, MPI_FLOAT, res, chunks, MPI_FLOAT, 0, MPI_COMM_WORLD );
        
        MPI_Scatter(res, chunks, MPI_FLOAT, res_sub, chunks, MPI_FLOAT, 0, MPI_COMM_WORLD);
        
        for (int j=0; j<chunks;j++)
        {
            sum_chunks += res_sub[j];
        }

        MPI_Allreduce(&sum_chunks, & sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

        MPI_Scatter(res, chunks, MPI_FLOAT, res_sub, chunks, MPI_FLOAT, 0, MPI_COMM_WORLD);
        for (int i=0; i<chunks; i++)
        {
            res_sub[i] = res_sub[i]/sum;
        }
        MPI_Gather(&res_sub, chunks, MPI_FLOAT, res, chunks, MPI_FLOAT, 0, MPI_COMM_WORLD );
    }

    if (rank ==0)
    {
        cout << "R vector" << endl;
        for (int j=0; j<n_actual;j++)
        {
            cout << res[j] << " "; 
        }
        cout << endl;
    }

    MPI_Finalize();

}
