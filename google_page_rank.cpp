/*
author: Hari 
email : harihararajanrm.cms@gmail.com
file created on 05/03/2022
Introduction to High Performance Computing and Optimisation
*/
#include <iostream>
#include <vector>
using namespace std; 
/* header files at the end move all these function to a header file and include the header file in this 
main file*/
void show_matrix(vector<vector<int>> matrix, int n);
void show_vector(vector<float> A, int n);
vector<int> f_get_n_links(vector<int> n_links, vector<vector<int>> matrix, int n);
vector<float> f_matrix_vector_multiplication(vector<vector<float>> mat, vector<float> r, int n);
vector<vector<float>> f_compute_Q_matrix(vector<vector<int>>L, vector<vector<float>>Q, vector<int> n);
vector<vector<float>> f_compute_P_matrix(vector<vector<float>>Q, vector<int>e, vector<int>d, vector<int> n_links) ; // returns matrix of nxn
vector<float> f_power_iteration(vector<vector<float>>P_matrix, vector<float> r, int num_power_iteration);
float f_compute_l1_norm(vector<float> q);
vector<float> f_vector_normalising(vector<float>q, float l1_norm);

int main(int argc, char *argv [])
{
    int n = 12;
    int num_power_iteration = 30;
    vector<vector<int>> matrix(n, vector<int>(n)); // matrix 
    vector<vector<float>> Q(n, vector<float>(n));
    vector<vector<float>> P(n, vector<float>(n));
    vector<int> n_links(n);
    vector<float> r (n);
    vector<float> r_f (n);
    vector<int> e (n);
    vector<int> d (n);
    /* Generating random L matrix*/
    for (int i=0; i<n ; i++)
    {
        for (int j=0; j<n; j++)
        {
            matrix[i][j] = rand() % 2 ;
        }
        r[i] =((double) rand() / (RAND_MAX)) ;
    }
    // matrix = {{0,1,0,0,0,0},{1,0,0,0,0,0},{0,1,0,0,1,1},{}};
    n_links = f_get_n_links(n_links, matrix, n);
    Q = f_compute_Q_matrix(matrix, Q, n_links);

    P = f_compute_P_matrix(Q, e, d, n_links);
        for (int i=0; i<n ; i++)
    {
        for (int j=0; j<n; j++)
        {
            cout <<P[i][j] << "  ";
        }
        cout << endl;
    }

    r_f = f_power_iteration(P, r, num_power_iteration);

    for (int i=0; i<n; i++)
    {
        cout << r_f[i] << " ";
    }
    cout<< endl;
}

void show_matrix(vector<vector<int>> matrix, int n)
{

    for (int i=0; i<n ; i++)
    {
        for (int j=0; j<n ; j++)
        {
            cout << matrix [i] [j] << "  ";
        }
        cout << endl;
    }
}

void show_vector(vector<float>A, int n)
{
    
    for (int i=0; i<n; i++)
    {
        cout << A[i] << " " ;
    }
    cout << endl;
}

vector<int> f_get_n_links(vector<int> n_links, vector<vector<int>> matrix, int n)
{
    for (unsigned int i = 0 ; i<n; i++)
    {
        for (unsigned int j = 0 ; j<n; j++)
        {
            n_links [i] += matrix[j][i];
        }
    }
    return n_links;
}

vector<vector<float>> f_compute_Q_matrix(vector<vector<int>>L, vector<vector<float>>Q, vector<int> n_links)
{
    int n = L.size();

    for (int i = 0; i<n ; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (n_links[j] == 0)
            {
                Q[i][j] = L[i][j]; 
            }
            else
            {
                Q[i][j] = (1/float(n_links[j])) * L[i][j];  
            }
        }
    }
    return Q;
}

vector<float> f_matrix_vector_multiplication(vector<vector<float>> mat, vector<float> r, int n)
{
    vector<float> res(n); 

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            res[i] += mat[i][j] * r[j];
        }
    }
    return res;
}

vector<vector<float>> f_compute_P_matrix(vector<vector<float>>Q, vector<int>e, vector<int>d, vector<int> n_links)  // returns matrix of nxn
{
    int n = e.size();
    vector<vector<float>> ed_matrix(n, vector<float>(n));
    vector<vector<float>> P_matrix(n, vector<float>(n));
    for (int i = 0; i < n; i++)
    {
        e[i] = 1;
        if (n_links[i]==0)
        {
            d[i] = 1;
        }
        else
        {
            d[i] = 0;
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            ed_matrix[i][j] = 1/float(n) * (e[i] * d [j]);
        }

    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            P_matrix[i][j] = Q[i][j] + ed_matrix[i][j];
        }
    }
    // P_matrix = Q + ed_matrix;
    return P_matrix;
}

float f_compute_l1_norm(vector<float> q)
{
    int n = q.size();
    float l1_norm ;
    for (int i=0; i<n; i++)
    {
        l1_norm += q[i];
    }

    return l1_norm; 
}
vector<float> f_power_iteration(vector<vector<float>>P_matrix, vector<float> r, int num_power_iteration)
{
    int n = r.size();
    float l1_norm;
    vector<float>q (r.size());
    vector<float> r_new(n);
    for (int i=0; i<num_power_iteration; i++)
    {
        q = f_matrix_vector_multiplication(P_matrix, r,n);
        l1_norm = f_compute_l1_norm(q);
        r_new = f_vector_normalising(q, l1_norm);
        r = r_new;
    }
    return r;
}
vector<float> f_vector_normalising(vector<float>q, float l1_norm)
{
    int n = q.size();
    vector<float> r_new(n);
    for (int i=0; i<n; i++)
    {
        r_new[i] = q[i]/l1_norm;
    }
    return r_new;
}