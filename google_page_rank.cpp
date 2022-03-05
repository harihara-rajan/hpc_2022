#include <iostream>
#include <vector>
using namespace std; 

void show_matrix(vector<vector<int>> matrix, int n);
void show_vector(vector<int> A, int n);
vector<int> f_get_n_links(vector<int> n_links, vector<vector<int>> matrix, int n);
vector<vector<float>> f_compute_Q_matrix(vector<vector<int>>L, vector<vector<float>>Q, vector<int> n);
int main()
{
    int n = 3;
    vector<vector<int>> matrix(n, vector<int>(n)); // matrix 
    vector<vector<float>> Q(n, vector<float>(n));
    vector<int> n_links(n);

    /* Generating random L matrix*/
    for (int i=0; i<n ; i++)
    {
        for (int j=0; j<n; j++)
        {
            matrix[i][j] = rand() % 2 ;
        }
    }
    show_matrix(matrix, n);
    n_links = f_get_n_links(n_links, matrix, n);
    show_vector(n_links, n);
    Q = f_compute_Q_matrix(matrix, Q, n_links);

        for (int i = 0; i <n ; i++)
    {
        for (int j=0; j<n; j++)
        {
            cout << Q[i][j] << " " ; 
        }
        cout << endl;
    }   
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

void show_vector(vector<int>A, int n)
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
    for (int i = 0; i <n ; i++)
    {
        for (int j=0; j<n; j++)
        {
            Q[i][j] = (1/float(n_links[j])) * L[i][j];
        }
    }
    return Q;
}
