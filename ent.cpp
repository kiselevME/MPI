#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "mpi.h"
#include "header.h"
using namespace std;

double f(int, int, int, int);

double f(int k, int n, int i, int j)
{
	if (k == 1)
    {
    	return(n - std::max (i+1,j+1) + 1);
    }
    
    else if (k == 2)
    {
    	return(std::max (i+1,j+1));
    }
    
    else if (k == 3)
    {
    	return(abs(i-j));
    }
    
    else if (k == 4)
    {
    	return(1/(double)(i+j+1));
    }
    return(0);
}

int entermatrix(double* M, double* A, char* input, int n, int form, int my_rank, int threads_count)
{
	int i, j;
    int error = 0;
    int k = 0;
    int size_i;
    int dest_rank;
	int m;
    
    MPI_Status status;	

    size_i = n / threads_count;
    if(my_rank < n % threads_count)
    {
    	++size_i;
    }
    int a;
        
        if (form == 0)
        {
        	FILE *in;
            if (my_rank == 0)
            {
                in = fopen(input, "r");
                if (in == NULL)
                {
                    error = 1;
                }
            }
		    MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if(error==1) return 0;
		    	
			for (i = 0; i < n; ++i)
			{	
				if(my_rank==0)
				{
					for (j = 0; j < n; ++j) 
					a = fscanf(in, "%lf", &A[j]);
		                if (a != 1)
		                {
							fclose(in);
							error = 1;
							break;
		                }
				}
				MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
				if(error==1) return 0;
				
				dest_rank = i % threads_count;
				m = (i-my_rank)/threads_count;
				
				if (my_rank == 0 && dest_rank == 0)
				{
		            for(j = 0; j < n; ++j)
		            {
		                M[m*n + j] = A[j];
		            }
		        }
		        else if (my_rank == 0)
		        {
		            MPI_Send(A, n, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
		        }
		        else if (dest_rank == my_rank)
		        {	
		            MPI_Recv(A, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		            for(j = 0; j < n; ++j)
		            {
		                M[m*n + j] = A[j];
		            }
		        }
        	}
        	if(my_rank==0) fclose(in);
        }
        else
        {
        	for (i = 0; i < n; ++i)
        	{
        		int dest_rank;
				dest_rank = i % threads_count;
				if (dest_rank == my_rank)
				{
					for (j = 0; j < n; ++j)
					{
						M[k*n + j] = f(form, n, i, j);
					}
					++k;
				}
        	}
        }
        
	return 1;
}

void entermatrix_E(double* E, int n, int threads_count, int my_rank)
{
	int size_i;
	int i;
	int j;
	
    size_i = n / threads_count;
    if(my_rank < n % threads_count)
    {
    	++size_i;
    }

    for (int i = 0; i < size_i; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (my_rank + threads_count*i == j)
            {
            	E[i*n + j] = 1;
            }
            else
            {    
                E[i*n + j] = 0;
            }
        }
    }
}
