#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <unistd.h>
#include "mpi.h"
#include "header.h"

using namespace std;

double norma(int n, int size, double* matrix)
{
	double result, sum = 0, max = 0;
	
	
	int i, j, h;
	
	for (i = 0; i < size; ++i)
	{
		sum = 0;
		for (j = 0; j < n; ++j)
		{
			sum += matrix[i*n+j]*matrix[i*n+j];
		}
		if (sum > max)
		{
			max = sum;
		}
	}
	
	return sqrt(max);
}

double norm_matrix(double* matrix, double* inverse, int n, int my_rank, int threads_count, double* A)
{
    int i;
    int j;
    int k;
    int size;
    int dest_rank;
    
    double b, max = 0;
    double sum = 0;
    double elem;

    size = n / threads_count;
    if(my_rank < n % threads_count)
    {
    	++size;
    }

	
    for(i = 0; i < n; ++i)
    {
        dest_rank = i % threads_count;
    	
    	sum = 0;
    	for (j = 0; j < n; ++j)
    	{	
    		elem = 0;
    		if (my_rank == i % threads_count)
    		{
    			for (k = 0; k < n; ++k) A[k] = matrix[(i/threads_count)*n + k];
    		}
    		MPI_Bcast(A, n, MPI_DOUBLE, i % threads_count, MPI_COMM_WORLD);
    		
    		
    		for (k = 0; k < n; ++k)
			{
				if (k%threads_count == my_rank) elem += A[k]*inverse[(k/threads_count)*n + j];
			}
			MPI_Allreduce(&elem, &b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			if (i == j) b -= 1;// -E
			sum += b*b;
    	}

        if(sum > max)
        {
        	max = sum;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    max = sqrt(max);
    return max;
}
