#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdexcept>
#include <fstream>
#include <string>
#include <unistd.h>
#include "header.h"
#include "mpi.h"

using namespace std;

void inverse(double* matrix, double* inverse, int n, int my_rank, int threads_count, double* A, double* B, int* return_flag)
{
    struct{
        double value;
        int   rank;
    } in, out;

	
    int i;
    int j;
    int k;
	int l;
    int rank_max;
    int loc_ind_max;
    int ind_max;
    int index;
    int begin;
    int* flag;
    double a;
    double x;
    double loc_max;
    int change_rank1;
    int change_rank2;
    int change_loc_i1;
    int change_loc_i2;
    int size = n / threads_count;
    if(my_rank < n % threads_count)
    {
    	++size;
    }
    
    double eps1 = norma(n, size, matrix)*1e-15;
	double eps;
	MPI_Allreduce(&eps1, &eps, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    MPI_Status status;
    
	
    for (j = 0; j < n; ++j)
    {
    	MPI_Barrier(MPI_COMM_WORLD);
        in.value = 0;
        in.rank = 0;
		
		begin = j / threads_count;
			if (my_rank < j % threads_count) ++begin;
		
        for (i = begin; i < size; ++i)
        {
            if (fabs(matrix[i*n+j]) > fabs(in.value))
            {
                in.value = fabs(matrix[i*n+j]);
                in.rank = i*threads_count + my_rank;
            }
        }

        MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        if (abs(out.value) < eps) //проверка матрицы на вырожденность
        {
		    *return_flag = 0;
		}
		if (*return_flag == 0){return;}
		
        change_rank1 = j % threads_count;
        change_rank2 = out.rank % threads_count;
        change_loc_i1 = j / threads_count;
        change_loc_i2 = out.rank / threads_count;


	    	if (my_rank == change_rank1)/////////////////////////////////////////// MATRIX
	    	{
				for (k = 0; k < n; ++k) 
				{
					A[k] = matrix[change_loc_i1*n + k];
				}
	    	}

			MPI_Bcast(A, n, MPI_DOUBLE, change_rank1, MPI_COMM_WORLD);

	    	if (my_rank == change_rank2)
	    	{
	    		for (k = 0; k < n; ++k) 
				{
					B[k] = matrix[change_loc_i2*n + k];
				}
				
				for (k = 0; k < n; ++k) 
				{
					matrix[change_loc_i2*n + k] = A[k];
				}
	    	}

			MPI_Bcast(B, n, MPI_DOUBLE, change_rank2, MPI_COMM_WORLD);

	    	if (my_rank == change_rank1)
	    	{
	    		for (k = 0; k < n; ++k)
				{
					matrix[change_loc_i1*n + k] = B[k];
				}
	    	}//////////////////////////////////////////
	  	
    	if (my_rank == change_rank1)/////////////////////////////////////////// INVERSE
    	{
				for (k = 0; k < n; ++k) 
				{
					A[k] = inverse[change_loc_i1*n + k];
				}
	    	}

			MPI_Bcast(A, n, MPI_DOUBLE, change_rank1, MPI_COMM_WORLD);

	    	if (my_rank == change_rank2)
	    	{
	    		for (k = 0; k < n; ++k) 
				{
					B[k] = inverse[change_loc_i2*n + k];
				}
				
				for (k = 0; k < n; ++k) 
				{
					inverse[change_loc_i2*n + k] = A[k];
				}
	    	}

			MPI_Bcast(B, n, MPI_DOUBLE, change_rank2, MPI_COMM_WORLD);

	    	if (my_rank == change_rank1)
	    	{
	    		for (k = 0; k < n; ++k)
				{
					inverse[change_loc_i1*n + k] = B[k];
				}
	    	}//////////////////////////////////////////
		
		
		if (my_rank == change_rank1) ///gauss metod
		{
			a = 1.0 / matrix[change_loc_i1*n + j];
			for (k = 0; k < n; ++k)
			{
				matrix[change_loc_i1*n + k] *= a;
				inverse[change_loc_i1*n + k] *= a;
				A[k] = matrix[change_loc_i1*n + k];
				B[k] = inverse[change_loc_i1*n + k];
			}
	}
	MPI_Bcast(A, n, MPI_DOUBLE, change_rank1, MPI_COMM_WORLD);
	MPI_Bcast(B, n, MPI_DOUBLE, change_rank1, MPI_COMM_WORLD);
        
        
		begin = j / threads_count;
		if(my_rank <= j % threads_count) ++begin;
	    for (i = begin ; i < size; ++i)
	    {
	    	x = matrix[i*n + j];
	    	for (k = 0; k < n; ++k)
	    	{
	    		matrix[i*n + k] -= x*A[k];
	    		inverse[i*n + k] -= x*B[k];
	    	}
	    }
    }

	
    for(j = n-1; j >= 0 ; --j)
    {
    	change_rank1 = j % threads_count;
    	change_loc_i1 = j / threads_count;
    	
		if (my_rank == change_rank1)
		{
			for (k = 0; k < n; ++k)
			{
				A[k] = matrix[change_loc_i1*n + k];
				B[k] = inverse[change_loc_i1*n + k];
			}
		}

		MPI_Bcast(A, n, MPI_DOUBLE, change_rank1, MPI_COMM_WORLD);
		MPI_Bcast(B, n, MPI_DOUBLE, change_rank1, MPI_COMM_WORLD);

		int end = j / threads_count;
		if (my_rank < j % threads_count) ++end;
		
        for(i = 0; i < end; ++i)
        {
        	x = matrix[i*n + j];
        	
        	
        	for (k = 0; k < n; ++k)
	    	{
	    		matrix[i*n + k] -= x*A[k];
	    		inverse[i*n + k] -= x*B[k];
	    	}
        }
    }
	
	MPI_Barrier(MPI_COMM_WORLD);
    
    return;
    
}

