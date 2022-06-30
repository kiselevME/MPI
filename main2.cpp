#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <pthread.h>
#include <stdexcept>
#include <fstream>
#include <string>
#include "header.h"
#include <unistd.h>
#include "mpi.h"

using namespace std;


int main(int argn, char **args_) 
{
	int n;
	int m;
	int form;
	int threads_count;
	int my_rank;
	int max_rows;
	int mem_error = 0;
	int mem = 0;
	
	int* i_mas;
	int* j_mas;

	long double total_time;
	
    int return_flag = 1;
    int i, j, l, error;
    double *matrix;
    double *inverse_matrix;
    double *Copymatrix;
    double *A;
    double *B;
	double time = 0;
    
    MPI_Init(&argn, &args_);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &threads_count);
	
	if (my_rank == 0)// Проверки на некорректный ввод
	{	
		if (argn < 4 || argn > 5)
		{
			cout << "Некорректные входные данные";
			MPI_Finalize();
			return -1;
		}
		
		if (sscanf(args_[1], "%d", &n) != 1 || sscanf(args_[2], "%d", &m) != 1 || sscanf(args_[3], "%d", &form) != 1)
		{
			cout << "Некорректные входные данные";
			MPI_Finalize();
			return -1;
		}
		
		if (n < 1 || m < 0 || m > n)
		{
			cout << "Некорректные входные данные";
			MPI_Finalize();
			return -1;
		}
		
		if (!(((form == 0) && (argn == 5)) || (((0 < form) && (form < 5)) && (argn == 4))))
		{
			cout << "Некорректные входные данные";
			MPI_Finalize();
			return -1;
		}
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&form, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	max_rows = n / threads_count + n % threads_count;


	try
	{
		matrix = new double[max_rows * n];
	}
	catch (bad_alloc& e)
	{
		cout << "error: " << e.what() << endl;
		mem_error = 1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &mem, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (mem == 1)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}
	
	
	try
	{
		inverse_matrix = new double[max_rows * n];
	}
	catch (bad_alloc& e)
	{
		cout << "error: " << e.what() << endl;
		mem_error = 1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &mem, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (mem == 1)
	{
		delete[] matrix;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}
	
	
	try
	{
		Copymatrix = new double[max_rows * n];
	}
	catch (bad_alloc& e)
	{
		cout << "error: " << e.what() << endl;
		mem_error = 1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &mem, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (mem == 1)
	{
		delete[] matrix;
		delete[] inverse_matrix;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}
	
	try
	{
		A = new double[n];
	}
	catch (bad_alloc& e)
	{
		cout << "error: " << e.what() << endl;
		mem_error = 1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &mem, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (mem == 1)
	{
		delete[] matrix;
		delete[] inverse_matrix;
		delete[] Copymatrix;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}
	
	try
	{
		B = new double[n];
	}
	catch (bad_alloc& e)
	{
		cout << "error: " << e.what() << endl;
		mem_error = 1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &mem, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (mem == 1)
	{
		delete[] matrix;
		delete[] inverse_matrix;
		delete[] Copymatrix;
		delete[] A;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}
    
    if(entermatrix(matrix, A, args_[4], n, form, my_rank, threads_count) == 0)
		{

			if(my_rank==0) cout<<"initialization error"<<endl;
			delete [] matrix;
			delete [] inverse_matrix;
			delete [] Copymatrix;
			delete [] A;
			delete [] B;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			return -1;
		}
	
	int size_i = n / threads_count;
    if(my_rank < n % threads_count)
    {
    	++size_i;
    }
    for (i = 0; i < size_i; ++i)
    {
    	for (j = 0; j < n; ++j)
    	{
    		Copymatrix[i*n + j] = matrix[i*n +j];
    	}
    }
    		
	entermatrix_E(inverse_matrix, n, threads_count, my_rank);
    
    for(l=0; l<m; l++)
	{
        if(my_rank==l % threads_count)
        outmatrix(n, l/threads_count+1, m, l/threads_count, matrix);
		sleep(0.1);
        MPI_Barrier(MPI_COMM_WORLD);
    }

	if (my_rank == 0) cout << endl;

	sleep(1);

	if (my_rank == 0) time = MPI_Wtime();

	inverse(matrix, inverse_matrix, n, my_rank, threads_count, A, B, &return_flag);
	if (return_flag == 0)
	{
		if (my_rank == 0) cout << "MATRIX DEGENERATE";
		delete [] matrix;
		delete [] inverse_matrix;
		delete [] Copymatrix;
		delete [] A;
		delete [] B;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return 1;
	}

	if (my_rank == 0) time = MPI_Wtime() - time;

	MPI_Barrier(MPI_COMM_WORLD);

	for(l=0; l<m; l++)
	{
        if(my_rank==l % threads_count)
        outmatrix(n, l/threads_count+1, m, l/threads_count, inverse_matrix);
		sleep(0.1);
        MPI_Barrier(MPI_COMM_WORLD);
	}

	sleep(1);
	double norma = norm_matrix(Copymatrix, inverse_matrix, n, my_rank, threads_count,A);
	
	if (my_rank == 0) cout << endl;

	if (my_rank == 0) printf("%s : residual = %e elapsed = %.2f s = %d n = %d m = %d p = %d\n", args_[0], norma, time, form, n, m, threads_count);
    	
    delete [] matrix;
	delete [] inverse_matrix;
	delete [] Copymatrix;
	delete [] A;
	delete [] B;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 1;
}

