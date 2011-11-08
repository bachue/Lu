#include <stdio.h>
#include <stdlib.h>
#include "lu.h"
#include "mpi.h"

void lu(float *A, float *X, int N, int rank)
{
	int i, j, k = N / 2;
	float t;

	if(rank == 0)
	{
		c(1) /= b(1);
		d(1) /= b(1);
		for(i = 2; i <= k - 1; i++)
		{
			t = b(i) - (a(i) * c(i - 1));
			c(i) /= t;
			d(i) = (d(i) - (a(i) * d(i - 1))) / t;
		}
	}
	else if(rank == 1)
	{
		a(N) /= b(N);
		d(N) /= b(N);
		for(i = N - 1; i >= k + 1; i--)
		{
			t = b(i) - (c(i) * a(i + 1));
			a(i) /= t;
			d(i) = (d(i) - (c(i) * d(i + 1))) / t;
		}
	}

	if(rank == 0)
	{
		float buf[2];
		MPI_Status status;
		MPI_Recv(buf, 2, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &status);
		a(k + 1) = buf[0];
		d(k + 1) = buf[1];

		t = b(k) - c(k) * a(k + 1) - a(k) * c(k - 1);
		d(k) = (d(k) - c(k) * d(k + 1) - a(k) * d(k - 1)) / t;
		x(k) = d(k);
	}
	else if(rank == 1)
	{
		float buf[2] = { a(k + 1), d(k + 1) };
		MPI_Request request;
		MPI_Isend(buf, 2, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &request);
	}


	if(rank == 0)
	{
		for(i = k - 1; i >= 1; --i)
		{
			x(i) = d(i) - c(i) * x(i + 1);
		}
	}
	else if(rank == 1)
	{
		for(i = k + 1; i <= N; ++i)
		{
			x(i) = d(i) - a(i) * x(i - 1);
		}
	}

	MPI_Bcast(X + k + 1, N - k, MPI_FLOAT, 1, MPI_COMM_WORLD);
}

int main(int argc, char *argv[])
{
	int size, rank, i;
	float a[3 * 4] = {
		1, 1, 1, 2,
		1, 1, 1, 3,
		1, 1, 1, 2
	};
	float x[3]; 

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(size != 2)
	{
		fprintf(stderr, "this program requires 2 processes.");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	lu(a, x, 3, rank);

	if(rank == 0)
	{
		printf("result = ");
		for(i = 0; i < 3; ++i)
		{
			printf("%d ", x[i]);
		}
		printf("\n");
	}

	MPI_Finalize();

	return 0;
}
