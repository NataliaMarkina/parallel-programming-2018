// Lab2_my_mpi_reduce.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <mpi.h>
#include <time.h>

using namespace std;


void mpi_maxTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			if (((int *)recvbuf)[i] < ((int *)sendbuf)[i])
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
		}
	}
	if (type == MPI_FLOAT)
	{
		for (int i = 0; i < count; i++)
		{
			if (((float *)recvbuf)[i] < ((float *)sendbuf)[i])
				((float *)recvbuf)[i] = ((float *)sendbuf)[i];
		}
	}
	if (type == MPI_DOUBLE)
	{
		for (int i = 0; i < count; i++)
		{
			if (((double *)recvbuf)[i] < ((double *)sendbuf)[i])
				((double *)recvbuf)[i] = ((double *)sendbuf)[i];
		}
	}
}

void mpi_minTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			if (((int *)recvbuf)[i] >((int *)sendbuf)[i])
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
		}
	}
	if (type == MPI_FLOAT)
	{
		for (int i = 0; i < count; i++)
		{
			if (((float *)recvbuf)[i] >((float *)sendbuf)[i])
				((float *)recvbuf)[i] = ((float *)sendbuf)[i];
		}
	}
	if (type == MPI_DOUBLE)
	{
		for (int i = 0; i < count; i++)
		{
			if (((double *)recvbuf)[i] >((double *)sendbuf)[i])
				((double *)recvbuf)[i] = ((double *)sendbuf)[i];
		}
	}
}

void mpi_sumTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			((int *)recvbuf)[i] = ((int *)recvbuf)[i] + ((int *)sendbuf)[i];
		}
	}
	if (type == MPI_FLOAT)
	{
		for (int i = 0; i < count; i++)
		{
			((float *)recvbuf)[i] = ((float *)recvbuf)[i] + ((float *)sendbuf)[i];
		}
	}
	if (type == MPI_DOUBLE)
	{
		for (int i = 0; i < count; i++)
		{
			((double *)recvbuf)[i] = ((double *)recvbuf)[i] + ((double *)sendbuf)[i];
		}
	}
}

void mpi_prodTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			((int *)recvbuf)[i] = ((int *)recvbuf)[i] * ((int *)sendbuf)[i];
		}
	}
	if (type == MPI_FLOAT)
	{
		for (int i = 0; i < count; i++)
		{
			((float *)recvbuf)[i] = ((float *)recvbuf)[i] * ((float *)sendbuf)[i];
		}
	}
	if (type == MPI_DOUBLE)
	{
		for (int i = 0; i < count; i++)
		{
			((double *)recvbuf)[i] = ((double *)recvbuf)[i] * ((double *)sendbuf)[i];
		}
	}
}

void mpi_landTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			((int *)recvbuf)[i] = ((int *)recvbuf)[i] && ((int *)sendbuf)[i];
		}
	}
	if (type == MPI_FLOAT)
	{
		for (int i = 0; i < count; i++)
		{
			((float *)recvbuf)[i] = ((float *)recvbuf)[i] && ((float *)sendbuf)[i];
		}
	}
	if (type == MPI_DOUBLE)
	{
		for (int i = 0; i < count; i++)
		{
			((double *)recvbuf)[i] = ((double *)recvbuf)[i] && ((double *)sendbuf)[i];
		}
	}
}

void mpi_bandTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			((int *)recvbuf)[i] = ((int *)recvbuf)[i] & ((int *)sendbuf)[i];
		}
	}
}

void mpi_lorTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			((int *)recvbuf)[i] = ((int *)recvbuf)[i] || ((int *)sendbuf)[i];
		}
	}
	if (type == MPI_FLOAT)
	{
		for (int i = 0; i < count; i++)
		{
			((float *)recvbuf)[i] = ((float *)recvbuf)[i] || ((float *)sendbuf)[i];
		}
	}
	if (type == MPI_DOUBLE)
	{
		for (int i = 0; i < count; i++)
		{
			((double *)recvbuf)[i] = ((double *)recvbuf)[i] || ((double *)sendbuf)[i];
		}
	}
}

void mpi_borTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			((int *)recvbuf)[i] = ((int *)recvbuf)[i] | ((int *)sendbuf)[i];
		}
	}
}

void mpi_lxorTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			((int *)recvbuf)[i] = ((int *)recvbuf)[i] != ((int *)sendbuf)[i];
		}
	}
	if (type == MPI_FLOAT)
	{
		for (int i = 0; i < count; i++)
		{
			((float *)recvbuf)[i] = ((float *)recvbuf)[i] != ((float *)sendbuf)[i];
		}
	}
	if (type == MPI_DOUBLE)
	{
		for (int i = 0; i < count; i++)
		{
			((double *)recvbuf)[i] = ((double *)recvbuf)[i] != ((double *)sendbuf)[i];
		}
	}
}

void mpi_bxorTree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type) {
	if (type == MPI_INT)
	{
		for (int i = 0; i < count; i++)
		{
			((int *)recvbuf)[i] = ((int *)recvbuf)[i] ^ ((int *)sendbuf)[i];
		}
	}
}

int log2_func(const int x)
{
	int y = x;
	int i = 0;
	while (y % 2 == 0)
	{
		i++;
		y = y / 2;
	}
	if ((int)(pow(2, (float)(i))) == x)
	{
		return i;
	}
	else
	{
		return 0;
	}
}

int Tree_reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm)
{
	int ProcRank, ProcNum, h;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	h = (int)(log((float)ProcNum) / log(2.0));
	if (log2_func(ProcNum) == 0) h++;
	
	for (int i = 0; i<h; i++)
	{
		int j = 0;
		while (j<ProcNum)
		{
			int k = j + (int)pow(2, (float)(i));

			if (k<ProcNum)
			{
				if (ProcRank == j)
				{
					MPI_Status status;
					MPI_Recv(recvbuf, count, type, k, 0, MPI_COMM_WORLD, &status);
					if (op == MPI_MIN) mpi_minTree(recvbuf, sendbuf, count, type);
					if (op == MPI_MAX) mpi_maxTree(recvbuf, sendbuf, count, type);
					if (op == MPI_SUM) mpi_sumTree(recvbuf, sendbuf, count, type);
					if (op == MPI_PROD) mpi_prodTree(recvbuf, sendbuf, count, type);
					if (op == MPI_LAND) mpi_landTree(recvbuf, sendbuf, count, type);
					if (op == MPI_BAND) mpi_bandTree(recvbuf, sendbuf, count, type);
					if (op == MPI_LOR) mpi_lorTree(recvbuf, sendbuf, count, type);
					if (op == MPI_BOR) mpi_borTree(recvbuf, sendbuf, count, type);
					if (op == MPI_LXOR) mpi_lxorTree(recvbuf, sendbuf, count, type);
					if (op == MPI_BXOR) mpi_bxorTree(recvbuf, sendbuf, count, type);
				}
				if (ProcRank == k)
				{
					MPI_Send(sendbuf, count, type, j, 0, MPI_COMM_WORLD);
				}
			}
			j = j + (int)pow(2, (float)(i + 1));

		}
	}
	for (int i = 0; i < count; i++)
		((int*)recvbuf)[i] = ((int*)sendbuf)[i];

	return 0;
}


int main(int argc, char *argv[])
{
	int n = atoi(argv[1]);
	int ProcRank, ProcNum;
	double TimeStartStandart, TimeFinishStandart, TimeStart, TimeFinish;
	int *mas = new int[n];
	int *mas_res = new int[n];

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	srand(time(0));

	for (int i = 0; i < n; i++)
	{
		mas[i] = rand() % 100 + ProcRank;
		if (n < 30)
		{
			cout << mas[i] << " ";
		}
	}
	cout << endl;

	for (int i = 0; i < n; i++)
	{
		mas_res[i] = 0;
	}

	TimeStartStandart = MPI_Wtime();

	MPI_Reduce(mas, mas_res, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	TimeFinishStandart = MPI_Wtime();

	if (ProcRank == 0)
	{
		if (n < 30)
		{
			cout << endl << "Standart realisation reduce - result: ";
			for (int i = 0; i < n; i++)
			{
				cout << mas_res[i] << " ";
			}
			cout << endl;
		}
	}

	TimeStart = MPI_Wtime();

	Tree_reduce(mas, mas_res, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	TimeFinish = MPI_Wtime();

	if (ProcRank == 0)
	{
		if (n < 30)
		{
			cout << endl << "Tree realisation reduce - result: ";
			for (int i = 0; i < n; i++)
			{
				cout << mas_res[i] << " ";
			}
			cout << endl;
		}
	}


	if (ProcRank == 0)
	{
		cout << endl;
		cout.setf(std::ios::fixed);
		cout.precision(10);
		cout << "Standart realisation reduce: t =  " << TimeFinishStandart - TimeStartStandart << endl;
		cout << "Tree realisation reduce: t =  " << TimeFinish - TimeStart << endl;
	}

	MPI_Finalize();

    return 0;
}

