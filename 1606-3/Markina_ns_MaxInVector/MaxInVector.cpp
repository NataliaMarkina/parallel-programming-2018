#include "mpi.h"
#include <iostream>
#include <time.h>

using namespace std;


int FindMax(int* vector, int size, int firstNumber) {
	int res = vector[firstNumber];
	for (int i = firstNumber + 1; i < firstNumber + size; i++) {
		if (res < vector[i]) {
			res = vector[i];
		}
	}
	return res;
}

int main(int argc, char **argv) {
	int ProcRank, ProcNum, n;
	int Max = -1, ProcMax = -1;
	int *vector;
	double startTime, endTime;

	if (argc == 2)
		n = atoi(argv[1]);
	else
		n = 1000;
	vector = new int[n];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int k = n / ProcNum;
	int *part = new int[k];

	if (ProcRank == 0)
	{
		srand(time(0));

		for (int i = 0; i < n; i++)
		{
			vector[i] = rand();
		}

		if (n < 50)
		{
			cout << endl;
			for (int i = 0; i < n; i++) {
				cout << vector[i] << " ";
			}
			cout << endl << endl;
		}

		startTime = MPI_Wtime();
		Max = FindMax(vector, n, 0);
		endTime = MPI_Wtime();
		cout << endl << endl << "One process: max = " << Max << ". Search time = " << endTime - startTime << endl << endl;
	}

	startTime = MPI_Wtime();

	MPI_Scatter(vector, k, MPI_INT, part, k, MPI_INT, 0, MPI_COMM_WORLD);

	ProcMax = FindMax(part, k, 0);

	MPI_Reduce(&ProcMax, &Max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
	{
		if (n % ProcNum != 0)
		{
			for (int i = k * ProcNum; i < n; i++)
			{
				if (Max < vector[i])
				{
					Max = vector[i];
				}
			}
		}

		endTime = MPI_Wtime();

		cout << endl << endl << ProcNum << " process: max = " << Max << ". Search time = " << endTime - startTime << endl << endl;
	}

	MPI_Finalize();

	return 0;
}