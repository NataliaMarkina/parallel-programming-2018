// My_MPI_RadixSort.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream> 
#include <ctime>  
#include <mpi.h> 

using namespace::std;

#define ROOT 0

int ProcRank;
int ProcNum;

int* CreateArray(int size)
{
	int* arr = new int[size];

	srand((unsigned)time(NULL));

	for (int i = 0; i < size; i++)
		arr[i] = rand() % 1000 - 500;
	return arr;
}

void PrintArray(int* arr, int size)
{
	for (int i = 0; i < size; i++)
		cout << arr[i] << " ";
	cout << endl;
}

// Сравнение и обмен элементов
void Swap(int& a1, int& a2)
{
	int tmp = a1;
	a1 = a2;
	a2 = tmp;
}

void Radix(int byte, int N, int* source, int* dest)
{
	// *source - входной массив
	// *dest - отсортированный
	int count[256];
	int offset[256];
	memset(count, 0, sizeof(count));

	for (int i = 0; i < N; ++i)
	{
		if (byte == 3)
			count[((source[i] >> (byte * 8)) + 128) & 0xff]++;
		else
			count[((source[i]) >> (byte * 8)) & 0xff]++;
	}

	offset[0] = 0;
	for (int i = 1; i < 256; ++i)
		offset[i] = offset[i - 1] + count[i - 1];

	for (int i = 0; i < N; ++i)
	{
		if (byte == 3)
			dest[offset[((source[i] >> (byte * 8)) + 128) & 0xff]++] = source[i];
		else
			dest[offset[((source[i]) >> (byte * 8)) & 0xff]++] = source[i];
	}
}

void RadixSort(int *source, int N)
{
	int* temp = new int[N];
	Radix(0, N, source, temp);
	Radix(1, N, temp, source);
	Radix(2, N, source, temp);
	Radix(3, N, temp, source);
	delete[] temp;
}

void CalcDisplsAndWork(int* displs, int* CountWork, int size)
{
	int mid_workload = size / ProcNum;
	int ost = size % ProcNum;

	for (int i = 0; i < ost; ++i)
	{
		displs[i] = i * (mid_workload + 1);
		CountWork[i] = mid_workload + 1;
	}

	for (int i = ost; i < ProcNum; ++i)
	{
		displs[i] = mid_workload * i + ost;
		CountWork[i] = mid_workload;
	}
}

int main(int argc, char* argv[])
{
	int* MasRadixSeq = nullptr;
	int* MasRadixPp = nullptr;

	int* displs; // Массив смещений относительно начала буфера Array
	int* CountWork; // Массив кол-ва работы для каждого процесса
	int* rbuf;

	int size = 0;

	double StartSequentTime = 0;
	double StartParallelTime = 0;

	double EndSequentTime = 0;
	double EndParallelTime = 0;

	double TotalSequentTime = 0;
	double TotalParallelTime = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == ROOT)
	{
		cout << endl << "Enter size of array: " << endl;
		cin >> size;

		MasRadixSeq = new int[size];

		srand((unsigned)time(NULL));

		for (int i = 0; i < size; i++)
			MasRadixSeq[i] = rand() % 1000 - 500;

		MasRadixPp = new int[size];

		for (int i = 0; i < size; i++)
		{
			MasRadixPp[i] = MasRadixSeq[i];
		}

		if (size < 50)
		{
			cout << endl << "Unsorted array: " << endl;
			for (int i = 0; i < size; i++)
				cout << MasRadixSeq[i] << " ";
			cout << endl;
		}

		// Сортировка последовательной версии
		StartSequentTime = MPI_Wtime();
		RadixSort(MasRadixSeq, size);
		EndSequentTime = MPI_Wtime();
		TotalSequentTime = EndSequentTime - StartSequentTime;

		if (size < 50)
		{
			cout << endl << "Sorted array - sequential version: " << endl;
			for (int i = 0; i < size; i++)
				cout << MasRadixSeq[i] << " ";
			cout << endl;
		}
	}

	//Параллельная версия
	if (ProcRank == ROOT)
		StartParallelTime = MPI_Wtime();

	CountWork = new int[ProcNum];
	displs = new int[ProcNum];

	// Передаем размер от корневого процесса всем процессам
	MPI_Bcast(&size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	// Подсчитываем массив смещения и кол-ва работы
	CalcDisplsAndWork(displs, CountWork, size);

	int N = CountWork[ProcRank];

	rbuf = new int[N];

	int count[256];
	int offset[256];
	int GCount[256]; //Global count

	MPI_Status status;

	for (int byte = 0; byte < 4; ++byte)
	{
		//	Делим массивы на блоки
		MPI_Scatterv(MasRadixPp, CountWork, displs, MPI_INT, rbuf, N, MPI_INT, ROOT, MPI_COMM_WORLD);

		for (int i = 0; i < 256; i++)
			count[i] = 0;

		for (int i = 0; i < N; ++i)
		{
			if (byte == 3)
				count[((rbuf[i] >> (byte * 8)) + 128) & 0xff]++;
			else
				count[((rbuf[i]) >> (byte * 8)) & 0xff]++;
		}

		MPI_Allreduce(count, GCount, 256, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		offset[0] = 0;
		for (int i = 1; i < 256; ++i)
			offset[i] = offset[i - 1] + GCount[i - 1];

		int* SizesOfMas = new int[ProcNum];

		// Зная offset, можно создать массивы для всех процессов определенного размера
		int** NewMas;
		int CurrSizeMas; // максимальное число элементов в массиве NewMas[i]
		NewMas = new int*[ProcNum];
		for (int i = 0; i < ProcNum; i++)
		{
			CurrSizeMas = 0;

			for (int j = i * 256 / ProcNum; j < (i + 1) * 256 / ProcNum; j++)
				CurrSizeMas += GCount[j];

			if (CurrSizeMas > 0)
			{
				NewMas[i] = new int[CurrSizeMas];
				SizesOfMas[i] = CurrSizeMas;
			}
			else
			{
				NewMas[i] = NULL;
				SizesOfMas[i] = 0;
			}
		}

		int* index = new int[ProcNum];  // значение index[i] определяет положение элемента rbuf[i] в массиве NewMas[j]/количество элементов, который вошли в NewMas[i](размер)
		for (int i = 0; i < ProcNum; i++)
			index[i] = 0;

		for (int i = 0; i < N; ++i)
		{
			int ind; // байт с номером byte у элемента rbuf[i]
			if (byte == 3)
				ind = ((rbuf[i] >> (byte * 8)) + 128) & 0xff;
			else
				ind = ((rbuf[i]) >> (byte * 8)) & 0xff;

			//И смотрим на элемент который поступил
			int j = 0;
			for (; j < ProcNum; ++j)
				if (j * 256 / ProcNum <= ind && ind < (j + 1) * 256 / ProcNum)
					break;

			//Если он в j диапазоне массива смещений то располагаем элементы в j новый массив
			// Потом смотрим следующий элемент и так далее распределяем по массивам которые мы создали
			if (SizesOfMas[j] > 0)
			{
				NewMas[j][index[j]] = rbuf[i];
				index[j]++;
			}
		}

		// В этот массив соберем все размеры буферов других процессов с соответствующим диапазоном
		int * Gindex = new int[ProcNum];
		MPI_Alltoall(index, 1, MPI_INT, Gindex, 1, MPI_INT, MPI_COMM_WORLD);

		// Полученные буферы. Например, RecvBuffers[i] - буфер, полученный от i-ого процесса
		int** RecvBuffers = new int*[ProcNum];
		for (int i = 0; i < ProcNum; i++)
		{
			if (Gindex[i] > 0)
				RecvBuffers[i] = new int[Gindex[i]];
			else
				RecvBuffers[i] = NULL;
		}

		// Выполняем отсылку буферов другим процессам, ожидая так же, что эти процессы пришлют буферы и отправляющему процессу
		for (int i = 0; i < ProcNum; i++)
		{
			if (i != ProcRank)
			{
				MPI_Sendrecv(NewMas[i], index[i], MPI_INT, i, 1, RecvBuffers[i], Gindex[i], MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			}
		}

		// i-ый процесс объединяем все свои полученные куски массива и свой кусок, который он создал(NewMas[i]) в 1 массив tmpbuf в той последовательности, в который они шли в MasRadixPP
		int *tmpbuf;
		if (SizesOfMas[ProcRank] > 0)
		{
			tmpbuf = new int[SizesOfMas[ProcRank]];
			for (int i = 0, k = 0; i < ProcNum; i++)
			{
				if (i != ProcRank)
				{
					for (int j = 0; j < Gindex[i]; j++)
						tmpbuf[k++] = RecvBuffers[i][j];
				}
				else
				{
					for (int j = 0; j < index[i]; j++)
						tmpbuf[k++] = NewMas[i][j];
				}
			}
		}
		else tmpbuf = NULL;

		/* Распределяем по нужным местам(сортируем)*/
		int lim = 0; // смещение, чтобы записать на нужную позицию свеого массива
		for (int i = 0; i < ProcRank; i++)
			lim += SizesOfMas[i];

		int *dest = new int[SizesOfMas[ProcRank]];
		for (int i = 0; i < SizesOfMas[ProcRank]; i++)
		{
			if (byte == 3)
				dest[offset[((tmpbuf[i] >> (byte * 8)) + 128) & 0xff]++ - lim] = tmpbuf[i];
			else
				dest[offset[((tmpbuf[i]) >> (byte * 8)) & 0xff]++ - lim] = tmpbuf[i];
		}

		// Смещения для Gatherv
		int* displs_1 = new int[ProcNum];

		memset(displs_1, 0, sizeof(int)* ProcNum);
		displs_1[0] = 0;
		for (int i = 1; i < ProcNum; i++)
			displs_1[i] = displs_1[i - 1] + SizesOfMas[i - 1];

		// Собираем все в один массив. MasRadixPp будет отсортирован по байту byte
		MPI_Gatherv(dest, SizesOfMas[ProcRank], MPI_INT, MasRadixPp, SizesOfMas, displs_1, MPI_INT, ROOT, MPI_COMM_WORLD);

		delete[] displs_1;
		delete[] dest;
		delete[] tmpbuf;

		for (int i = 0; i < ProcNum; i++)
			delete[] RecvBuffers[i];

		delete[] RecvBuffers;
		delete[] Gindex;
		delete[] index;

		for (int i = 0; i < ProcNum; i++)
			delete[] NewMas[i];

		delete[] NewMas;
		delete[] SizesOfMas;
	}

	if (ProcRank == ROOT)
	{
		EndParallelTime = MPI_Wtime();
		TotalParallelTime = EndParallelTime - StartParallelTime;

		if (size < 50)
		{
			cout << endl << "Sorted array - parallel version:" << endl;
			PrintArray(MasRadixPp, size);
			cout << endl;
		}

		cout << endl << "Time sequence version: " << TotalSequentTime << endl;
		cout << "Time parallel version: " << TotalParallelTime << endl;

		cout << endl;

		// Проверка результатов
		int check = true;
		for (int i = 0; i < size; i++)
		{
			if (MasRadixPp[i] != MasRadixSeq[i])
			{
				cout << "MasRadixPp !=  MasRadixSeq" << endl;
				check = false;
				break;
			}
		}

		if (check)
			cout << "MasRadixPp ==  MasRadixSeq" << endl;

		delete[] MasRadixSeq;
	}

	delete[] MasRadixPp;
	delete[] CountWork;
	delete[] displs;
	delete[] rbuf;

	MPI_Finalize();

	return 0;
}