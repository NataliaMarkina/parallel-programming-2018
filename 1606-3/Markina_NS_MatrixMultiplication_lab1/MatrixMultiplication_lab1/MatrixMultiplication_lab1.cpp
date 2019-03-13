// MatrixMultiplication_lab1.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <time.h>
#include <ctime>
#include <iomanip>

using namespace std;

/* Номера матриц (для srand, чтобы м-цы генерировались разными) */
#define FIRST_CMT 1
#define SECOND_CMT 2

/* Максимальный размер плотной м-цы */
#define MAX_SIZE_DENSE_MATRIX 13000


struct CRSMatrix
{
	int Size; // Размер матрицы Size x Size
	int NonzeroSize; // Число ненулевых эл-тов в матрице

	double* Values; // Массив ненулевых эл-тов м-цы (размер NonzeroSize)
	int* Columns; // Массив номеров столбцов, в которых находятся элементы массива Value (размер NonzeroSize)
	int* RowIndex; // Массив индексов строк (размер Size+1)
};

// Создать CRS матрицу
void CreateCRSMatrix(int Size, int NonzeroSize, CRSMatrix& CMT)
{
	CMT.Size = Size;
	CMT.NonzeroSize = NonzeroSize;

	CMT.Values = new double[NonzeroSize];
	CMT.Columns = new int[NonzeroSize];
	CMT.RowIndex = new int[Size + 1];
}

// Удалить CRS матрицу
void DeleteCRSMatrix(CRSMatrix& CMT)
{
	delete[] CMT.Values;
	delete[] CMT.Columns;
	delete[] CMT.RowIndex;
}

// Заполнение разреженной CRS матрицы
// (3 массива, индексация с 0)
// В каждой строке NumElInRow ненулевых элементов
void InitCRSMatr(int NumOfMatr, int Size, int NumElInRow, CRSMatrix& CMT)
{
	int NonzeroSize, tmp;
	bool flag = false;

	NonzeroSize = NumElInRow * Size;
	srand(NumOfMatr);

	CreateCRSMatrix(Size, NonzeroSize, CMT);

	for (int i = 0; i < Size; i++)
	{
		for (int j = 0; j < NumElInRow; j++)
		{
			do
			{
				flag = false;
				CMT.Columns[i * NumElInRow + j] = rand() % Size;
				for (int k = 0; k < j; k++)
					if (CMT.Columns[i * NumElInRow + j] == CMT.Columns[i * NumElInRow + k])
						flag = true;
			} 
			while (flag == true);
		}

		for (int j = 0; j < NumElInRow - 1; j++)
			for (int k = 0; k < NumElInRow - 1; k++)
				if (CMT.Columns[i * NumElInRow + k] > CMT.Columns[i * NumElInRow + k + 1])
				{
					tmp = CMT.Columns[i * NumElInRow + k];
					CMT.Columns[i * NumElInRow + k] = CMT.Columns[i * NumElInRow + k + 1];
					CMT.Columns[i * NumElInRow + k + 1] = tmp;
				}
	}

	for (int i = 0; i < NonzeroSize; i++)
	{
		CMT.Values[i] = rand() % 100 - 50;

		if (CMT.Values[i] == 0)
			CMT.Values[i] += 1;
	}

	CMT.RowIndex[0] = 0;
	for (int i = 1; i < Size + 1; i++)
		CMT.RowIndex[i] = CMT.RowIndex[i - 1] + NumElInRow;
}

// Транспонирование матрицы
CRSMatrix Transposing(CRSMatrix CMT)
{
	CRSMatrix TCMT;
	int TSize = CMT.Size;
	int TNonzeroSize = CMT.NonzeroSize;

	CreateCRSMatrix(TSize, TNonzeroSize, TCMT);

	memset(TCMT.RowIndex, 0, sizeof(int) * (TSize + 1));

	for (int i = 0; i < TNonzeroSize; i++)
		TCMT.RowIndex[CMT.Columns[i] + 1]++;

	int offset = 0;
	int tmp;
	for (int i = 1; i <= TSize; i++)
	{
		tmp = TCMT.RowIndex[i];
		TCMT.RowIndex[i] = offset;
		offset += tmp;
	}

	for (int i = 0; i < TSize; i++)
	{
		int start = CMT.RowIndex[i];
		int end = CMT.RowIndex[i + 1];
		int Column = i;
		int Value;
		int CMTcolumn, TCMTindex;

		for (int j = start; j < end; j++)
		{
			Value = CMT.Values[j];
			CMTcolumn = CMT.Columns[j];
			TCMTindex = TCMT.RowIndex[CMTcolumn + 1];
			TCMT.Values[TCMTindex] = Value;
			TCMT.Columns[TCMTindex] = Column;
			TCMT.RowIndex[CMTcolumn + 1]++;
		}
	}

	return TCMT;
}

// Оптимизированная последовательная версия умножения матриц
CRSMatrix seq_MultiplicationMatrCRS(CRSMatrix CMT_1, CRSMatrix TCMT_2)
{
	CRSMatrix CMT_Rez;
	vector<double> Value;
	vector<int> Column, Row;

	int Size = CMT_1.Size;
	int *tmp = new int[Size];
	int NzSize = 0;

	Row.push_back(0);

	for (int i = 0; i < Size; i++)
	{
		memset(tmp, -1, Size * sizeof(int));

		int start = CMT_1.RowIndex[i];
		int end = CMT_1.RowIndex[i + 1];

		for (int j = start; j < end; j++)
		{
			int col = CMT_1.Columns[j];
			tmp[col] = j;
		}

		for (int j = 0; j < Size; j++)
		{
			double ScalarMult = 0;

			int tstart = TCMT_2.RowIndex[j];
			int tend = TCMT_2.RowIndex[j + 1];

			for (int k = tstart; k < tend; k++)
			{
				int bcol = TCMT_2.Columns[k];
				int aind = tmp[bcol]; 
				if (aind != -1)
					ScalarMult += CMT_1.Values[aind] * TCMT_2.Values[k];
			}

			if (ScalarMult != 0)
			{
				Column.push_back(j);
				Value.push_back(ScalarMult);
				NzSize++;
			}
		}
		Row.push_back(NzSize);
	}

	CreateCRSMatrix(Size, Value.size(), CMT_Rez);

	for (int j = 0; j < Value.size(); j++)
	{
		CMT_Rez.Values[j] = Value[j];
		CMT_Rez.Columns[j] = Column[j];
	}

	for (int j = 0; j < Size + 1; j++)
		CMT_Rez.RowIndex[j] = Row[j];

	delete[] tmp;
	return CMT_Rez;
}

// Отобразить массив
void Show_arr(int* arr, int size_arr)
{
	if (arr == NULL || size_arr < 1)
		return;

	for (int i = 0; i < size_arr; i++)
		cout << arr[i] << " ";

	cout << endl;
}

// Отобразить массив
void Show_arr(double* arr, int size_arr)
{
	if (arr == NULL || size_arr < 1)
		return;

	for (int i = 0; i < size_arr; i++)
		cout << arr[i] << " ";

	cout << endl;
}

// Отобразить матрицу
void Show_matr(double** matr, int size)
{
	for (int i = 0; i < size; i++)
		Show_arr(matr[i], size);
}

// Отобразить CRS матрицу
void ShowCRSMatr(CRSMatrix CMT)
{
	cout << endl;
	Show_arr(CMT.Values, CMT.NonzeroSize);
	Show_arr(CMT.RowIndex, CMT.Size + 1);
	Show_arr(CMT.Columns, CMT.NonzeroSize);
	cout << endl;
}



// Сравнение плотных м-ц
bool DenseMatrixsAreEqual(double** M_1, double** M_2, int Size)
{
	for (int i = 0; i < Size; i++)
		for (int j = 0; j < Size; j++)
			if (M_1[i][j] != M_2[i][j])
				return false;

	return true;
}

// Создать и проинициализировать матрицу
double** Create_and_init_matr(int size)
{
	if (size < 1)
		return NULL;

	double** matr;

	matr = new double*[size];
	for (int i = 0; i < size; i++)
		matr[i] = new double[size];

	srand((unsigned)time(NULL));

	if (size < 3)
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			{
				cout << "Enter element at the [" << i << ", " << j << "] position:" << endl;
				cin >> matr[i][j];
				cout << endl;
			}
	else
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				matr[i][j] = rand() % 100 - 50;

	return matr;
}

// Умножение плотных м-ц
double** seq_MultiplicationDenseMatrix(double** M_1, double** M_2, int Size)
{
	double** M_Rez = Create_and_init_matr(Size);

	for (int i = 0; i < Size; i++)
		for (int j = 0; j < Size; j++)
		{
			M_Rez[i][j] = 0;
			for (int k = 0; k < Size; k++)
				M_Rez[i][j] += M_1[i][k] * M_2[k][j];
		}

	return M_Rez;
}

// Преобразовать CRS разреженную матрицу к обычному виду (не работает для больших значений размера матриц)
double** ConvertCRSToSimpleMatrix(CRSMatrix CMT)
{
	double** MX = Create_and_init_matr(CMT.Size);
	int k = 0, l = 0;
	int SizeElOnRow = 0;

	for (int i = 0; i < CMT.Size; i++)
		for (int j = 0; j < CMT.Size; j++)
			MX[i][j] = 0;


	for (int i = 0; i < CMT.Size; i++)
	{
		SizeElOnRow = CMT.RowIndex[i + 1] - CMT.RowIndex[i];
		for (int j = 0; j < CMT.Size; j++)
			if (CMT.Columns[k] == j)
			{
				MX[i][j] = CMT.Values[k];
				k++;
				l++;
				if (l == SizeElOnRow)
					break;
			}

		l = 0;
	}

	return MX;
}

// Сравнить результат умножения CRS матриц с обычным умножением плотных матриц (Это работает при условии, если Size < 13000)
bool CheckSeqAndSimpleRezults(CRSMatrix CMT_1, CRSMatrix CMT_2, CRSMatrix CMT_Rez)
{
	if (CMT_1.Size > MAX_SIZE_DENSE_MATRIX)
	{
		cout << "Size of matrix must be less than " << MAX_SIZE_DENSE_MATRIX << endl;
		return false;
	}

	if (CMT_1.Size != CMT_2.Size)
	{
		cout << "Sizes of CMT_1 and CMT_2 differ " << endl;
		return false;
	}

	double** CMT_1_simple = ConvertCRSToSimpleMatrix(CMT_1);
	double** CMT_2_simple = ConvertCRSToSimpleMatrix(CMT_2);
	double** CMT_Rez_simple = ConvertCRSToSimpleMatrix(CMT_Rez);
	double** Dense_Rez;

	Dense_Rez = seq_MultiplicationDenseMatrix(CMT_1_simple, CMT_2_simple, CMT_1.Size);

	if (CMT_1.Size < 10)
	{
		Show_matr(CMT_1_simple, CMT_1.Size);
		cout << endl;
		Show_matr(CMT_2_simple, CMT_1.Size);
		cout << endl;
		Show_matr(Dense_Rez, CMT_1.Size);
		cout << endl;
		Show_matr(CMT_Rez_simple, CMT_1.Size);
	}
	

	if (DenseMatrixsAreEqual(Dense_Rez, CMT_Rez_simple, CMT_1.Size))
	{
		cout << endl << "Multiplication of CRS Matrix is correct " << endl << endl;
		delete[] CMT_1_simple;
		delete[] CMT_2_simple;
		delete[] CMT_Rez_simple;
		return true;
	}
	else
	{
		cout << endl << "Multiplication of CRS Matrix is not correct " << endl << endl;
		delete[] CMT_1_simple;
		delete[] CMT_2_simple;
		delete[] CMT_Rez_simple;
		return false;
	}
}


void main(int argc, char* argv[])
{
	CRSMatrix CMT_1, CMT_2, TCMT_2, CMT_Rez;
	int Size;
	int NonzeroSize;
	double seq_start_time = 0.0;
	double seq_end_time = 0.0;
	double seq_work_time = 0.0;

	cout << "Enter size of matrix:   ";
	cin >> Size;
	cout << endl;

	cout << "Enter the number of nonzero elements in the string:   ";
	cin >> NonzeroSize;
	cout << endl << endl;

	if (NonzeroSize > Size || Size < 1 || NonzeroSize < 0)
	{
		cout << "Error" << endl;
		return;
	}

	/* Начало инициализации данных */

	InitCRSMatr(FIRST_CMT, Size, NonzeroSize, CMT_1);
	InitCRSMatr(SECOND_CMT, Size, NonzeroSize, CMT_2);

	TCMT_2 = Transposing(CMT_2);

	/* Конец инициализации данных */

	/* Начало последовательной версии */

	seq_start_time = clock();
	CMT_Rez = seq_MultiplicationMatrCRS(CMT_1, TCMT_2);
	seq_end_time = clock();

	seq_work_time = (double)(seq_end_time - seq_start_time) / CLOCKS_PER_SEC;

	if (Size <= 1000)
		CheckSeqAndSimpleRezults(CMT_1, CMT_2, CMT_Rez);

	cout << "Sequence version of multiplication matrix is worked: " << seq_work_time << endl;

	/* Конец последовательной версии */

	DeleteCRSMatrix(CMT_1);
	DeleteCRSMatrix(CMT_2);
	DeleteCRSMatrix(CMT_Rez);

	system("pause");

}

