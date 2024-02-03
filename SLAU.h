#ifndef SLAU_H_
#define SLAU_H_

#include <vector>

using namespace std;

struct Matrix
{
    int N = 0; // Размерность матрицы 
    int m = 0; // Количество не нулевых 
    int Nd = 2;
    vector<vector<double>> ggu;
    vector<vector<double>> ggl;
    vector<int> ig;
    vector<double> di;

    vector<int> jOffset;

    void SaveMemory();
	//  умножение матрицы на вектор  
	std::vector<double> MulMatrVec(const std::vector<double>& x);

    /* Вставка элемента на нужно место */
    /* 
        true - успешно вставили элемент 
        false - промазали с индексами 
     */
    bool insert(int i, int j, double val);

    double GetElement(int i, int j);

    /* Добавить значение в заданную ячейку */
    bool add(int i, int j, double val);

    /* Распечатка матрицы в плотном вормате */
    void PrintDenseFormatMatrix();
};



// Алгоритм работает за константное время так, как нужные индексы однозначно расчитываются 
// и все операции не зависят от размера входных данных 
double SumRow(Matrix& matr, std::vector<double>& xk, int i);


std::vector<double> subVec(std::vector<double>& x, const std::vector<double>& y);

// Норма вектора в евклидовом пространстве 
double NormVec(const std::vector<double>& vec);



/* Структура для генерации СЛАУ */
struct SLAUData
{
	int MAX_ITER = 10e6; // Базовое количество максимального числа итераций 
	double eps = 1e-10; // Базовое значение невязки 
	std::vector<double> f;
};


// data - Содержит в себе вектор правой части, макисимальное количество итераций, заданную точность 
// xk1 - начальное приближение и результирующий вектор наше решение кароче 
// matr - матрица СЛАУ 
void Jakobi(Matrix& matr, std::vector<double>& xk1,SLAUData &data, const double w = 1.0);

// data - Содержит в себе вектор правой части, макисимальное количество итераций, заданную точность 
// xk1 - начальное приближение и результирующий вектор наше решение кароче 
// matr - матрица СЛАУ 
void Zeidel(Matrix& matr, std::vector<double> &xk1, SLAUData &data, const double w = 1.0);


#endif