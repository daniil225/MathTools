#include "SLAU.h"
#include <algorithm>
#include <cmath>
#include <iostream>


/* Matrix */
void Matrix::SaveMemory()
{

    di.resize(N);
    ig.resize(Nd);

    ggu.resize(Nd);
    ggl.resize(Nd);

    jOffset.resize(2 * Nd);

    //  Нижний треугольник

    for (int i = 0; i < Nd; i++)
    {
        ggl[i].resize(N - 1); // Сразу после главной диагонали

        // Верхний треугольник
        ggu[i].resize(N - 1); // Сразу после главной диагонали
    }

    // Составим массив ig
    ig[0] = 1;
    ig[1] = m + 2;
    // Составляем массив оффетов
    jOffset = {-2 - m, -1, 1, 2 + m};
}

std::vector<double> Matrix::MulMatrVec(const std::vector<double> &x)
{
    std::vector<double> res(N);

    for (size_t i = 0; i < (size_t)N; i++)
        res[i] = x[i] * di[i];

    for (size_t i = 0; i < (size_t)Nd; i++)
    {
        for (size_t j = 0; j < (size_t)N - ig[i]; j++)
        {
            size_t ir = j + ig[i];
            res[ir] += ggl[i][j] * x[j];
            res[j] += ggu[i][j] * x[ir];
        }
    }

    return res;
}

bool Matrix::add(int i, int j, double val)
{
bool res = true;
    std::vector<int>::iterator idx; // Итератор для поиска в списке ig
    int k; // Номер строки в массивах ggl, ggu

    if( i > j)
    {
        int ij = i - j;
        idx = std::find(ig.begin(), ig.end(), ij);
        if (idx != ig.end())
        {
            k = idx - ig.begin();
            ggl[k][j] += val;
        }
        else
        {
            res = false;
        }
    }
    else if(i < j)
    {
        int ji = j - i;
        idx = std::find(ig.begin(), ig.end(), ji);
        if (idx != ig.end())
        {
            k = idx - ig.begin();
            ggu[k][i] += val;
        }
        else
        {
            res = false;
        }

    }
    else
    {
        di[i] += val;
    }

    return res;
}

bool Matrix::insert(int i, int j, double val)
{
    bool res = true;
    std::vector<int>::iterator idx; // Итератор для поиска в списке ig
    int k; // Номер строки в массивах ggl, ggu

    if( i > j)
    {
        int ij = i - j;
        idx = std::find(ig.begin(), ig.end(), ij);
        if (idx != ig.end())
        {
            k = idx - ig.begin();
            ggl[k][j] = val;
        }
        else
        {
            res = false;
        }
    }
    else if(i < j)
    {
        int ji = j - i;
        idx = std::find(ig.begin(), ig.end(), ji);
        if (idx != ig.end())
        {
            k = idx - ig.begin();
            ggu[k][i] = val;
        }
        else
        {
            res = false;
        }

    }
    else
    {
        di[i] = val;
    }

    return res;
}


double Matrix::GetElement(int i, int j)
{
    double res = 0;
    std::vector<int>::iterator idx; // Итератор для поиска в списке ig
    int k; // Номер строки в массивах ggl, ggu
    
    std::vector<int> j_ = jOffset;
    
    if(i == j)
    {
        res = di[i];
    }
    else if( i > j)
    {
        int ij = i - j;
        idx = std::find(ig.begin(), ig.end(), ij);
        if (idx != ig.end())
        {
            k = idx - ig.begin();
            res = ggl[k][j];
        }
    }
    else if(i < j)
    {
        int ji = j - i;
        idx = std::find(ig.begin(), ig.end(), ji);
        if (idx != ig.end())
        {
            k = idx - ig.begin();
            res = ggu[k][i];
        }
    }

    return res;
}

void Matrix::PrintDenseFormatMatrix()
{   
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            std::cout << GetElement(i,j) << " ";
        }
        std::cout << "\n";
    }
}

/***********************************************************************************************/
double SumRow(Matrix &matr, std::vector<double> &xk, int i)
{
    double res = matr.di[i] * xk[i];
    std::vector<int>::iterator idx; // Итератор для поиска в списке ig

    // Определеяем псевдонимы (ссылки)
    std::vector<int> &ig = matr.ig;
    std::vector<std::vector<double>> &ggl = matr.ggl;
    std::vector<std::vector<double>> &ggu = matr.ggu;

    int k; // Номер строки в массивах ggl, ggu

    std::vector<int> j_ = matr.jOffset;
    for (int k = 0; k < matr.Nd * 2; k++)
        j_[k] += i;

    for (int j : j_)
    {
        if (j >= 0)
        {
            if (i > j)
            {
                int ij = i - j;
                idx = std::find(ig.begin(), ig.end(), ij);
                if (idx != ig.end())
                {
                    k = idx - ig.begin();
                    res += ggl[k][j] * xk[j];
                }
            }
            else if (i < j)
            {
                int ji = j - i;
                idx = std::find(ig.begin(), ig.end(), ji);
                if (idx != ig.end())
                {
                    k = idx - ig.begin();
                    res += ggu[k][i] * xk[j];
                }
            }
        }
    }

    return res;
}


std::vector<double> subVec(std::vector<double>& x, const std::vector<double>& y)
{
	std::vector<double> res(x.size());
	for (size_t i = 0; i < x.size(); i++)
		res[i] = x[i] - y[i];

	return res;
}


double NormVec(const std::vector<double>& vec)
{
	double res = 0.0;
	for (size_t i = 0; i < vec.size(); i++)
		res += vec[i] * vec[i];

	return std::sqrt(res);
}

void Jakobi(Matrix& matr, std::vector<double>& xk1,SLAUData &data, const double w)
{
	std::vector<double> xk(matr.N, 0);
	int MAX_ITER = data.MAX_ITER;
	std::vector<double>& f = data.f;
	double eps = data.eps; // Заданная точность 
	double F_norm = NormVec(data.f); // Норма вектора правой части 
	double NonRepan; // Относительная невязка 
	// Итерации идут пока не будет достигнута максимальное число 

	int n = matr.N;
	std::vector<double>& di = matr.di;
	for (int k = 0; k < MAX_ITER; k++)
	{
		for (int i = 0; i < n; i++)
		{
			xk1[i] = xk[i] + w*(f[i] - SumRow(matr, xk, i)) / matr.di[i];
		}
		xk = xk1;
		// Вычисляем относительную невязку для выхода из цикла по ней 
		NonRepan = NormVec(subVec(f, matr.MulMatrVec(xk1))) / F_norm;
		std::cout << "Iteration = " << k + 1 << "  Non-repan = " << NonRepan << "\n";

		if (NonRepan < eps) break;
	}
}

void Zeidel(Matrix& matr, std::vector<double> &xk1, SLAUData &data, const double w)
{
	int n = matr.N;
	int MAX_ITER = data.MAX_ITER;
	std::vector<double>& f = data.f;
	double NonRepan;

	double eps = data.eps; // Заданная точность 
	double F_norm = NormVec(data.f); // Норма вектора правой части
	std::vector<double>& di = matr.di;
	for (int k = 0; k < MAX_ITER; k++)
	{
		for (int i = 0; i < n; i++)
		{
			xk1[i] = xk1[i] + w*(f[i] - SumRow(matr, xk1, i)) / matr.di[i];
		}
		NonRepan = NormVec(subVec(f, matr.MulMatrVec(xk1))) / F_norm;
		std::cout << "Iteration = " << k+1 << "  Non-repan = " << NonRepan << "\n";
		if (NonRepan < eps) break; // Выход по неявязке 
	}
}