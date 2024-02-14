#include "Grid_1D.h"
#include <iostream>
#include <cmath>

/* Private */

void Grid_1D::Load(const string &filename)
{
    fin.open(filename);

    if(!fin.is_open())
    {
        cout << "Файл можна????\n";
    }

    fin >> Nx;
    BaseGridX.resize((uint64_t)Nx);

    for(int32_t i = 0; i < Nx; i++)
        fin >> BaseGridX[(uint64_t)i];

    // Интервалов на 1 меньше чем опроных узлов 
    DivideParam.resize(uint64_t(Nx-1));
    for(int32_t i = 0; i < Nx-1; i++)
        fin >> DivideParam[(uint64_t)i].num >> DivideParam[(uint64_t)i].coef;
    
    cout << "File is load done\n";
    fin.close();
}

/* Public */

Grid_1D::Grid_1D(const string &filename)
{
    Load(filename);

    /* Рачет общего числа узлов */
    int32_t GlobalNx = 0;
    for(int32_t i = 0; i < Nx; i++)
        GlobalNx+=DivideParam[(uint64_t)i].num;

    GlobalNx++;
    Grid.resize((uint64_t)GlobalNx);
}

void Grid_1D::GenerateGrid()
{
     struct SettingForDivide
    {
        double step; // Шаг на отрезке
        double coef; // Коэффициент увеличения шага
        int num;     // Количество интервалов идем то num-1 и потом явно вставляем элемент
    };

    /* Расчитываем шаг для сетки  */
    /*
        @param
        int j - Номер элемента в массиве
        double left - левая грани отрезка
        double right - правая граница отрезка
        ret: SettingForDivide -  структура с вычесленными параметрами деления сетки
    */
    auto CalcSettingForDivide = [&](int j, double left, double right) -> SettingForDivide
    {   
        SettingForDivide res;
        int num = DivideParam[j].num;
        double coef = DivideParam[j].coef;

        if (coef > 1.0)
        {
            double coefStep = 1.0 + (coef * (std::pow(coef, num - 1) - 1.0)) / (coef - 1.0);

            res.step = (right - left) / coefStep;
        }
        else
        {
            res.step = (right - left) / num;
        }

        //  Убираем погрешность
        if (std::abs(res.step) < eps)
            res.step = 0.0;

        res.num = num;
        res.coef = coef;
        return res;
    };

    /* Генерация разбиения по X  с учетом разбиения */
    /*
        @param
        SettingForDivide &param - параметр разбиения
        double left - левая граница отрезка
        double right - правая граница отрезка
        double *Line - генерируемый массив
        int &idx - индекс в массиве на какую позицию ставить элемент
    */
    auto GenerateDivide = [](SettingForDivide &param, double left, double right, vector<double>& Line, int &idx) -> void
    {
        int num = param.num;
        double coef = param.coef;
        double step = param.step;

        Line[idx] = left;
        idx++;
        double ak = left;
        for (int k = 0; k < num - 1; k++)
        {
            ak = ak + step * std::pow(coef, k);
            Line[idx] = ak;
            idx++;
        }
        Line[idx] = right;
    };

    int idx = 0;
    for(int32_t j = 0; j < Nx-1; j++)
    {
        double left = BaseGridX[j];
        double right = BaseGridX[j+1];
        SettingForDivide param = CalcSettingForDivide(j, left, right);
        GenerateDivide(param, left, right, Grid, idx);
    }

}

void Grid_1D::DivideGrid(const int coef)
{
}

void Grid_1D::ReGenerateGrid()
{
}


void Grid_1D::PrintGrid() const
{
    cout << "Print 1D Grid: \n";
    for(int32_t i = 0; i < Grid.size(); i++)
    {
        cout << Grid[i] << " ";
    }
    cout << "\n";
}