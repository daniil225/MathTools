#ifndef FDM_H_
#define FDM_H_

#include "Grid2D_Quad.h"
#include "SLAU.h"
#include <tuple>
#include <functional>


/* Шаблон расчетной области - 5-и точечный крест => получается, что краевые узлы не входят в расчетную область */

struct BoundCond
{
    /* 1 КУ */

    /* 3 КУ */
};

struct DEquation
{
    double lambda = 0;
    double gamma = 0;
    /* Правая часть ДУ */
    function<double(double, double)> f = [](double x, double y) -> double { return 0.0; }; // f = f(x,y)

    /* Истинное решение ДУ */
    bool IsInitTrue = false; // Имеется ли Истинное решение 
    function<double(double, double)> u_true;

    /* КУ */
    BoundCond Bound;
};


class FDM
{

private:
    bool IsravGrid = false;
    Grid2D_Quad Grid;
    Matrix matr;
    SLAUData slau;
    vector<double> uij; // Результирующий вектор решений 
    int N = 0; // размерность СЛАУ
    int m = 0; // Количество нулей до диагонали 

    /*  Параметры ДУ */
    DEquation deq;
    
    /* Генерация Матрицы в которой уже учитывается все необходимое тоесть все КУ */
    bool GenerateMatrix();

    /* Угловые узлы они в расчетную область не входят их всего 5 штук */
    pair<int, int> AngularPoint[5];

public: 

    FDM() = delete;
    FDM(const FDM&) = delete;

    FDM(const string &filename, DEquation deq_, bool IsravGrid_ = false);

    // Решение ДУ
    void Solve();

    // Расчет расстояния между истинным и расчитанным рещениями считаем как сумма Квадратов разности 
    // под Корнем в общем евклидовское расстояние 
    double Norma() const;


    inline void DivideGrid(const int coef) { Grid.DivideGrid(coef); Grid.ReGenerateGrid(); }
    

    FDM& operator=(const FDM&) = delete;

    ~FDM() = default;

};


#endif