#ifndef FDM_H_
#define FDM_H_

#include "Grid2D_Quad.h"
#include "SLAU.h"
#include <tuple>
#include <functional>


/* Шаблон расчетной области - 5-и точечный крест => получается, что краевые узлы не входят в расчетную область */


struct BoundCond
{
    int num = 0; // Порядковый номер границы
    int typeBound = 0; // Тип КУ
    vector<int> idx; // Индексы на гарнице с учетом которых и будет производиться расчет 
    function<double(double)> func[2];
};

/* Эту структуру инициализируем руками для получения тестов. */
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
    BoundCond Bound[6];
};


class FDM
{

private:
    // false - не равномерная сетка
    // true - равномерная сетка
    bool IsravGrid = false;
    Grid2D_Quad Grid;
    Matrix matr;
    SLAUData slau;
    vector<double> uij; // Результирующий вектор решений 
    int N = 0; // размерность СЛАУ
    int m = 0; // Количество нулей до диагонали 

    /*  Параметры ДУ */
    DEquation deq;
    
    /* Массив внутренних точек Расчетной области */
    vector<int>InternalNode;

    vector<bool> CheckNode;
    
    /* Генерация Матрицы в которой уже учитывается все необходимое тоесть все КУ */
    bool GenerateMatrix();

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