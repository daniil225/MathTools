#ifndef GRID_H_
#define GRID_H_

#include <string>
#include <fstream>
#include "GridStatus.h"

using namespace std;


class GridI
{

private:


protected:
    double eps = 1e-7; //  машинный ноль 
    ofstream fout; // Файловый поток на запись
    ifstream fin; // Файловый поток на чтение 

public:

    GridI() = default;
    
    /* Загрузка базовой сетки из файла и по ней строится уже все */
    virtual GridStatus Load(const string& filename) noexcept;

    /* Генерация сетки */
    virtual GridStatus GenerateGrid() noexcept;

    /* Дробление сетки в заданное количество раз */
    virtual GridStatus DivideGrid(const int coef) noexcept;

    /* Перегенерация сетки при изменении ее параметров */
    virtual GridStatus ReGenerateGrid() noexcept;
    virtual ~GridI();

};

#endif