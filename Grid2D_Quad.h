#ifndef GRID2D_QUAD_H_
#define GRID2D_QUAD_H_
#include "Grid.h"
#include "PointInfo.h"
#include <vector>

using namespace std;

// Сетка представляет из себя решулярное разбиение области при помощи 4-х угольников произвольной формы
/* Чтение и сетки из файла. Разбиение сетки.  */

/* Формат файла
* Nx Ny
* x1,1 y1,1 x2,1 y1,2 ... xn,1 y1,n
* x1,2 y2,1 y2,2 y2,2 ... x2,n y2,n
*
* L - целое число - количество подобластей и далее L наборов чисел по 5 штук
* 1 число - номер формул определяющих параметры ДУ в подобласти
* 2 число - первая вертикальная ломанная определеяет правую границу
* 3 число - вторая вертикальная ломанная определяет левую границу
  их номера определяются в соответствии с правилом для k узла их номера в массивах точек будут: k, (Nx + k), (2*Nx + k), ... , ( (Ny-1)*Nx + k )
* 4 число - опеределяет горизонтальную границу снизу
* 5 число - определяет горизонтальную границу сверху
* Дальше 3 сторики задающие разбиение сетки
* P - целое число - количество подоблостей граничных условий и далее P наборов чисел по 6 штук
* 1 число - номер формул определяющих параметры граничных условий
* 2 число - Тип краевого условия 0 - соответсвует не заданным краевым
* 3 число - первая вертикальная ломанная определеяет правую границу
* 4 число - вторая вертикальная ломанная определяет левую границу
* 5 число - опеределяет горизонтальную границу снизу
* 6 число - определяет горизонтальную границу сверху
* !!! Отличительная особенность одна из координатных линий должна быть фиксированной
* N1 coef1 N2 coef 2 ...
* ...
* ...
* Первый набор задает разбиение каждого отрезка по х
* Второй набора задает разбиение по y
*/

struct BaseGrid2D
{
    int CountOfDivision = 1; // Количество делений базовой настройки сетки; 1 - соответсвует тому что делений не было 2 тому что исходная область разделена в двое и.т.д
    const int SizeOfCalculationAreaElemet = 5;
    const int SizeOfBoundsAreaElement = 6;
    const int SizeOfDivideParam = 2;

    /* Точка */
    struct PointXY
    {
        double x;
        double y;
    };

    struct DivideParamS
    {
        int num;     // количество интервалов на которое нужно разделить отрезок
        double coef; // Коэффициент растяжения или сжатия
    };

    int Nx = 0; // количество узлов вдоль горизонатального направления
    int Ny = 0; // количество узлов вдоль вертикального направления
    int L = 0;  // Количество подоблостей
    int P = 0;  // Количество видов границ

    vector<vector<PointXY>> BaseGridXY; //  Базовая сетка в плосткости XY

    // | номер формул | граница по x(start) | граница по x(end) | | граница по y(start) | граница по y(end) |
    vector<vector<int>> CalculationArea; // Массив расчетных областей

    //| номер формул | | тип КУ | граница по x(start) | граница по x(end) | граница по y(start) | граница по y(end) |
    vector<vector<int>> BoundsArea;

    /*
        2 массива
        DivideParam[0] - массив для разбиения по оси x размер Nx-1
        DivideParam[1] - массив для разбиения по оси y размер Ny-1
    */
    vector<vector<DivideParamS>> DivideParam;

    /* Флаг того что структура инициализированная и готова к использованию */
    // false - Структура не инициализированная
    // true - Структура проинициализированна и готова к использованию
    bool isReadyToUse = false;
};

struct Point
{
    // Информация о границе конкретный набор формул + информация об узле фиктивный или нет + информация о том граничный ли узел или нет
    // + информация о конечном элементе + тип КУ + точка в пространтве
    Info info;
    double x = 0.0;
    double y = 0.0;
};

/* Класс сетки */
class Grid2D_Quad : public GridI
{

private:
    BaseGrid2D baseGrid;

    int Dim = 0; // Размерность сетки и СЛАУ в то же время
    int Nx = 0;  // Сумарное количество узлов по оси Х
    int Ny = 0;  // Сумарное количество узлов по оси У все величины получаются после генерации сетки

    vector<Point> Grid; // Массив точек получающийся при генерации конечных элементов

protected:
public:
    Grid2D_Quad() = default;

    // Инициализация класса при помощи файла - формат описан выше
    /*
        @param
        const string &filename - Текстовый файл с разметкой
        @return void
        Создание объекта
    */
    explicit Grid2D_Quad(const string &filename);

    // Инициализация "в ручную"  при помощи задания структуры
    /*
        @param
       const BaseGrid2D &baseGrid_ - Структура с базовой разметкой области
        @return void
        Создание объекта
    */
    explicit Grid2D_Quad(const BaseGrid2D &baseGrid_);

    // Копирование класса отключено
    Grid2D_Quad(const Grid2D_Quad &) = delete;
    Grid2D_Quad &operator=(const Grid2D_Quad &) = delete;

    /* Методы интерфейса обязательны к реализации */

    /* Загрузка базовой сетки из файла и по ней строится уже все */
    /*
        @param
        const string &filename - Текстовый файл с разметкой
        @return GridStatus
    */
    GridStatus Load(const string &filename) noexcept;

    /* Генерация сетки */
    /*
        @param:
        void
        @return: GridStatus
    */
    GridStatus GenerateGrid() noexcept;

    /* Дробление сетки в заданное количество раз */
    /*
        @param:
        const int coef - Коэффициент дробления
        @return: GridStatus
    */
    GridStatus DivideGrid(const int coef) noexcept;

    /* Перегенерация сетки при изменении ее параметров */
    /*
        @param:
        void
        @return: GridStatus
    */
    GridStatus ReGenerateGrid() noexcept;

    /* Гетеры и сеттеры */
    /*
        @param:
        void
        @return: BaseGrid2D
    */
    BaseGrid2D GetBaseGrid() noexcept;

    /* Установка параметров базовой сетки */
    /*
        @param:
        const BaseGrid2D& baseGrid_ - Базовая сетка области
        @return: GridStatus
    */
    GridStatus SetBaseGrid(const BaseGrid2D &baseGrid_) noexcept;

    ~Grid2D_Quad() = default;
};

#endif