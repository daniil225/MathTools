#include "Grid2D_Quad.h"
#include "PointInfo.h"
#include <iostream>
#include <cmath>
#include <iomanip>

/* Private section */
int Grid2D_Quad::Getlevel(int i, int axis) const noexcept
{
    int res = 0;
    for (int k = 0; k < i; k++)
        res += baseGrid.DivideParam[axis][k].num;
    return res;
};

void Grid2D_Quad::GetTotalNumberOfNodes() noexcept
{
    for (int i = 0; i < baseGrid.Nx - 1; i++)
        GlobalNx += baseGrid.DivideParam[0][i].num;

    for (int i = 0; i < baseGrid.Ny - 1; i++)
        GlobalNy += baseGrid.DivideParam[1][i].num;

    GlobalNx++;
    GlobalNy++;
}

void Grid2D_Quad::GenerateBaseGrid(GridStatus &status) noexcept
{
    if (status.GetState() != State::OK)
        return;

    struct SettingForDivide
    {
        double step; // Шаг на отрезке
        double coef; // Коэффициент увеличения шага
        int num;     // Количество интервалов идем то num-1 и потом явно вставляем элемент
    };
    /* Расчитываем шаг для сетки  */
    /*
        @param
        int i - Номер массива от 0 до 2
        int j - Номер элемента в массиве
        double left - левая грани отрезка
        double right - правая граница отрезка
        ret: SettingForDivide -  структура с вычесленными параметрами деления сетки
    */
    auto CalcSettingForDivide = [&](int i, int j, double left, double right) -> SettingForDivide
    {
        SettingForDivide res;
        int num = baseGrid.DivideParam[i][j].num;
        double coef = baseGrid.DivideParam[i][j].coef;

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
    /* Генерация разбиения по X или Y( гороизонтальная линия или вертикальная ) с учетом разбиения */
    /*
        @param
        SettingForDivide &param - параметр разбиения
        double left - левая граница отрезка
        double right - правая граница отрезка
        double *Line - генерируемый массив
        int &idx - индекс в массиве на какую позицию ставить элемент
    */
    auto GenerateDivide = [](SettingForDivide &param, double left, double right, double *Line, int &idx) -> void
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

    try
    {
        Dim = GlobalNx * GlobalNy;
        Grid.resize(Dim);
        // Псевдоним для быстрого обращения
        int Nx = baseGrid.Nx;
        int Ny = baseGrid.Ny;
        vector<vector<BaseGrid2D::PointXY>> &BaseGridXY = baseGrid.BaseGridXY;

        double *LineX = new double[GlobalNx]; // Массив элементов в строке по Х
        double *LineY = new double[GlobalNy]; // Массив элементов в строке по Y
        // Сгенерируем одну плоскость xy

        /* Разбиение по x и y */
        /* Расстановка элементов основных линий с учетом их расположения в сетке ( Опорная сетка ) */

        /* Расстановка элементов по x ( Опорных ) */

        /* Нужно значения х расставить в соответсвующие строки они соответсвуют разбиению по z */

        for (int i = 0; i < Ny; i++)
        {
            int idx = 0;

            for (int j = 0; j < Nx - 1; j++)
            {
                double left = BaseGridXY[i][j].x;
                double right = BaseGridXY[i][j + 1].x;
                SettingForDivide param = CalcSettingForDivide(0, j, left, right);
                GenerateDivide(param, left, right, LineX, idx);
            }
            /* Заносим соответствующие значения x на свои позиции */

            int startIdx = Getlevel(i, 1) * GlobalNx;
            int endIdx = startIdx + GlobalNx;
            for (int k = startIdx, kk = 0; k < endIdx; k++, kk++)
                Grid[k].x = LineX[kk];
        }

        /* Расстановка элементов по z ( Опорных ) */
        for (int i = 0; i < Nx; i++)
        {
            int idx = 0;
            for (int j = 0; j < Ny - 1; j++)
            {
                double left = BaseGridXY[j][i].y;
                double right = BaseGridXY[j + 1][i].y;
                SettingForDivide param = CalcSettingForDivide(1, j, left, right);
                GenerateDivide(param, left, right, LineY, idx);
            }

            /* Процедура расстановки узлов в глобальный массив */
            int startIdx = Getlevel(i, 0); // Стартовый индекс для прохода по массиву
            for (int k = 0; k < GlobalNy; k++)
            {
                // Скачки будут ровно на величину GlobalNx
                Grid[startIdx].y = LineY[k];
                startIdx += GlobalNx;
            }
        }

        /****************************************************/

        /* Генерация вспомогательных вертикальных линий  */
        /*
            Кратко Алгоритм:
            Работаем с осью Y соответственно индексация будет происходить по этой оси
            в цикле идем по всем столбцам массива сетки

            Нужно получить левую и правую границу на каждом интервале
            Сформировать массив отрезков по  данной координате
            Занести полученный массив в Глобальную сетку
        */
        /* Цикл по всем горизонтальным линиям */
        for (int i = 0; i < GlobalNx; i++)
        {
            int idx = 0;
            /* Цикл по интеравалам оси Z */
            for (int j = 0; j < Ny - 1; j++)
            {
                int startIdx = i + GlobalNx * Getlevel(j, 1);
                int endIdx = i + GlobalNx * Getlevel(j + 1, 1);
                double left = Grid[startIdx].x; // Левая граница по х
                double right = Grid[endIdx].x;  // Правая граница по х

                // Разбиение интервала подчиняется разбиению по оси y
                SettingForDivide param = CalcSettingForDivide(1, j, left, right);
                GenerateDivide(param, left, right, LineY, idx);
            }
            /* Занесение результата в Итоговый массив */

            int startIdx = i; // Стартовая позиция
            for (int k = 0; k < GlobalNy; k++)
            {
                Grid[startIdx].x = LineY[k];
                startIdx += GlobalNx;
            }
        }
        /****************************************************/

        /* Генерация вспомогательных горизонтальных линий  */
        /* Цикл по всем горизонтальным линиям */
        for (int i = 0; i < GlobalNy; i++)
        {
            int idx = 0;
            for (int j = 0; j < Nx - 1; j++)
            {
                int startIdx = Getlevel(j, 0) + i * GlobalNx;
                int endIdx = Getlevel(j + 1, 0) + i * GlobalNx;
                double left = Grid[startIdx].y;
                double right = Grid[endIdx].y;
                // Разбиение интервала подчиняется разбиению по оси x
                SettingForDivide param = CalcSettingForDivide(0, j, left, right);
                GenerateDivide(param, left, right, LineX, idx);
            }
            /* Занесение результатов в Глобальную сетку */
            int startIdx = i * GlobalNx;
            for (int k = 0; k < GlobalNx; k++)
            {
                Grid[startIdx].y = LineX[k];
                startIdx++;
            }
        }
        /****************************************************/

        /* Очистка памяти */
        delete[] LineX;
        delete[] LineY;
    }
    catch (std::bad_alloc &e)
    {
        status.SetStatus(State::MEMORY_ALLOC_ERROR, "Ошибка при выделении памяти в Grid2D_Quad::GenerateBaseGrid(GridStatus &status)\n");
        return;
    }
    catch (const std::exception &e)
    {
        status.SetStatus(State::UNKNOWN_ERROR, e.what());
        return;
    }
}

void Grid2D_Quad::DivisionIntoSubAreas(GridStatus &status) noexcept
{
    if (status.GetState() != State::OK)
        return;

    struct BoundArea
    {
        int AreaNum = -1;    // номер подобласти (определяет набор формул отвечающих параметрам ДУ)
        int PlaneXZSize = 0; // Размер массива PlaneXZ
        int *PlaneXZ;        // Массив точек многоугольника (номера)
    };
    int *BoundGrid = nullptr;
    try
    {
        /* Размер массива берем с запасом его размер равен Nx*Ny
            контур номеров точек многоугольника соттветствующих каким то координатам
        */
        BoundGrid = new int[GlobalNx * GlobalNy];

        /* Получает подобласть в соответсвии с ее номером */
        /*
        @param
         int i - Номер подобласти (Порядковый) в массиве определяется порядок
         ret BoundArea - сформированный массив подобласти
        */
        auto GetBound = [&](int i) -> BoundArea
        {
            BoundArea Bound;

            Bound.AreaNum = baseGrid.CalculationArea[i][0]; // Выставили номер подобласти

            /* Вычислим все номера принадлежащие данной области */
            /* Инедексы границ многоугольника в глобальной нумерации */
            /* XY */

            int leftStartX = Getlevel(baseGrid.CalculationArea[i][1], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[i][3], 1);
            int leftEndX = Getlevel(baseGrid.CalculationArea[i][1], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[i][4], 1);
            int rightStartX = Getlevel(baseGrid.CalculationArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[i][3], 1);
            int rightEndX = Getlevel(baseGrid.CalculationArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[i][4], 1);

            int nX = (rightEndX - rightStartX) / GlobalNx;
            int nY = (rightStartX - leftStartX);

            int idxBoundGrid = 0;
            for (int i = 0; i <= nX; i++)
            {
                int Idx = leftStartX + i * GlobalNx;
                for (int j = 0; j <= nY; j++)
                {
                    BoundGrid[idxBoundGrid] = Idx;
                    idxBoundGrid++;
                    Idx++;
                }
            }
            Bound.PlaneXZSize = idxBoundGrid;
            Bound.PlaneXZ = BoundGrid;

            return Bound;
        };
        /* Бинарный поиск по массиву */
        /*
            @param
            const BoundArea& Bound - Массив границ
            int numPointGlobal - Значение для поиска
            ret bool: true - элемент в массиве есть false в противном случае

        */
        auto BinarySerch = [](const BoundArea &Bound, int numPointGlobal) -> bool
        {
            int *arr = Bound.PlaneXZ;
            int left = 0;
            int right = Bound.PlaneXZSize - 1;
            int midd = 0;
            while (1)
            {
                midd = (left + right) / 2;

                if (numPointGlobal < arr[midd])      // если искомое меньше значения в ячейке
                    right = midd - 1;                // смещаем правую границу поиска
                else if (numPointGlobal > arr[midd]) // если искомое больше значения в ячейке
                    left = midd + 1;                 // смещаем левую границу поиска
                else                                 // иначе (значения равны)
                    break;                           // функция возвращает индекс ячейки

                // если границы сомкнулись
                if (left > right)
                {
                    midd = -1;
                    break;
                }
            }
            if (midd == -1)
                return false;
            else
                return true;
        };

        /* Проверяет принадлежит ли точка Заданной области */
        /*
            @param
            const BoundArea &Bound - Заданная граница
            Point &point - Точка для проверки
            int numPointGlobal - индекс точки в глобальной нумерации
            ret bool: true - Точка принадлежит заданной области, false в противном случае
        */
        auto IsInArea = [&](const BoundArea &Bound, Point &point, int numPointGlobal) -> bool
        {
            auto le = [&](double x1, double x2) -> bool
            {
                if (x1 < x2)
                    return true;
                if (std::abs(x1 - x2) < eps)
                    return true;
                return false;
            };

            auto ge = [&](double x1, double x2) -> bool
            {
                return le(x2, x1);
            };

            bool arg1 = BinarySerch(Bound, numPointGlobal);
            return arg1;
        };

        /* В цикле по всем областям */
        for (int i = 0; i < baseGrid.L; i++)
        {
            BoundArea Bound = GetBound(i);

            /* По той части области где распологаются элементы (включая фиктивный)*/
            /* Для этого получим минимальный и максимальный значений индексов */
            // int leftAreaBound = ;
            // int RightAreaBound = ;

            /* Этот цикл тупой он по всем элементам идет */
            for (int j = 0; j < Dim; j++)
            {
                if (IsInArea(Bound, Grid[j], j))
                {
                    /* Мы в заданной области устанавливаем нужные параметры */
                    InfoManeger::SetFictitious(Grid[j].info, 1);
                    InfoManeger::SetAreaInfo(Grid[j].info, Bound.AreaNum);
                }
            }
        }

        /* Очистка локальной переменной  */
        delete[] BoundGrid;
    }
    catch (std::bad_alloc &e)
    {
        status.SetStatus(State::MEMORY_ALLOC_ERROR, "Ошибка при выделении памяти в Grid2D_Quad::GenerateBaseGrid(GridStatus &status)\n");
        if (BoundGrid != nullptr)
            delete[] BoundGrid;
        return;
    }
    catch (const std::exception &e)
    {
        status.SetStatus(State::UNKNOWN_ERROR, e.what());
        if (BoundGrid != nullptr)
            delete[] BoundGrid;
        return;
    }
}

void Grid2D_Quad::DivisionIntoSubBounds(GridStatus &status) noexcept
{
    if (status.GetState() != State::OK)
        return;

    struct Bound
    {
        int Size = 0;
        int *Line = nullptr; // Массив точек границы
        int BoundType = 0;
        int BoundFormula = -1;
    };

    /* Бинарный поиск по массиву */
    /*
        @param
        const BoundArea& Bound - Массив границ
        int numPointGlobal - Значение для поиска
        ret bool: true - элемент в массиве есть false в противном случае

    */
    auto BinarySerch = [](const Bound &bound, int numPointGlobal) -> bool
    {
        int *arr = bound.Line;
        int left = 0;
        int right = bound.Size - 1;
        int midd = 0;
        while (1)
        {
            midd = (left + right) / 2;

            if (numPointGlobal < arr[midd])      // если искомое меньше значения в ячейке
                right = midd - 1;                // смещаем правую границу поиска
            else if (numPointGlobal > arr[midd]) // если искомое больше значения в ячейке
                left = midd + 1;                 // смещаем левую границу поиска
            else                                 // иначе (значения равны)
                break;                           // функция возвращает индекс ячейки

            // если границы сомкнулись
            if (left > right)
            {
                midd = -1;
                break;
            }
        }
        if (midd == -1)
            return false;
        else
            return true;
    };

    auto IsInBound = [&](const Bound &bound, int numPointGlobal)
    { return BinarySerch(bound, numPointGlobal); };

    int *MemoryPool = nullptr;

    try
    {

        MemoryPool = new int[max(GlobalNx, GlobalNy)];
        /* Краевые условия задаются практически так же как и области с тем лишь исключением, что одна из координат фиксируется */
        auto GetBound = [&](int i) -> Bound
        {
            Bound bound;
            bound.BoundType = baseGrid.BoundsArea[i][1];
            bound.BoundFormula = baseGrid.BoundsArea[i][0];

            /* х - фиксирован */
            if (baseGrid.BoundsArea[i][2] == baseGrid.BoundsArea[i][3])
            {
                /* Определим базовые узлы на прямой OY а потом растиражируем узлы по границе ( Дробление фиксированное )  шаг по массиву GlobalNx*/
                int StartPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1); // Стартовая позиция для  Y
                int EndPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][5], 1);   // Конечная точка для оси Y

                int nY = (EndPositionY - StartPositionY) / GlobalNx; //  Количество узлов по оси Y

                int Idx = StartPositionY;
                for (int i = 0; i <= nY; i++)
                {
                    MemoryPool[i] = Idx;
                    Idx += GlobalNx;
                }
                bound.Line = MemoryPool;
                bound.Size = nY;
            }
            /* y - фиксирован */
            else if (baseGrid.BoundsArea[i][4] == baseGrid.BoundsArea[i][5])
            {
                /* Определяем базовые узлы по оси OX, а потом аналогично растиражируем узлы по границе (Дробление фиксированное) шаг по массиву + 1 */
                int StartPositionX = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1); // Стартовая позиция по оси Х
                int EndPositionX = Getlevel(baseGrid.BoundsArea[i][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1);   // Конечная позиция по оси Х
                int nX = EndPositionX - StartPositionX;                                                                          //  Количество узлов по оси X

                int idxBoundGrid = 0;
                int Idx = StartPositionX;
                for (int j = 0; j <= nX; j++)
                {
                    MemoryPool[j] = Idx;
                    Idx++;
                }

                bound.Line = MemoryPool;
                bound.Size = nX;
            }

            return bound;
        };

        for (int i = 0; i < baseGrid.P; i++)
        {
            /* Получить список точек границы  */
            Bound bound = GetBound(i);
            for (int j = 0; j < Dim; j++)
            {
                /* Расставляем нужные значения в соответствующие точки сетки */
                if (IsInBound(bound, j))
                {
                    /* Мы в заданной области устанавливаем нужные параметры */
                    InfoManeger::SetBoundInfo(Grid[j].info, bound.BoundFormula, bound.BoundType + 1);
                }
            }
        }
    }
    catch (std::bad_alloc &e)
    {
        status.SetStatus(State::MEMORY_ALLOC_ERROR, "Ошибка при выделении памяти в Grid2D_Quad::GenerateBaseGrid(GridStatus &status)\n");
        if (MemoryPool != nullptr)
            delete[] MemoryPool;
        return;
    }
    catch (const std::exception &e)
    {
        status.SetStatus(State::UNKNOWN_ERROR, e.what());
        if (MemoryPool != nullptr)
            delete[] MemoryPool;
        return;
    }
}

/*******************************************************/
/* Public section */
Grid2D_Quad::Grid2D_Quad(const string &filename)
{
    GridStatus status = Load(filename);
    if (status.GetState() != State::OK)
        throw status;
}

Grid2D_Quad::Grid2D_Quad(const BaseGrid2D &baseGrid_) : baseGrid(baseGrid_) {}

GridStatus Grid2D_Quad::Load(const string &filename) noexcept
{
    GridStatus status;

    fin.open(filename);
    if (!fin.is_open())
    {
        status.SetStatus(State::FILE_OPEN_ERROR, "Ошибка открытия файла " + filename + " в Grid2D_Quad::Grid2D_Quad(const string &filename)\n");
        fout.close();
        return status;
    }

    /* Файл корректно открылся и можно читать
        Корректность входных данных не проверяется
    */
    /* Базовая сетка по XY */
    fin >> baseGrid.Nx >> baseGrid.Ny;
    baseGrid.BaseGridXY = vector(baseGrid.Ny, vector<BaseGrid2D::PointXY>(baseGrid.Nx));
    for (int i = 0; i < baseGrid.Ny; i++)
    {
        for (int j = 0; j < baseGrid.Nx; j++)
        {
            fin >> baseGrid.BaseGridXY[i][j].x >> baseGrid.BaseGridXY[i][j].y;
        }
    }

    /* Расчетные подобласти */
    fin >> baseGrid.L;
    baseGrid.CalculationArea = vector(baseGrid.L, vector<int>(baseGrid.SizeOfCalculationAreaElemet));

    for (int i = 0; i < baseGrid.L; i++)
    {
        for (int j = 0; j < baseGrid.SizeOfCalculationAreaElemet; j++)
        {
            fin >> baseGrid.CalculationArea[i][j];
            baseGrid.CalculationArea[i][j]--; // Приведение нумерации с нуля
        }
    }
    /***********************/

    /* Описание Границ */
    fin >> baseGrid.P;
    baseGrid.BoundsArea = vector(baseGrid.P, vector<int>(baseGrid.SizeOfBoundsAreaElement));
    for (int i = 0; i < baseGrid.P; i++)
    {
        for (int j = 0; j < baseGrid.SizeOfBoundsAreaElement; j++)
        {
            fin >> baseGrid.BoundsArea[i][j];
            baseGrid.BoundsArea[i][j]--; // Приведение к нумерации с нуля
        }
    }
    /***********************/

    /* Правила дробления базовой сетки */
    baseGrid.DivideParam.resize(baseGrid.SizeOfDivideParam);
    baseGrid.DivideParam[0].resize(baseGrid.Nx - 1);
    baseGrid.DivideParam[1].resize(baseGrid.Ny - 1);

    for (int i = 0; i < baseGrid.Nx - 1; i++)
    {
        fin >> baseGrid.DivideParam[0][i].num >> baseGrid.DivideParam[0][i].coef;
    }

    for (int i = 0; i < baseGrid.Ny - 1; i++)
    {
        fin >> baseGrid.DivideParam[1][i].num >> baseGrid.DivideParam[1][i].coef;
    }

    baseGrid.isReadyToUse = true;

    return status;
}

GridStatus Grid2D_Quad::GenerateGrid() noexcept
{
    GridStatus status;

    // Расчет общего количества узлов в сетке 
	// Для этого пройдемся по массиву разбиения каждого отрезка и вычислим общее число узлов
    GetTotalNumberOfNodes();
    
    // Генерация базовой сетки 
    GenerateBaseGrid(status);
    
    // Учет фиктивных узлов 
    DivisionIntoSubAreas(status);

    // Учет КУ и расстановка границ 
    DivisionIntoSubBounds(status);

    return status;
}

GridStatus Grid2D_Quad::DivideGrid(const int coef) noexcept
{
    GridStatus status;
    return status;
}

GridStatus Grid2D_Quad::ReGenerateGrid() noexcept
{
    GridStatus status;
    return status;
}

/* Гетеры сеттеры */
BaseGrid2D Grid2D_Quad::GetBaseGrid() const noexcept { return baseGrid; }
Cross_TypeTemplate Grid2D_Quad::GetElement(int idx) const noexcept
{
    Cross_TypeTemplate Element;
    return Element;
}
GridStatus Grid2D_Quad::SetBaseGrid(const BaseGrid2D &baseGrid_) noexcept
{
    GridStatus status;
    return status;
}


Point& Grid2D_Quad::operator[](int idx) noexcept { return Grid[idx]; }

/* Debug functions */
void Grid2D_Quad::PrintBaseGrid() const noexcept
{
    cout << "Base Grid Print\n";
    cout << "Nx = " << baseGrid.Nx << " Ny = " << baseGrid.Ny << "\n";
    for (int i = 0; i < baseGrid.Ny; i++)
    {
        cout << i + 1 << ": ";
        for (int j = 0; j < baseGrid.Nx; j++)
        {
            cout << "(" << baseGrid.BaseGridXY[i][j].x << ";" << baseGrid.BaseGridXY[i][j].y << ") ";
        }
        cout << "\n";
    }
    cout << "\n L = " << baseGrid.L << "\n";
    for (int i = 0; i < baseGrid.L; i++)
    {
        cout << i + 1 << ": ";
        for (int j = 0; j < baseGrid.SizeOfCalculationAreaElemet; j++)
        {
            cout << baseGrid.CalculationArea[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n P = " << baseGrid.P << "\n";
    for (int i = 0; i < baseGrid.P; i++)
    {
        cout << i + 1 << ": ";
        for (int j = 0; j < baseGrid.SizeOfBoundsAreaElement; j++)
        {
            cout << baseGrid.BoundsArea[i][j] << " ";
        }
        cout << "\n";
    }

    cout << "\nDivide Parametrs\n";
    cout << "\nDivide X\n";
    for (int i = 0; i < baseGrid.Nx - 1; i++)
    {
        cout << "(num = " << baseGrid.DivideParam[0][i].num << "; coef = " << baseGrid.DivideParam[0][i].coef << ")\n";
    }
    cout << "\nDivide Y\n";
    for (int i = 0; i < baseGrid.Ny - 1; i++)
    {
        cout << "(num = " << baseGrid.DivideParam[1][i].num << "; coef = " << baseGrid.DivideParam[1][i].coef << ")\n";
    }
}

void Grid2D_Quad::PrintGrid() const noexcept
{
    int idx = 0;
    std::cout << "Start idx: " << idx << "  End idx: " << GlobalNx * GlobalNy - 1 << " Step Row: " << GlobalNx << "\n";

    for (int i = 0; i < GlobalNy; i++)
    {
        for (int j = 0; j < GlobalNx; j++)
        {
            std::cout << std::fixed << std::setprecision(2) << "(" << Grid[idx].x << ";" << Grid[idx].y << ") ";
            idx++;
        }
        std::cout << "\n";
    }
}

/*******************************************************/