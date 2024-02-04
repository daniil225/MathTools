#include "FDM.h"
#include <cmath>
#include <iostream>

/* Private section */
bool FDM::GenerateMatrix()
{
    return true;
}

/*********************************************************************************/

/* Public section */
FDM::FDM(const string &filename, DEquation deq_, bool IsravGrid_): deq(deq_), IsravGrid(IsravGrid_)
{
    Grid.Load(filename);
    Grid.GenerateGrid();

    /* Определим сумарное количество не фиктивных узлов - это будет соответсвовать размернсти матрицы 
       А так же определим угловые узлы, которые не входят в расчетную область 
    */
    
    /* Выставляем Границу вместе с нумерацией  всего 6 границ с разными КУ*/
    int start;
    int end;
    int size;
    int idx = 0;
    
    /* S1 */
    /* 0;0 не включается  и верхную тоже не вкулючаем */
    start = Grid.GlobalNx;
    end = (Grid.GlobalNy-1)*Grid.GlobalNx;
    size = (end-start)/Grid.GlobalNx;

    deq.Bound[0].idx.resize(size);
    for(int i = start; i < end; i+= Grid.GlobalNx)
    {
        deq.Bound[0].idx[idx] = i;
        N++;
        idx++;
    }
    /***************************************************/

    /* S2 */
    idx = 0;
    start = Grid.GlobalNx*(Grid.GlobalNy-1);
    end = Grid.GlobalNx*Grid.GlobalNy;
    deq.Bound[1].idx.resize(Grid.GlobalNx);
    for(int i = start+1; i < end; i++)
    {
        if(InfoManeger::IsFiFictitious(Grid[i+1].info))
        {
            deq.Bound[1].idx[idx] = i;
            idx++;
            N++;
        }
        else break;
    }
    deq.Bound[1].idx.resize(idx);
    /***************************************************/

    /* S3 */
    start = Grid.GlobalNx*(Grid.GlobalNy-1) + idx + 1 - Grid.GlobalNx;
    end = start % Grid.GlobalNx;
    idx = 0;
    deq.Bound[2].idx.resize(Grid.GlobalNy);
    for(int i = start; i >= end; i -= Grid.GlobalNx)
    {
        if( InfoManeger::IsBound(Grid[i].info))
        {
            deq.Bound[2].idx[idx] = i;
            idx++;
            N++;
        }
        else break;
    }
    deq.Bound[2].idx.resize(idx);
    /***************************************************/

    /* S4 */
    idx = 0;
    start = deq.Bound[2].idx[deq.Bound[2].idx.size()-1];
    end = start + (Grid.GlobalNx -start%Grid.GlobalNx - 2); // -2 так как точка с краю не включается в расчетную область 
    size = end - start + 1;
    deq.Bound[3].idx.resize(size);

    for(int i = start; i <= end; i++)
    {
        deq.Bound[3].idx[idx] = i;
        idx ++;
        N++;
    }
    /***************************************************/

    /* S5 */
    idx = 0;
    start = deq.Bound[3].idx[deq.Bound[3].idx.size()-1] + 1 - Grid.GlobalNx;
    end = Grid.GlobalNx - 1;
    deq.Bound[4].idx.resize(Grid.GlobalNy);
    for(int i = start; i > end; i -= Grid.GlobalNx)
    {
        deq.Bound[4].idx[idx] = i;
        idx++;
        N++;
    }
    deq.Bound[4].idx.resize(idx);
    /***************************************************/

    /* S6 */
    start = 1;
    end = Grid.GlobalNx - 2;
    size = end - start + 1;
    deq.Bound[5].idx.resize(size);
    idx = 0;
    for(int i = start; i <= end; i++)
    {
        deq.Bound[5].idx[idx] = i;
        idx++;
        N++;
    }
    /***************************************************/


    /* Внутренние узлы в расчетной области */
    InternalNode.resize(Grid.Dim);
    idx = 0;
    for(int i = 0; i < Grid.Dim; i++)
    {
        Point point = Grid[i];
        if(InfoManeger::IsFiFictitious(point.info) && !InfoManeger::IsBound(point.info))
        {
            InternalNode[idx] = i;
            idx++;
            N++;
        }
    }

    InternalNode.resize(idx);
    N--; // При расчетах один узел учитывался 2 раза поэтому от общего числа отнимаем 1

    // Расчет m
    m = Grid.GlobalNx - 2;
    
    // Сохранение памяти под матрицу
    matr.m = m;
    matr.N = Grid.Dim;
    matr.SaveMemory(); // Есть матрица уже нужной структуры 
    
    CheckNode.resize(Grid.Dim);

    //Grid.PrintGrid();
}

void FDM::Solve()
{


    uij.clear();
    // Очистка того что было, если было 
    uij.resize(Grid.Dim);
    slau.f.resize(Grid.Dim);

    // Собираем матрицу c учетом равномерности сетки  
    if(IsravGrid)
    {
        cout << "Work with rav Grid\n";
        // Внутренние узлы сначала 
        
        double hx = Grid[1].x - Grid[0].x;
        double hy = Grid[Grid.GlobalNx].y - Grid[0].y;

        for(int k = 0; k < InternalNode.size(); k++)
        {
            // Выведем индексы для точки шаблона крест а после конвертируем эти числа в натации (i;j) индексов 
            int kij = InternalNode[k]; // Узел с которого стартуем 
            int kim1j = kij-1;
            int kip1j = kij+1;
            int kijm1 = kij - Grid.GlobalNx;
            int kijp1 = kij + Grid.GlobalNx;
            
            CheckNode[kij] = true; // Уравнение на этот узел есть все тип топ 
            //cout << "kij = " << kij << " kip1j = " << kip1j << " kim1j = " << kim1j << " kijmi = " << kijm1 << " kijp1 = " << kijp1 << "\n";
            // u(i-1, j)/hx2
            matr.add(kim1j, kij, -1.0/(hx*hx));
            // u(i+1, j)/hx2
            matr.add(kip1j, kij, -1.0/(hx*hx));
            // -2u(i, j)/hx2 - 2u(i,j)/hy2
            matr.add(kij, kij, 2.0/(hx*hx) +2.0/(hy*hy) );
            
            // u(i, j-1)/hy2
            matr.add(kij, kijm1, -1.0/(hy*hy));
            
            // u(i, j+1)/hy2
            matr.add(kij, kijp1, -1.0/(hy*hy));
        
            // Учет вектора f
        }
        
        // Учет 3-х КУ

        // Учет 1-х КУ
        // Все Ку 1 ого типа значит на строке удаляем все значения кроме диагональной 
        // S1
        if(deq.Bound[0].typeBound == 1)
        {
            // Подставляем значения по оси Y
            vector<int> &idx = deq.Bound[0].idx;
            function<double(double)> phi = deq.Bound[0].func[0];
            for(int k = 0; k < idx.size(); k++)
            {
                int kij = idx[k]; // Диагональный элемент
                double y = Grid[kij].y;
                matr.insert(kij, kij, 1.0); // На Диагональ ставим 1
                matr.insert(kij, kij-1, 0.0);
                matr.insert(kij, kij+1, 0.0);
                matr.insert(kij, kij+2+m, 0.0);
                matr.insert(kij, kij-2-m, 0.0);

                CheckNode[kij] = true; // Уравнение на этот узел есть все тип топ 
                slau.f[kij] = phi(y);
            }
        }

        if(deq.Bound[1].typeBound == 1)
        {
            // Подставляем значения по оси Y
            vector<int> &idx = deq.Bound[1].idx;
            function<double(double)> phi = deq.Bound[1].func[0];
            for(int k = 0; k < idx.size(); k++)
            {
                int kij = idx[k]; // Диагональный элемент
                double x = Grid[kij].x;
                matr.insert(kij, kij, 1.0); // На Диагональ ставим 1
                matr.insert(kij, kij-1, 0.0);
                matr.insert(kij, kij+1, 0.0);
                matr.insert(kij, kij+2+m, 0.0);
                matr.insert(kij, kij-2-m, 0.0);
                CheckNode[kij] = true; // Уравнение на этот узел есть все тип топ 

                slau.f[kij] = phi(x);
            }
        }
        
        if(deq.Bound[2].typeBound == 1)
        {
            // Подставляем значения по оси Y
            vector<int> &idx = deq.Bound[2].idx;
            function<double(double)> phi = deq.Bound[2].func[0];
            for(int k = 0; k < idx.size(); k++)
            {
                int kij = idx[k]; // Диагональный элемент
                double y = Grid[kij].y;
                matr.insert(kij, kij, 1.0); // На Диагональ ставим 1
                matr.insert(kij, kij-1, 0.0);
                matr.insert(kij, kij+1, 0.0);
                matr.insert(kij, kij+2+m, 0.0);
                matr.insert(kij, kij-2-m, 0.0);

                CheckNode[kij] = true; // Уравнение на этот узел есть все тип топ 
                slau.f[kij] = phi(y);
            }
        }
        
        if(deq.Bound[3].typeBound == 1)
        {
            // Подставляем значения по оси Y
            vector<int> &idx = deq.Bound[3].idx;
            function<double(double)> phi = deq.Bound[3].func[0];
            for(int k = 0; k < idx.size(); k++)
            {
                int kij = idx[k]; // Диагональный элемент
                double x = Grid[kij].x;
                matr.insert(kij, kij, 1.0); // На Диагональ ставим 1
                matr.insert(kij, kij-1, 0.0);
                matr.insert(kij, kij+1, 0.0);
                matr.insert(kij, kij+2+m, 0.0);
                matr.insert(kij, kij-2-m, 0.0);

                CheckNode[kij] = true; // Уравнение на этот узел есть все тип топ 
                slau.f[kij] = phi(x);
            }
        }
        
        if(deq.Bound[4].typeBound == 1)
        {
            // Подставляем значения по оси Y
            vector<int> &idx = deq.Bound[4].idx;
            function<double(double)> phi = deq.Bound[4].func[0];
            for(int k = 0; k < idx.size(); k++)
            {
                int kij = idx[k]; // Диагональный элемент
                double y = Grid[kij].y;
                matr.insert(kij, kij, 1.0); // На Диагональ ставим 1
                matr.insert(kij, kij-1, 0.0);
                matr.insert(kij, kij+1, 0.0);
                matr.insert(kij, kij+2+m, 0.0);
                matr.insert(kij, kij-2-m, 0.0);
                CheckNode[kij] = true; // Уравнение на этот узел есть все тип топ 

                slau.f[kij] = phi(y);
            }
        }
        
        if(deq.Bound[5].typeBound == 1)
        {
            // Подставляем значения по оси Y
            vector<int> &idx = deq.Bound[5].idx;
            function<double(double)> phi = deq.Bound[5].func[0];
            for(int k = 0; k < idx.size(); k++)
            {
                int kij = idx[k]; // Диагональный элемент
                double x = Grid[kij].x;
                matr.insert(kij, kij, 1.0); // На Диагональ ставим 1
                matr.insert(kij, kij-1, 0.0);
                matr.insert(kij, kij+1, 0.0);
                matr.insert(kij, kij+2+m, 0.0);
                matr.insert(kij, kij-2-m, 0.0);
                CheckNode[kij] = true; // Уравнение на этот узел есть все тип топ 

                slau.f[kij] = phi(x);
            }
            
        }

        /* Делаем СЛАУ разрешимой */
        for(int i = 0; i < CheckNode.size(); i++)
        {
            if(!CheckNode[i])
            {
                matr.insert(i,i, 1.0);
            }
        }

        matr.PrintDenseFormatMatrix();

        Zeidel(matr, uij, slau);

        int idx = 0;
        for(auto u: uij)
        {
            if(idx == Grid.GlobalNx)
            {
                cout << "\n";
                idx = 0;
            }
            cout << u << " ";
            idx++;
        }

    }
    else
    {

    }
    
}

 double FDM::Norma() const
 {
    double res = 0;

    return res;
 }