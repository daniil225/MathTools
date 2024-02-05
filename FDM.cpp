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
    
    // Расчет m
    m = Grid.GlobalNx - 2;
    
    // Сохранение памяти под матрицу
    matr.m = m;
    matr.N = Grid.Dim;
    matr.SaveMemory(); // Есть матрица уже нужной структуры 
    matr.PrintDenseFormatMatrix();
    // Размеры массивов под граничные массивы 
    BoundsNodes.resize(4*max(Grid.GlobalNx, Grid.GlobalNy));

    Grid.PrintGrid();
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
        double hx, hy;
        double lambda = deq.lambda;
        double gamma = deq.gamma;
        for(int i = 0; i < Grid.GlobalNy; i++)
        {
            for(int j = 0; j < Grid.GlobalNx; j++)
            {
                int m = i*Grid.GlobalNx + j;
                Point centralNode = Grid[m];
                if(InfoManeger::IsFiFictitious(centralNode.info))
                {
                    // Узел не фиктивный 
                    if(InfoManeger::IsBound(centralNode.info))
                    {
                        // Это граница
                        BoundInfo bif = InfoManeger::GetBoundInfo(centralNode.info);
                        int typeCond = (int32_t)bif.TypeCond[0];
                        // 1 КУ
                        if(typeCond == 1)
                        {
                            BoundsNodes[BoundsNodesIdx] = m;
                            BoundsNodesIdx++;
                        }
                        else if(typeCond == 2)
                        {
                            // 2КУ Заносим в матрицу нужные знаения 
                        }
                        else
                        {
                            // 3КУ Заносим в матрицу нужные значения 
                        
                        }
                        
                        //InfoManeger::PrintInfo(centralNode.info);
                        //cout << "Node = " << m << " Formula nums = " << (int32_t)bif.Cond[0] << " ";
                        //cout << "Type Cond = " <<  << "\n"; 
                    }
                    else
                    {
                        // Внутренний узел 
                        int kim1j = m-1;
                        int kip1j = m + 1;
                        int kijm1 = m - Grid.GlobalNx;
                        int kijp1 = m + Grid.GlobalNx;
                        
                        hx = Grid[kip1j].x - Grid[m].x;
                        hy = Grid[kijp1].y - Grid[m].y;
                        
                        // Строка фиксирована это m ая строка 
                        matr.add(m, kim1j, -lambda/(hx*hx));
                        matr.add(m, kip1j, -lambda/(hx*hx));
                        matr.add(m, kijm1, -lambda/(hy*hy));
                        matr.add(m, kijp1, -lambda/(hy*hy));
                        matr.add(m, m, lambda*( (2)/(hx*hx) + (2)/(hy*hy) ) + gamma);

                        // Вектор f
                        slau.f[m] += deq.f(Grid[m].x, Grid[m].y);
                    }
                }
                else
                {
                    // Фиктивный узел 
                    matr.insert(m,m,1.0);
                }
            }
        }

        /* А теперь учитываем 1 - ое краевое условие */
        for(int i = 0; i < BoundsNodesIdx; i++)
        {
            int m = BoundsNodes[i];
            matr.insert(m, m, 1.0);
            BoundInfo bif = InfoManeger::GetBoundInfo(Grid[m].info);
            // Вносим значения в праую часть 
            slau.f[m] = deq.Bound[(int)bif.Cond[0]].func[0](Grid[m].x, Grid[m].y);
        }

        Zeidel(matr, uij, slau);
        int k = 0;
        for(double u: uij)
        {
            if(k == Grid.GlobalNx)
            {
                cout << "\n";
                k = 0;
            }
            k++;
            cout << u << " ";
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