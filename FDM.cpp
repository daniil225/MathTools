#include "FDM.h"
#include <cmath>
#include <iostream>
#include <stdio.h>

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
    // Размеры массивов под граничные массивы 
    BoundsNodes.resize(4*max(Grid.GlobalNx, Grid.GlobalNy));

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
                        if(typeCond == 1 || (int)bif.TypeCond[1] == 1)
                        {
                            // В рамках задачи на границе может быть только 2 и 1 или 3 и 1 => т.к 1КУ перекрывают все заносим индекс в массив 
                            BoundsNodes[BoundsNodesIdx] = m;
                            BoundsNodesIdx++;
                            //cout << "\nNode num = " << m << " Bound Info\n";
                            //cout << "Node = " << m << "\n";
                        }
                        else if(typeCond == 2)
                        {
                            
                            hy = Grid[m + Grid.GlobalNx].y - Grid[m].y;
                            // Учет 2 КУ везде это первый элемент на границе S6 
                            matr.add(m, m, lambda/hy);
                            matr.add(m, m+Grid.GlobalNx, -lambda/hy);
                            slau.f[m] = deq.Bound[(int)bif.Cond[0]].func[0](Grid[m].x, Grid[m].y);
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
        //matr.PrintDenseFormatMatrix();
        Zeidel(matr, uij, slau);
        //int k = 0;
        // for(double u: uij)
        // {
        //     if(k == Grid.GlobalNx)
        //     {
        //         cout << "\n";
        //         k = 0;
        //     }
        //     k++;
        //     cout << u << " ";
        // }

    }
    else
    {

    }
    
}

 double FDM::Norma()
 {
    double res = 0;

    /* В цикле по узлам */
    for(int i = 0; i < Grid.GlobalNy; i++)
    {
        for(int j = 0; j < Grid.GlobalNx; j++)
        {
            int m = i*Grid.GlobalNx + j;
            
            if(InfoManeger::IsFiFictitious(Grid[m].info) && !InfoManeger::IsBound(Grid[m].info))
            {
                res += pow(deq.u_true(Grid[m].x, Grid[m].y) - uij[m], 2);
            }
        }
    }

    return sqrt(res);
 }

 void FDM::PrintTable()
 {
     /* В цикле по узлам */
    printf("\n\n Расчетная таблица. Символом * отмечены внутренние узлы сетки\n");
    printf("-------------------------------------------------------------------------------------------\n");
    printf("|N             |X             |Y             |U             |U*            |   |U* - U|   |\n");
    printf("|--------------|--------------|--------------|--------------|--------------|--------------|\n");
    for(int i = 0; i < Grid.GlobalNy; i++)
    {
        for(int j = 0; j < Grid.GlobalNx; j++)
        {
            int m = i*Grid.GlobalNx + j;
            if(InfoManeger::IsFiFictitious(Grid[m].info))
            {
                double u_true = deq.u_true(Grid[m].x, Grid[m].y);
                if(InfoManeger::IsBound(Grid[m].info))
                {
                    printf("|%14d|%14.3f|%14.3f| %e | %e | %e |\n", m,Grid[m].x, Grid[m].y, uij[m], u_true, abs(u_true - uij[m]));
                }
                else
                {
                    printf("|*%13d|%14.3f|%14.3f| %e | %e | %e |\n", m,Grid[m].x, Grid[m].y, uij[m], u_true, abs(u_true - uij[m]));
                }
               
               
            }   
        }
    }
    printf("-------------------------------------------------------------------------------------------\n");
 }