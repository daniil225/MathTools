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
FDM::FDM(const string &filename, DEquation deq_, bool IsravGrid_) : deq(deq_), IsravGrid(IsravGrid_)
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
    BoundsNodes.resize(4 * max(Grid.GlobalNx, Grid.GlobalNy));

    //Grid.PrintGrid();
}

void FDM::Solve()
{

    uij.clear();
    // Очистка того что было, если было
    uij.resize(Grid.Dim);
    slau.f.resize(Grid.Dim);

    // Собираем матрицу для произвольной сетки
    if (IsravGrid)
    {
        double hxRight, hxLeft, hyTop, hyBottom;
        double lambda = deq.lambda;
        double gamma = deq.gamma;
        for (int i = 0; i < Grid.GlobalNy; i++)
        {
            for (int j = 0; j < Grid.GlobalNx; j++)
            {
                int m = i * Grid.GlobalNx + j;
                Point centralNode = Grid[m];
                if (InfoManeger::IsFiFictitious(centralNode.info))
                {
                    // Узел не фиктивный
                    if (InfoManeger::IsBound(centralNode.info))
                    {
                        // Это граница
                        BoundInfo bif = InfoManeger::GetBoundInfo(centralNode.info);
                        int typeCond = (int32_t)bif.TypeCond[0];
                        // 1 КУ
                        if (typeCond == 1 || (int)bif.TypeCond[1] == 1)
                        {
                            // В рамках задачи на границе может быть только 2 и 1 или 3 и 1 => т.к 1КУ перекрывают все заносим индекс в массив
                            BoundsNodes[BoundsNodesIdx] = m;
                            BoundsNodesIdx++;
                            // cout << "\n----------------------------------\n";
                            // cout << "\nNode num = " << m << " Bound Info\n";
                            // InfoManeger::PrintBoundInfo(bif);
                            // cout << "---------------------------------\n";
                        }
                        else if (typeCond == 2)
                        {
                            //  Граница S6
                            if (int(m / Grid.GlobalNx) == 0)
                            {
                                hyTop = Grid[m + Grid.GlobalNx].y - Grid[m].y;
                                // Учет 2 КУ везде это первый элемент на границе S6
                                matr.add(m, m, lambda / hyTop);
                                matr.add(m, m + Grid.GlobalNx, -lambda / hyTop);
                                slau.f[m] += deq.Bound[(int)bif.Cond[0]].func[0](Grid[m].x, Grid[m].y);
                                //cout << "\nNode num(S6) = " << m << " Bound Info\n";
                            }
                            // Граница S3
                            else
                            {
                                hxLeft = Grid[m].x - Grid[m-1].x;
                                matr.add(m, m, lambda/hxLeft);
                                matr.add(m, m-1, -lambda/hxLeft);
                                slau.f[m] += deq.Bound[(int)bif.Cond[0]].func[0](Grid[m].x, Grid[m].y);
                                //cout << "\nNode num(S3) = " << m << " Bound Info\n";
                            }
                            
                        }
                        else
                        {
                            // 3КУ Заносим в матрицу нужные значения
                        }

                        // InfoManeger::PrintInfo(centralNode.info);
                        // cout << "Node = " << m << " Formula nums = " << (int32_t)bif.Cond[0] << " ";
                        // cout << "Type Cond = " <<  << "\n";
                    }
                    else
                    {
                        // Внутренний узел
                        int kim1j = m - 1;
                        int kip1j = m + 1;
                        int kijm1 = m - Grid.GlobalNx;
                        int kijp1 = m + Grid.GlobalNx;

                        hxRight = Grid[kip1j].x - Grid[m].x;
                        hxLeft = Grid[m].x - Grid[kim1j].x;
                        hyTop = Grid[kijp1].y - Grid[m].y;
                        hyBottom = Grid[m].y - Grid[kijm1].y;
                        // Строка фиксирована это m ая строка

                        matr.add(m, kim1j, -(2 * lambda) / (hxLeft * (hxLeft + hxRight)));
                        matr.add(m, kip1j, -(2 * lambda) / (hxRight * (hxLeft + hxRight)));
                        matr.add(m, kijm1, -(2 * lambda) / (hyBottom * (hyBottom + hyTop)));
                        matr.add(m, kijp1, -(2 * lambda) / (hyTop * (hyBottom + hyTop)));
                        matr.add(m, m, lambda * ((2.0) / (hxRight * hxLeft) + (2.0) / (hyTop * hyBottom)) + gamma);

                        // Вектор f
                        slau.f[m] += deq.f(Grid[m].x, Grid[m].y);
                    }
                }
                else
                {
                    // Фиктивный узел
                    slau.f[m] = 0;
                    matr.insert(m, m, 1.0);
                }
            }
        }

        /* А теперь учитываем 1 - ое краевое условие */
        for (int i = 0; i < BoundsNodesIdx; i++)
        {
            int m = BoundsNodes[i];
            matr.insert(m, m, 1.0);
            BoundInfo bif = InfoManeger::GetBoundInfo(Grid[m].info);
            // Вносим значения в праую часть
            // cout << "Node Num = " << m << " x = " << Grid[m].x << " y = " << Grid[m].y << " val = " << deq.Bound[0].func[0](Grid[m].x, Grid[m].y) << "\n";

            // Берем от нулевой потому что это всегда 1-КУ
            if((int)bif.TypeCond[0] == 1)
            {
                slau.f[m] += deq.Bound[(int)bif.Cond[0]].func[0](Grid[m].x, Grid[m].y);
            }
            else
            {
                cout << "Node m = " << m << "\n";
                slau.f[m] += deq.Bound[(int)bif.Cond[1]].func[0](Grid[m].x, Grid[m].y);
            }
            
        }
        //matr.PrintDenseFormatMatrix();
        //Grid.PrintGrid();
        Zeidel(matr, uij, slau, 1.0);
        // int k = 0;
        //  for(double u: uij)
        //  {
        //      if(k == Grid.GlobalNx)
        //      {
        //          cout << "\n";
        //          k = 0;
        //      }
        //      k++;
        //      cout << u << " ";
        //  }
    }
    else
    {
    }
}

double FDM::Norma()
{
    double res = 0;

    /* В цикле по узлам */
    for (int i = 0; i < Grid.GlobalNy; i += coef)
    {
        for (int j = 0; j < Grid.GlobalNx; j += coef)
        {
            int m = i * Grid.GlobalNx + j;

            if (InfoManeger::IsFiFictitious(Grid[m].info))
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
    printf("-------------------------------------------------------------------------------------------------\n");
    printf("|N             |X             |Y             |U               |U*              |     |U* - U|   |\n");
    printf("|--------------|--------------|--------------|----------------|----------------|----------------|\n");
    for (int i = 0; i < Grid.GlobalNy; i+=coef)
    {
        for (int j = 0; j < Grid.GlobalNx; j+=coef)
        {
            int m = i * Grid.GlobalNx + j;
            if (InfoManeger::IsFiFictitious(Grid[m].info))
            {
                double u_true = deq.u_true(Grid[m].x, Grid[m].y);
                if (InfoManeger::IsBound(Grid[m].info))
                {
                    printf("|%14d|%14.3f|%14.3f| %14e | %14e | %14e |\n", m, Grid[m].x, Grid[m].y, uij[m], u_true, abs(u_true - uij[m]));
                }
                else
                {
                    printf("|*%13d|%14.3f|%14.3f| %14e | %14e | %14e |\n", m, Grid[m].x, Grid[m].y, uij[m], u_true, abs(u_true - uij[m]));
                }
            }
        }
        if (i < Grid.GlobalNy - 1)
            printf("|--------------|--------------|--------------|----------------|----------------|----------------|\n");
    }
    printf("-------------------------------------------------------------------------------------------------\n\n\n");
}