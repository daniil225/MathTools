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
    /* S2 */
    idx = 0;
    start = Grid.GlobalNx*(Grid.GlobalNy-1);
    end = Grid.GlobalNx*Grid.GlobalNy;
    cout << end << "\n";
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
    N--; // При расчетах один узел учитывался 2 раза пожтому от общего числа отнимаем 1
    cout << N << "\n";
    Grid.PrintGrid();
}

void FDM::Solve()
{
    //uij.clear();

}

 double FDM::Norma() const
 {
    double res = 0;

    return res;
 }