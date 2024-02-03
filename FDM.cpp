#include "FDM.h"
#include <cmath>

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