#include "Grid_1D.h"

int main()
{
    Grid_1D grid("CalcArea.txt");
    grid.GenerateGrid();
    grid.PrintGrid();

    return 0;
}