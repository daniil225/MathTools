#include "Grid2D_Quad.h"
#include <iostream>
int main()
{

    Grid2D_Quad Grid("CalcArea1.txt");
    GridStatus status = Grid.GenerateGrid();


    // cout << status.GetMsg();
    // Grid.PrintGrid();

    for(int i = 20; i <=24; i++)
    {
        cout << "Num = " << i << "\n";
        InfoManeger::PrintInfo(Grid[i].info);
        cout << "-------------------------------------------\n";
    }
    return 0;
}