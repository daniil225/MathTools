#include "Grid2D_Quad.h"
#include <iostream>
int main()
{

    Grid2D_Quad Grid("CalcArea1.txt");
    GridStatus status = Grid.GenerateGrid();

    //

    cout << status.GetMsg();
    Grid.PrintGrid();
    InfoManeger::PrintInfo(Grid[0].info);
    InfoManeger::PrintInfo(Grid[1].info);
    InfoManeger::PrintInfo(Grid[2].info);
    InfoManeger::PrintInfo(Grid[3].info);
    InfoManeger::PrintInfo(Grid[4].info);
    InfoManeger::PrintInfo(Grid[5].info);
    InfoManeger::PrintInfo(Grid[6].info);
    return 0;
}