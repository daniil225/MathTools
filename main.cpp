#include "PointInfo.h"

int main()
{

    Info info;

    InfoManeger::ClearInfo(info);
    InfoManeger::SetFictitious(info, 1);
    
    InfoManeger::SetAreaInfo(info, 12);
    InfoManeger::SetAreaInfo(info, 11);
    InfoManeger::SetAreaInfo(info, 10);
    InfoManeger::SetAreaInfo(info, 10);
    InfoManeger::SetAreaInfo(info, 10);
    InfoManeger::SetAreaInfo(info, 10);

    InfoManeger::SetBoundInfo(info, 1, 3);
    InfoManeger::SetBoundInfo(info, 2, 1);
    InfoManeger::SetBoundInfo(info, 2, 1);
    InfoManeger::SetBoundInfo(info, 2, 1);

    InfoManeger::PrintInfo(info);

    BoundInfo binf = InfoManeger::GetBoundInfo(info);
    InfoManeger::PrintBoundInfo(binf);

    return 0;
}