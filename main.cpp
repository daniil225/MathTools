
#include <iostream>
#include "FDM.h"
#include <cmath>


int main()
{
    DEquation deq;
    BoundCond bound1;
   
    bound1.typeBound = 1;
    bound1.func[0] = [](double x,double y) -> double { return 2*x+4*y; };

    BoundCond bound2;
    
    bound2.typeBound = 1;
    bound2.func[0] = [](double x,double y) -> double { return 2*x+4*y; };

    BoundCond bound3;
    
    bound3.typeBound = 1;
    bound3.func[0] = [](double x,double y) -> double { return 2*x+4*y; };

    BoundCond bound4;
    
    bound4.typeBound = 1;
    bound4.func[0] = [](double x,double y) -> double { return 2*x+4*y; };
    
    BoundCond bound5;
   
    bound5.typeBound = 1;
    bound5.func[0] = [](double x,double y) -> double { return 2*x+4*y; };

    BoundCond bound6;
    
    bound6.typeBound = 2;
    bound6.func[0] = [](double x,double y) -> double { return -8.0; };

    deq.f = [](double x, double y) -> double { return 6*x + 12*y; };
    deq.gamma = 3;
    deq.lambda = 2;
    deq.u_true = [](double x, double y) -> double { return 2*x + 4*y; };

    deq.Bound.resize(6);
    deq.Bound[0] = bound1;
    deq.Bound[1] = bound2;
    deq.Bound[2] = bound3;
    deq.Bound[3] = bound4;
    deq.Bound[4] = bound5;
    deq.Bound[5] = bound6;
    

    

    


    FDM fdm("CalcArea1.txt", deq, true);
    fdm.Solve();
    double norma = fdm.Norma();

    cout << "\n\n Норма разности || U* - U || = " << norma << "\n";
    fdm.PrintTable();
    
    return 0;
}