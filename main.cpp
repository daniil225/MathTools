
#include <iostream>
#include "FDM.h"



int main()
{
    DEquation deq;
    BoundCond bound1;
   
    bound1.typeBound = 1;
    bound1.func[0] = [](double x,double y) -> double { return x*x*x+y*y*y; };

    BoundCond bound2;
    
    bound2.typeBound = 1;
    bound2.func[0] = [](double x,double y) -> double { return x*x*x+y*y*y; };

    BoundCond bound3;
    
    bound3.typeBound = 1;
    bound3.func[0] = [](double x,double y) -> double { return x*x*x+y*y*y; };

    BoundCond bound4;
    
    bound4.typeBound = 1;
    bound4.func[0] = [](double x,double y) -> double { return x*x*x+y*y*y; };
    
    BoundCond bound5;
   
    bound5.typeBound = 1;
    bound5.func[0] = [](double x,double y) -> double { return x*x*x+y*y*y; };

    BoundCond bound6;
    
    bound6.typeBound = 1;
    bound6.func[0] = [](double x,double y) -> double { return x*x*x+y*y*y; };

    deq.f = [](double x, double y) -> double { return -3.0*x - 3*y; };
    deq.gamma = 0;
    deq.lambda = 1;
    deq.Bound.resize(6);
    deq.Bound[0] = bound1;
    deq.Bound[1] = bound2;
    deq.Bound[2] = bound3;
    deq.Bound[3] = bound4;
    deq.Bound[4] = bound5;
    deq.Bound[5] = bound6;
    

    


    FDM fdm("CalcArea1.txt", deq, true);
    fdm.Solve();
    return 0;
}