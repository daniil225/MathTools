
#include <iostream>
#include "FDM.h"
#include <cmath>


int main()
{
    DEquation deq;
    BoundCond bound1;
   
    bound1.typeBound = 1;
    bound1.func[0] = [](double x,double y) -> double {return 2*x*x + 3*y*y; };

    BoundCond bound2;
    
    bound2.typeBound = 1;
    bound2.func[0] = [](double x,double y) -> double {return 2*x*x + 3*y*y; };

    BoundCond bound3;
    
    bound3.typeBound = 1;
    bound3.func[0] = [](double x,double y) -> double {return 2*x*x + 3*y*y; };

    BoundCond bound4;
    
    bound4.typeBound = 1;
    bound4.func[0] = [](double x,double y) -> double {return 2*x*x + 3*y*y; };
    
    BoundCond bound5;
   
    bound5.typeBound = 1;
    bound5.func[0] = [](double x,double y) -> double {return 2*x*x + 3*y*y; };

    BoundCond bound6;
    
    bound6.typeBound = 1;
    bound6.func[0] = [](double x,double y) -> double {return 2*x*x + 3*y*y; };

    BoundCond bound7;
    
    bound7.typeBound = 1;
    bound7.func[0] = [](double x,double y) -> double {return 2*x*x + 3*y*y; };

    BoundCond bound8;
    
    bound8.typeBound = 1;
    bound8.func[0] = [](double x,double y) -> double {return 2*x*x + 3*y*y; };

    // BoundCond bound9;
    
    // bound9.typeBound = 1;
    // bound9.func[0] = [](double x,double y) -> double {return x*x*x + y*y*y; };

    // BoundCond bound10;
    
    // bound10.typeBound = 1;
    // bound10.func[0] = [](double x,double y) -> double {return x*x*x + y*y*y; };


    deq.f = [](double x, double y) -> double { return -5.0;};
    //deq.gamma = 1;
    deq.lambda = 1;
    deq.u_true = [](double x, double y) -> double { return 2*x*x + 3*y*y; };



    deq.Bound.resize(8);
    deq.Bound[0] = bound1;
    deq.Bound[1] = bound2;
    deq.Bound[2] = bound3;
    deq.Bound[3] = bound4;
    deq.Bound[4] = bound5;
    deq.Bound[5] = bound6;
    deq.Bound[6] = bound7;
    deq.Bound[7] = bound8;
    //deq.Bound[8] = bound9;
    // deq.Bound[9] = bound10;
    FDM fdm("CalcArea1.txt", deq, true);
    fdm.Solve();
    double norma = fdm.Norma();

    cout << "\n\n Норма разности || U* - U || = " << norma << "\n";
    fdm.PrintTable();
    
    return 0;
}