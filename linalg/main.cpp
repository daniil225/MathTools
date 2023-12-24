#include "matrix.h"

using namespace mathtools::matrix;

int main()
{
    DenseMatrix matrix(10);
    DenseMatrix matrix2(move(matrix));
    
    double x[10] = {1,1,1,1,1,1,1,1,1,1};

    double res[10];

    double elem;
    StatusMatrix stat;
    for(size_t i = 0; i < matrix2.GetSize(); i++)
    {
        for(size_t j = 0; j < matrix2.GetSize(); j++)
        {
            stat = matrix2.SetElement(i, j, (i+j)%5);
            if( stat.GetStatusMatrixCode() != StatusMatrixCode::Ok)
            {
                stat.what();
                exit(-1);
            }
        }
    }

    cout << matrix2;


    return 0;
}