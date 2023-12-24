#include "matrixStatus.h"

namespace status
{
    void StatusMatrix::what()
    {
        if (code == StatusMatrixCode::Ok)
        {
            cout << "\nStatus: OK\n";
            cout << msg << "\n";
        }
        else if (code == StatusMatrixCode::InvalideMatrix)
        {
            cerr << "\nStatus code: InvalideMatrix\n";
            cerr << msg << "\n";
        }
        else if(code == StatusMatrixCode::UnknownError)
        {
            cerr << "\nStatus code: error\n";
            cerr << msg << "\n";
        }
        else if(code == StatusMatrixCode::OutOfRange)
        {
            cerr << "\nStatus code: out of range error\n";
            cerr << msg << "\n";
        }
    }

    void StatusMatrix::SetStatus(StatusMatrixCode status, const string& m)
    {
        this->code = status;
        this->msg = m;
    }

    StatusMatrixCode StatusMatrix::GetStatusMatrixCode() const { return code; }
};