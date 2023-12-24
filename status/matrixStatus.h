#include <exception>
#include <iostream>

namespace status
{
    using namespace std;
    

    enum class StatusMatrixCode
    {
        Ok,
        InvalideMatrix,
        UnknownError,
        OutOfRange
    };

    class [[nodiscard]]  StatusMatrix: public exception
    {
        private:
            StatusMatrixCode code;
            string msg;
        public:
            StatusMatrix() {}
            StatusMatrix(StatusMatrixCode N, const string& m):code(N), msg(move(m)) {}

            StatusMatrixCode GetStatusMatrixCode() const;
            void what();
            void SetStatus(StatusMatrixCode status, const string& m);

            ~StatusMatrix() {}
    };
}