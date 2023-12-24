#ifndef MATRIX_H_
#define MATRIX_H_

#include "../status/matrixStatus.h"

#include <memory>
#include <string>
#include <iostream>

namespace mathtools::matrix
{

    using namespace std;
    using namespace status;

    /* Базовый интерфейс матриц */
    class Matrix
    {
    protected:
        StatusMatrix status;
        size_t N = 0; // Размер
    public:
        Matrix() {}

        /* Гетеры сетеры */
        virtual StatusMatrix GetElement(const size_t i, const size_t j, double &elem) noexcept = 0;
        virtual StatusMatrix SetElement(const size_t i, const size_t j, const double elem) noexcept = 0;

        /* Операции проверки состояний и прочее */
        StatusMatrix GetStatus() const noexcept;
        bool IsValidMatrix() const noexcept;
        size_t GetSize() const noexcept;

        /* Умножение матрицы на вектор */
        virtual StatusMatrix MultVec(const unique_ptr<double[]> &x, unique_ptr<double[]> &res, const size_t size) noexcept = 0;
        virtual StatusMatrix MultVec(const double *x, double *res, const size_t size) noexcept = 0;

        virtual ~Matrix() {}
    };

    /***********************************************************************************************************************************/
    /* Плотный формат хранения матрицы  */

    class DenseMatrix : public Matrix
    {
    private:
        unique_ptr<unique_ptr<double[]>[]> matr; // Матрица

    public:
        DenseMatrix(const size_t N);
        DenseMatrix(const DenseMatrix &m);
        DenseMatrix(DenseMatrix &&m);

        /* Гетеры сетеры */
        virtual StatusMatrix GetElement(const size_t i, const size_t j, double &elem)  noexcept final;
        virtual StatusMatrix SetElement(const size_t i, const size_t j, const double elem) noexcept final;


        /* Умножение матрицы на вектор */
        virtual StatusMatrix MultVec(const unique_ptr<double[]> &x, unique_ptr<double[]> &res, const size_t size) noexcept final;
        virtual StatusMatrix MultVec(const double *x, double *res, const size_t size) noexcept final;

        /* Удаленные операции */
        DenseMatrix &operator=(const DenseMatrix &) = delete;
        DenseMatrix &operator=(DenseMatrix &&) = delete;
        DenseMatrix &operator=(DenseMatrix &) = delete;

        /* Дружественные функции */
        friend ostream& operator<<(ostream& cnt, const DenseMatrix &matrix);

        virtual ~DenseMatrix() {}
    };

    /***********************************************************************************************************************************/

    /* Разряженный строчно столбцовый формат хранения квадртных матриц */
    class SparseMatrix: public Matrix
    {
        private:
            
        public:

    };

};

#endif