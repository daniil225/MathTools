#include "matrix.h"



namespace mathtools::matrix
{
    
        StatusMatrix Matrix::GetStatus() const noexcept { return status; }
        bool Matrix::IsValidMatrix() const noexcept { return status.GetStatusMatrixCode() == StatusMatrixCode::Ok ? true: false; }
        size_t Matrix::GetSize() const noexcept { return N; }

};







namespace mathtools::matrix
{

    DenseMatrix::DenseMatrix(size_t N)
    {
        try
        {
            this->N = N;
            matr = std::unique_ptr<std::unique_ptr<double[]>[]>(new std::unique_ptr<double[]>[N]);

            for (size_t i = 0; i < N; i++)
                matr[i] = make_unique<double[]>(N);

            /* Все хорошо корректное создание */
            this->status.SetStatus(StatusMatrixCode::Ok, "Create success");
        }
        catch (bad_alloc &e)
        {
            cerr << "Memory error\n";
            cerr << e.what();
            exit(-1); // Тормозим если что то не так 
        }
        catch (...)
        {
            this->status.SetStatus(StatusMatrixCode::UnknownError, "Unknown Error from  DenseMatrix::DenseMatrix(const DenseMatrix &m)");
            exit(-1); // Тормозим если что то не так 
        }
    }

    DenseMatrix::DenseMatrix(const DenseMatrix &m)
    {

        try
        {
            this->N = m.N;
            matr = std::unique_ptr<std::unique_ptr<double[]>[]>(new std::unique_ptr<double[]>[N]);

            /* Копируем матрицу */
            for (size_t i = 0; i < N; i++)
                matr[i] = make_unique<double[]>(N);

            for (size_t i = 0; i < N; i++)
            {
                for (size_t j = 0; j < N; j++)
                {
                    matr[i][j] = m.matr[i][j];
                }
            }
        }
        catch (bad_alloc &e)
        {
            cerr << "Memory error\n";
            cerr << e.what();
            exit(-1); // Тормозим если что то не так 
        }
        catch (...)
        {
            this->status.SetStatus(StatusMatrixCode::UnknownError, "Unknown Error from  DenseMatrix::DenseMatrix(const DenseMatrix &m)");
            exit(-1); // Тормозим если что то не так 
        }
    }

    DenseMatrix::DenseMatrix(DenseMatrix &&m)
    {
        try
        {
            this->N = m.N;
            /* Перенос значения матрица более будт не валидна */
            matr.swap(m.matr);
            
            /* Обнуляем матрицу */
            m.N = 0;
            m.status.SetStatus(StatusMatrixCode::InvalideMatrix, "Matrix was copied");
        }
        catch (...)
        {
            this->status.SetStatus(StatusMatrixCode::UnknownError, "Unknown Error from  DenseMatrix::DenseMatrix(const DenseMatrix &m)");
            exit(-1); // Тормозим если что то не так 
        }
    }
    
    StatusMatrix DenseMatrix::GetElement(const size_t i, const size_t j, double &elem) noexcept
    {

        if(i < N && j < N)
        {
            elem = matr[i][j];
            status.SetStatus(StatusMatrixCode::Ok, "Success");
        }
        else
        {
            status.SetStatus(StatusMatrixCode::OutOfRange, "Acess error in DenseMatrix::GetElement()");
        }

        return status;
    }
    
    StatusMatrix DenseMatrix::SetElement(const size_t i, const size_t j, const double elem) noexcept
    {
        if(i < N && j < N)
        {
            matr[i][j] = elem;
            status.SetStatus(StatusMatrixCode::Ok, "Success");
        }
        else
        {
            status.SetStatus(StatusMatrixCode::OutOfRange, "Set error in DenseMatrix::SetElement()");
        }
        return status;
    }
    

    StatusMatrix DenseMatrix::MultVec(const unique_ptr<double[]> &x, unique_ptr<double[]> &res, const size_t size) noexcept
    {
        if(size <= N)
        {
            for(size_t i = 0; i < N; i++)
            {
                double s = 0.0;
                for(size_t j = 0; j < N; j++)
                {
                    s += matr[i][j]*x[i];
                }
                res[i] = s;
            }
            status.SetStatus(StatusMatrixCode::Ok, "Success");
        }
        else
        {
            status.SetStatus(StatusMatrixCode::OutOfRange, "Set error in DenseMatrix::MultVec()");
        }

        return status;
    }   
    
    StatusMatrix DenseMatrix::MultVec(const double *x, double *res, const size_t size) noexcept 
    {
        if(size <= N)
        {
            for(size_t i = 0; i < N; i++)
            {
                double s = 0.0;
                for(size_t j = 0; j < N; j++)
                {
                    s += matr[i][j]*x[i];
                }
                res[i] = s;
            }
            status.SetStatus(StatusMatrixCode::Ok, "Success");
        }
        else
        {
             status.SetStatus(StatusMatrixCode::OutOfRange, "Set error in DenseMatrix::MultVec()");
        }

        return status;
    }


    ostream& operator<<(ostream& cnt, const DenseMatrix &matrix)
    {

        cnt << "Class Dense Matrix\n";
        cnt << "Size: " << matrix.N << "*" << matrix.N << "\n";
        cnt << "Data: \n";

        for(size_t i = 0; i < matrix.N; i++)
        {
            for(size_t j = 0; j < matrix.N; j++)
            {
                cnt << matrix.matr[i][j] << " ";
            }
            cnt << "\n";
        }
        cnt << "\n";
        return cnt;
    }
};