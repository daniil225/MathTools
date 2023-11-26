#pragma once

#include "config.h"
#include "Support.h"
#include <cmath>

namespace SolverSLAU
{

	/****************************************** INTERFACE  *****************************************************/

	// Определелим интерфейс решателя через него будем наследоваться к различным решателям 
	template<typename OutData,typename Matrix ,template<typename Matrix> class SLAU>
	class ISolverSLAU
	{
	protected:
		// Инициализируем Данные для СЛАУ здесь => производный класс уже будет иметь все необходимые данные 
		
		SLAU<Matrix> Data;

		// Умножение матрицы на вектор Для всех типов матриц функция одна и та же 
		#ifdef _SPARSE_
		
		#ifdef _SPARSE_ROW_COL_SYMETRIC_
		std::vector<double> MultA(const std::vector<double>& x);

		/* Прямой и обратный ход Гауса для нахождения произведения обратной матрицы на заданный вектор */

		template<typename FactMatr>
		std::vector<double> normal(FactMatr& matr, std::vector<double> b)
		{
			int N = Data.matr.di.size();
			std::vector<int>& ig = Data.matr.ig;
			std::vector<int>& jg = Data.matr.jg;
			std::vector<double>& gg_l = matr.gg;
			std::vector<double>& di_l = matr.di;

			for (int i = 0; i < N; i++) {
				for (int j = ig[i]; j < ig[i + 1]; j++)
					b[i] -= gg_l[j] * b[jg[j]];

				b[i] = b[i] / di_l[i];
			}

			return b;
		}

		template <typename FactMatr>
		std::vector<double> reverse(FactMatr& matr, std::vector<double> x) 
		{
			int N = Data.matr.di.size();
			std::vector<int>& ig = Data.matr.ig;
			std::vector<int>& jg = Data.matr.jg;
			std::vector<double>& gg_l = matr.gg;
			std::vector<double>& di_l = matr.di;

			for (int j = N-1; j >= 0; j--) {
				x[j] = x[j] / di_l[j];

				for (int i = ig[j]; i < ig[j + 1]; i++)
					x[jg[i]] -= gg_l[i] * x[j];
			}
			return x;
		}

		/* Операция умножения M^-1*x */
		template <typename FactMatr>
		std::vector<double> M_1x(FactMatr& matr, std::vector<double> x) { return reverse(matr, normal(matr, x)); }

		#endif

		#ifdef _SPARSE_ROW_COL_
		std::vector<double> MultA(const std::vector<double>& x);
		#endif
		
		#endif
	
	public:
		ISolverSLAU(){}

		// Передача параметра по ссылке
		explicit ISolverSLAU(SLAU<Matrix>& data): Data(data) { }
		
		// Перемещающее присваивание
		template<typename Matrix>
		explicit ISolverSLAU(SLAU<Matrix>&& data): Data(std::move(data)) {}

		// Решение СЛАУ
		// Выходные данные решателя 
		virtual void solve(OutData& out) = 0;

		virtual ~ISolverSLAU() {}
	};

	#ifdef _SPARSE_

	#ifdef _SPARSE_ROW_COL_SYMETRIC_
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	std::vector<double> ISolverSLAU<OutData, Matrix, SLAU>::MultA(const std::vector<double>& x)
	{
		std::vector<double> res;
		std::vector<int>& ig = Data.matr.ig;
		std::vector<int>& jg = Data.matr.jg;
		std::vector<double>& gg = Data.matr.gg;
		std::vector<double>& di = Data.matr.di;
		int N = Data.f.size();
		res.resize(N);
		
		for (int i = 0; i < N; i++)
		{
			res[i] = di[i] * x[i];
			for (int j = ig[i]; j < ig[i + 1]; j++)
			{
				int k = jg[j];
				res[i] += gg[j] * x[k];
				res[k] += gg[j] * x[i];
			}
		}
		return res;
	}
	#endif

	#ifdef _SPARSE_ROW_COL_
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	std::vector<double> ISolverSLAU<OutData, Matrix, SLAU>::MultA(const std::vector<double>& x)
	{
		std::vector<double> res;
		std::vector<int>& ig = Data.matr.ig;
		std::vector<int>& jg = Data.matr.jg;
		std::vector<double>& ggu = Data.matr.ggu;
		std::vector<double>& ggl = Data.matr.ggl;
		std::vector<double>& di = Data.matr.di;
		int N = Data.f.size();
		res.resize(N);

		for (int i = 0; i < N; i++)
		{
			res[i] = di[i] * x[i];
			for (int j = ig[i]; j < ig[i + 1]; j++)
			{
				int k = jg[j];
				res[i] += ggl[j] * x[k];
				res[k] += ggu[j] * x[i];
			}
		}
		return res;
	}
	#endif
	#endif



	/****************************************************************************************************************************/

};