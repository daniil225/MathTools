#pragma once

#include "DataStruct.h"
#include "ISolverSLAU.h"
#include "Factorization.h"



namespace SolverSLAU {

	//using namespace Factorization;

	//  Модификации алгоритмов для решений слау  
	//	NONE - Без модификаций стандартный алгоритм 
	//	DIAGONALFACTORIZATION - Диагональное предобуславливание 
	//	LLTFACTORIZATION - Неполная факторизация холецкого 
	//	LUSQFACTORIZATION - Неполная факторизация с применением алгориитма разложения матрицы LU(sq)
	enum class MSGResolverModification
	{
		NONE,
		DIAGONALFACTORIZATION,
		LLTFACTORIZATION,
		//LUSQFACTORIZATION
	};


	/****************************************************************************************************************************/
	/******************************************************* MSG *********************************************************************/

	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	class MSG : public ISolverSLAU<OutData, Matrix, SLAU>
	{
	private:

		MSGResolverModification modification;

		// Определим ряд вспомогательных векторов для реализации методов 
		//MSGResolverModification::NONE
		std::vector<double> rk;
		std::vector<double> Azk;
		std::vector<double> zk;

		/* Для факторизации LLT */
		Factorization::SPARSE_ROW_COL_SYMETRIC_FACT FactMatr;

		// xk - это вектор в выходных данных 
		double f_norm = 0.0; // Норма вектора f
		double rk_norm = 0.0; // Норма вектора невязки
		double non_repan = 0.0; // Невязка 

		// Инициализацию векторов будем проиизводить непосредственно перед решением задачи 
		void MSG_Classic(OutData& x);

		// Диагональное предобуславливание для рещения сисетмы 
		void MSG_Diagonal(OutData& x);

		void MSG_LLTFact(OutData& x);


	public:

		MSG(SLAU<Matrix>& data, MSGResolverModification modification_ = MSGResolverModification::NONE) :
			ISolverSLAU<OutData, Matrix, SLAU>(data), modification(modification_) {}


		// Реализация интерфейса Вызываем функцию решения СЛАУ с заданными параметрами 
		void solve(OutData& out);

		~MSG() {}

	};


	/* Private section  */
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void  MSG<OutData, Matrix, SLAU>::MSG_Classic(OutData& outData)
	{
		// Реализуем классический алгоритм Метода сопряженных градиентов 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double eps = this->Data.eps;
		// дадим ссылки на все элементы данных 
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = outData.x;
		f_norm = std::sqrt(f * f); // Норма вектора f

		if (this->Data.x.size() == 0)
		{
			// Если начального приближения нет то 
			xk.resize(N);
		}
		else
		{
			//  Если начальное приближение все же есть, то
			xk.resize(N);
			std::copy(this->Data.x.begin(), this->Data.x.end(), xk.begin());
		}

		rk.resize(N);
		Azk.resize(N);
		zk.resize(N);

		// Если начального приближения нет 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Azk = this->MultA(xk); // A*x0
			rk = f - Azk; // r0 = f - A*x0
		}

		std::copy(rk.begin(), rk.end(), zk.begin()); // z0 = r0
		non_repan = std::sqrt(rk * rk) / f_norm; // Относительная невязка 

		for (int k = 1; k < maxiter; k++)
		{
			Azk = this->MultA(zk);

			double dot_rk = rk * rk;
			double alphak = dot_rk / (Azk * zk);
			xk = xk + alphak * zk;
			rk = rk - alphak * Azk;

			double bettak = rk * rk / dot_rk;

			zk = rk + bettak * zk;

			non_repan = std::sqrt(rk * rk) / f_norm;
			std::cout << "Iteration: " << k << " Non-repan: " << non_repan << "\n";

			if (non_repan < eps) break;
		}

		// Пересчитываем истинную невязку 
		rk = f - this->MultA(xk);
		non_repan = std::sqrt(rk * rk) / f_norm;

		std::cout << "Истинная невязка r = f - Ax: " << non_repan << "\n";

	}


	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void  MSG<OutData, Matrix, SLAU>::MSG_Diagonal(OutData& outData)
	{
		// Реализуем алгоритм Метода сопряженных градиентов с диагональным предобуславливанием 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double eps = this->Data.eps;
		// дадим ссылки на все элементы данных 
		std::vector<double>& di = this->Data.matr.di;
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = outData.x;
		f_norm = std::sqrt(f * f); // Норма вектора f

		if (this->Data.x.size() == 0)
		{
			// Если начального приближения нет то 
			xk.resize(N);
		}
		else
		{
			//  Если начальное приближение все же есть, то
			xk.resize(N);
			std::copy(this->Data.x.begin(), this->Data.x.end(), xk.begin());
		}

		rk.resize(N);
		Azk.resize(N);
		zk.resize(N);
		// Если начального приближения нет 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Azk = this->MultA(xk); // A*x0
			rk = f - Azk; // r0 = f - A*x0
		}

		// Хочу определить лямбда функцию для специфического умножения вектора на вектор 
		auto MultD = [](const std::vector<double> &di,const std::vector<double>& vec)
		{
			int size = vec.size();
			std::vector<double> res(size);

			for (int i = 0; i < size; i++)
				res[i] = vec[i] / di[i];
			return res;
		};

		zk = MultD(di,rk);

		for (int k = 1; k < maxiter; k++)
		{
			Azk = this->MultA(zk);
			double dot_rk = MultD(di, rk) * rk;
			double alphak = dot_rk /( Azk * zk);

			xk = xk + alphak * zk;
			rk = rk - alphak * Azk;

			std::vector<double> Mrk = MultD(di, rk);
			double bettak = (Mrk * rk) / dot_rk;
			zk = Mrk + bettak * zk;


			non_repan = std::sqrt(rk * rk) / f_norm;
			std::cout << "Iteration: " << k << " Non-repan: " << non_repan << "\n";

			if (non_repan < eps) break;
		}

		// Пересчитываем истинную невязку 
		rk = f - this->MultA(xk);
		non_repan = std::sqrt(rk * rk) / f_norm;

		std::cout << "Истинная невязка r = f - Ax: " << non_repan << "\n";
		
	}

	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void  MSG<OutData, Matrix, SLAU>::MSG_LLTFact(OutData& outData)
	{
		// Реализуем алгоритм Метода сопряженных градиентов с диагональным предобуславливанием 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double eps = this->Data.eps;
		// дадим ссылки на все элементы данных 
		std::vector<double>& di = this->Data.matr.di;
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = outData.x;
		f_norm = std::sqrt(f * f); // Норма вектора f

		if (this->Data.x.size() == 0)
		{
			// Если начального приближения нет то 
			xk.resize(N);
		}
		else
		{
			//  Если начальное приближение все же есть, то
			xk.resize(N);
			std::copy(this->Data.x.begin(), this->Data.x.end(), xk.begin());
		}

		rk.resize(N);
		Azk.resize(N);
		zk.resize(N);
		// Если начального приближения нет 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Azk = this->MultA(xk); // A*x0
			rk = f - Azk; // r0 = f - A*x0
		}

		
		// Неполное разложение 
		FactMatr = Factorization::LLTFact(this->Data.matr);
		

		zk = this->M_1x(FactMatr, rk);


		for (int k = 1; k < maxiter; k++)
		{
			Azk = this->MultA(zk);
			double dot_rk = this->M_1x(FactMatr, rk) * rk;
			double alphak = dot_rk / (Azk * zk);

			xk = xk + alphak * zk;
			rk = rk - alphak * Azk;

			std::vector<double> tmp = this->M_1x(FactMatr, rk);

			double bettak = (tmp * rk) / dot_rk;

			zk = tmp + bettak * zk;


			non_repan = std::sqrt(rk * rk) / f_norm;
			std::cout << "Iteration: " << k << " Non-repan: " << non_repan << "\n";

			if (non_repan < eps) break;
		}


		// Пересчитываем истинную невязку 
		rk = f - this->MultA(xk);
		non_repan = std::sqrt(rk * rk) / f_norm;

		std::cout << "Истинная невязка r = f - Ax: " << non_repan << "\n";

	}

	/* Public section */
	// Реализация интерфейса Вызываем функцию решения СЛАУ с заданными параметрами 
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void MSG<OutData, Matrix, SLAU>::solve(OutData& out)
	{
		switch (modification)
		{
		case MSGResolverModification::NONE:
			MSG_Classic(out);
			break;
		case MSGResolverModification::DIAGONALFACTORIZATION:
			MSG_Diagonal(out);
			break;
		case MSGResolverModification::LLTFACTORIZATION:
			MSG_LLTFact(out);
			break;
		default:
			std::cout << "Неопределенный тип модификации\n";
			break;
		}
	}
};