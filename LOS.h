#pragma once
#include"ISolverSLAU.h"
#include "Factorization.h"
#include "DataStruct.h"

namespace SolverSLAU {
	/****************************************************************************************************************************/
	/******************************************************* LOS ***************************************************************/
	using namespace Factorization;

	enum class LOSResolverModification
	{
		NONE,
		DIAGONALFACTORIZATION,
		LLTFACTORIZATION,
		//LUSQFACTORIZATION
	};


	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	class LOS : public ISolverSLAU<OutData, Matrix, SLAU>
	{
	private:

		/* Вспомогаткльные структуры данных */
		std::vector<double> rk;
		std::vector<double> Ark;
		std::vector<double> zk;
		std::vector<double> pk;

		/* Диагональ для факторизации диагональной */
		std::vector<double> D;
		std::vector<double> SAQrk;
		/* Для предобуславливания неполным разложением Холецкого */
		SPARSE_ROW_COL_SYMETRIC_FACT FactMatr;
		std::vector<double> Qrk; 


		// Определим ряд вспомогательных векторов для реализации методов 
		void LOS_Classic(OutData& out);

		void LOS_Diagonal(OutData& out);

		void LOS_LLTFact(OutData& out);

		LOSResolverModification modification;
	public:

		LOS(SLAU<Matrix>& data, LOSResolverModification modification_ = LOSResolverModification::NONE) : 
			ISolverSLAU<OutData, Matrix, SLAU>(data), modification(modification_) {}

		// Реализация интерфейса Вызываем функцию решения СЛАУ с заданными параметрами 
		void solve(OutData& out);

	};


	/* Private section  */
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void LOS<OutData, Matrix, SLAU>::LOS_Classic(OutData& out)
	{
		// Реализуем классический алгоритм Метода сопряженных градиентов 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double eps = this->Data.eps;
		// дадим ссылки на все элементы данных 
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = out.x;

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
		Ark.resize(N);
		zk.resize(N);
		pk.resize(N);


		// Если начального приближения нет 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Ark = this->MultA(xk); // A*x0
			rk = f - Ark; // r0 = f - A*x0
		}

		std::copy(rk.begin(), rk.end(), zk.begin()); // z0 = r0

		pk = this->MultA(rk); // p0 = Ar0
		std::copy(pk.begin(), pk.end(), Ark.begin());

		double dot_rk = rk * rk; //  Невязка через скалярное произведение 
		for (int k = 1; k < maxiter; k++)
		{
			double dot_pk = pk * pk;
			double alphak = (pk * rk) / dot_pk;
			// xk
			xk = xk + alphak * zk;
			// rk
			rk = rk - alphak * pk;
			// Bettak
			Ark = this->MultA(rk);
			double bettak = -1.0 * (pk * Ark) / dot_pk;
			zk = rk + bettak * zk;
			pk = Ark + bettak * pk;

			// (rk; rk)
			dot_rk -= alphak * alphak * dot_pk;
			std::cout << "Iteration: " << k << " (rk; rk) = " << dot_rk << "\n";
			if (dot_rk < eps) break;
		}
		// Пересчитываем истинную невязку 
		rk = f - this->MultA(xk);
		dot_rk = rk * rk;

		std::cout << "Истинная невязка r = f - Ax: " << dot_rk << "\n";


	}
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void LOS<OutData, Matrix, SLAU>::LOS_Diagonal(OutData& out)
	{
		// Реализуем классический алгоритм Метода сопряженных градиентов 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double r0_norm2 = 0.0;
		double eps = this->Data.eps* this->Data.eps;
		// дадим ссылки на все элементы данных 
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = out.x;
		std::vector<double>& di = this->Data.matr.di;

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
		Ark.resize(N);
		zk.resize(N);
		pk.resize(N);

		// Матрица диагонального предобуславливания 
		D = DiagonalFact(di);
		// Хочу определить лямбда функцию для специфического умножения вектора на вектор 
		auto MultD = [](const std::vector<double>& di, const std::vector<double>& vec)
		{
			int size = vec.size();
			std::vector<double> res(size);

			for (int i = 0; i < size; i++)
				res[i] = vec[i] / di[i];
			return res;
		};

		// Если начального приближения нет 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Ark = this->MultA(xk); // A*x0
			rk = f - Ark; // r0 = f - A*x0
		}

		rk = MultD(D, rk); // r0 с волночкой
		r0_norm2 = rk * rk;

		zk = MultD(D, rk);

		pk = MultD(D,(this->MultA(zk)));
		SAQrk.resize(N);
		
		double dot_rk = rk * rk; //  Невязка через скалярное произведение 
		for (int k = 1; k < maxiter; k++)
		{
			double dot_pk = pk * pk;
			double alphak = (pk * rk) / dot_pk;
			// xk
			xk = xk + alphak * zk;
			// rk
			rk = rk - alphak * pk;
			// Bettak
			SAQrk = MultD(D, this->MultA( MultD(D, rk) ) );
			double bettak = -1.0 * (pk*SAQrk) / dot_pk;
			zk = MultD(D, rk) + bettak * zk;
			pk = SAQrk + bettak * pk;

			// (rk; rk)
			dot_rk -= alphak * alphak * dot_pk;
			std::cout << "Iteration: " << k << " (rk; rk) = " << dot_rk/r0_norm2 << "\n";
			if (dot_rk/r0_norm2 < eps) break;
		}
		
		// Пересчитываем истинную невязку 
		rk = f - this->MultA(xk);
		dot_rk = rk * rk;

		std::cout << "Истинная невязка r = f - Ax: " << dot_rk/r0_norm2 << "\n";

	}

	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void LOS<OutData, Matrix, SLAU>::LOS_LLTFact(OutData& out)
	{
		// Реализуем классический алгоритм Метода сопряженных градиентов 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double r0_norm2 = 0.0;
		double eps = this->Data.eps * this->Data.eps;
		// дадим ссылки на все элементы данных 
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = out.x;

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
		Ark.resize(N);
		zk.resize(N);
		pk.resize(N);

		// Если начального приближения нет 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Ark = this->MultA(xk); // A*x0
			rk = f - Ark; // r0 = f - A*x0
		}

		// Факторизация неполным разложением Холецкого
		FactMatr = LLTFact(this->Data.matr);

		/* Измененная часть специфичная для этой функции */
		rk = this->normal(FactMatr, rk);

		zk = this->reverse(FactMatr, rk);

		pk = this->normal(FactMatr,this->MultA(zk));

		r0_norm2 = rk * rk;
		double dot_rk = rk * rk; //  Невязка через скалярное произведение 
		SAQrk.resize(N);
		Qrk.resize(N);

		for (int k = 1; k < maxiter; k++)
		{
			double dot_pk = pk * pk;
			double alphak = (pk * rk) / dot_pk;
			// xk
			xk = xk + alphak * zk;
			// rk
			rk = rk - alphak * pk;
			// Bettak
			Qrk = this->reverse(FactMatr, rk);
			SAQrk = this->normal(FactMatr,this->MultA(Qrk));

			double bettak = -1.0 * (pk * SAQrk) / dot_pk;
			zk = Qrk + bettak * zk;
			pk = SAQrk + bettak * pk;

			// (rk; rk)
			dot_rk -= alphak * alphak * dot_pk;
			std::cout << "Iteration: " << k << " (rk; rk) = " << dot_rk/r0_norm2 << " < " << eps << "\n";
			if (dot_rk / r0_norm2 < eps) break;
		}

		// Пересчитываем истинную невязку 
		rk = f - this->MultA(xk);
		dot_rk = rk * rk;

		std::cout << "Истинная невязка r = f - Ax: " << dot_rk / r0_norm2 << "\n";
	}

	/* Public section */

	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void LOS<OutData, Matrix, SLAU>::solve(OutData& out)
	{
		switch (modification)
		{
		case LOSResolverModification::NONE:
			LOS_Classic(out);
			break;
		case LOSResolverModification::DIAGONALFACTORIZATION:
			LOS_Diagonal(out);
			break;
		case LOSResolverModification::LLTFACTORIZATION:
			LOS_LLTFact(out);
			break;
		default:
			std::cout << "Неопределенный тип модификации\n";
			break;
		}
	}

};